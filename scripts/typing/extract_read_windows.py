#!/usr/bin/env python3
"""Per-read signal windows around labelled sites, for per-read modification typing.

For every (read, site) pair where the read covers the site ON THE STRAND THAT
CARRIES THE MODIFICATION, cut a window of 2*half_window+1 bases centred on the
site, in pore frame (the read's own 5'->3' order; see strand_orientation_check:
a reverse-mapped read's signal follows the complement). Each base's signal
segment (Dorado move table) is resampled to L samples.

Inputs
  --bam    aligned BAM with mv/ts tags (dorado --emit-moves --reference)
  --pod5   pod5 file or directory
  --sites  group=BED[,group=BED...]   BED: chrom, pos0, [end, name, score, strand]
           if a BED has no strand column the strand is inferred from the
           reference base and --mod-base for that group (needs --ref)
  --label-map group=class[,...]  what a read of THIS sample carries at sites of
           each group (e.g. wild type: dam=6mA; dam- mutant: dam=none)

Output npz: sig (N,W,L) f16, dwell (N,W) f16, base (N,W) i8 [0..3, 4=N],
  label (N,) str, group (N,) str, site (N,) i32 -> rows of site_chrom/site_pos/
  site_strand, read (N,) i32 -> read_ids, sample (str).

Never raises on a bad read: problems are counted and reported.
"""
import argparse, bisect, os, sys, time
from collections import Counter, defaultdict

import numpy as np
import pod5
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True); ap.add_argument("--pod5", required=True)
ap.add_argument("--sites", required=True); ap.add_argument("--label-map", required=True)
ap.add_argument("--sample", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--ref", default=None); ap.add_argument("--mod-base", default="",
                help="group=A|C[,...] for groups whose BED has no strand column")
ap.add_argument("--half-window", type=int, default=10); ap.add_argument("--L", type=int, default=10)
ap.add_argument("--max-sites", type=int, default=8000, help="per group")
ap.add_argument("--max-reads-per-site", type=int, default=12)
ap.add_argument("--max-bam-reads", type=int, default=0, help="stop after this many usable reads (0 = all)")
ap.add_argument("--min-mapq", type=int, default=10); ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
t0 = time.time()
HW, L, W = a.half_window, a.L, 2 * a.half_window + 1
rng = np.random.default_rng(a.seed)
kv = lambda s: dict(x.split("=", 1) for x in s.split(",") if x)
label_map, mod_base = kv(a.label_map), kv(a.mod_base)
fa = pysam.FastaFile(a.ref) if a.ref else None

# ---------------------------------------------------------------- sites
site_chrom, site_pos, site_strand, site_group = [], [], [], []
for grp, path in kv(a.sites).items():
    if grp not in label_map: sys.exit(f"group {grp} has no entry in --label-map")
    rows = []
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 2 or not c[1].lstrip("-").isdigit(): continue
            chrom, pos = c[0], int(c[1])
            strand = c[5] if len(c) > 5 and c[5] in "+-" and c[5] else None
            if strand is None:
                mb = mod_base.get(grp)
                if fa is None or mb is None: sys.exit(f"group {grp}: BED has no strand; give --ref and --mod-base {grp}=A|C")
                try: rb = fa.fetch(chrom, pos, pos + 1).upper()
                except (KeyError, ValueError): continue
                comp = {"A": "T", "C": "G"}[mb]
                strand = "+" if rb == mb else "-" if rb == comp else None
                if strand is None: continue
            rows.append((chrom, pos, strand))
    rows = sorted(set(rows))
    if len(rows) > a.max_sites:
        rows = [rows[i] for i in sorted(rng.choice(len(rows), a.max_sites, replace=False))]
    for chrom, pos, strand in rows:
        site_chrom.append(chrom); site_pos.append(pos); site_strand.append(strand); site_group.append(grp)
    print(f"[{a.sample}] sites {grp}: {len(rows):,} -> label '{label_map[grp]}'", flush=True)
if not site_pos: sys.exit("no sites loaded")

index = defaultdict(lambda: ([], []))           # (chrom, strand) -> (sorted positions, site ids)
for sid in np.argsort(site_pos, kind="stable"):
    k = (site_chrom[sid], site_strand[sid]); index[k][0].append(site_pos[sid]); index[k][1].append(int(sid))

# ---------------------------------------------------------------- pass 1: BAM
per_site = Counter(); plan = {}; skip = Counter(); n_used = 0
with pysam.AlignmentFile(a.bam, check_sq=False) as bam:
    for r in bam:
        if r.is_unmapped or r.is_secondary or r.is_supplementary: skip["not primary"] += 1; continue
        if r.mapping_quality < a.min_mapq: skip["mapq"] += 1; continue
        if not r.has_tag("mv"): skip["no mv tag"] += 1; continue
        if r.has_tag("pi"): skip["split read"] += 1; continue
        ct = r.cigartuples
        if not ct or ct[0][0] == 5 or ct[-1][0] == 5: skip["hard clip"] += 1; continue
        strand = "-" if r.is_reverse else "+"
        k = (r.reference_name, strand)
        if k not in index: skip["no sites on contig/strand"] += 1; continue
        poss, sids = index[k]
        lo, hi = bisect.bisect_left(poss, r.reference_start), bisect.bisect_left(poss, r.reference_end)
        if lo == hi: skip["covers no site"] += 1; continue
        n = r.query_length
        mv = r.get_tag("mv")
        if len(mv) < 2: skip["empty mv"] += 1; continue
        refpos = np.array([-1 if p is None else p for p in r.get_reference_positions(full_length=True)], dtype=np.int64)
        if len(refpos) != n: skip["length mismatch"] += 1; continue
        hits = []
        for j in range(lo, hi):
            sid = sids[j]
            if per_site[sid] >= a.max_reads_per_site: continue
            q = np.flatnonzero(refpos == poss[j])
            if len(q) != 1: continue                       # deleted in this read
            fi = int(n - 1 - q[0]) if r.is_reverse else int(q[0])
            if fi - HW < 0 or fi + HW >= n: continue
            hits.append((sid, fi)); per_site[sid] += 1
        if not hits: skip["no usable site"] += 1; continue
        plan[r.query_name] = (int(mv[0]), np.flatnonzero(np.asarray(mv[1:], dtype=np.int8)),
                              int(r.get_tag("ts")) if r.has_tag("ts") else 0, n, r.get_forward_sequence().upper(), hits)
        n_used += 1
        if a.max_bam_reads and n_used >= a.max_bam_reads: break
print(f"[{a.sample}] pass 1: {n_used:,} reads carry {sum(len(v[5]) for v in plan.values()):,} windows "
      f"({time.time()-t0:.0f}s); skipped: {dict(skip)}", flush=True)
if not plan: sys.exit("no read covers any site (wrong reference, strand, or BAM without mv tags?)")

# ---------------------------------------------------------------- pass 2: signal
BASE = {"A": 0, "C": 1, "G": 2, "T": 3}
N = sum(len(v[5]) for v in plan.values())
sig = np.zeros((N, W, L), np.float16); dwell = np.zeros((N, W), np.float16); base = np.full((N, W), 4, np.int8)
site = np.zeros(N, np.int32); read = np.zeros(N, np.int32); read_ids = []; k = 0; bad = Counter()
xs = np.linspace(0, 1, L)
ids = list(plan)
paths = [a.pod5] if os.path.isfile(a.pod5) else sorted(os.path.join(dp, f) for dp, _, fs in os.walk(a.pod5) for f in fs if f.endswith(".pod5"))
seen = set()
for path in paths:
    with pod5.Reader(path) as rd:
        for i in range(0, len(ids), 2000):
            for rec in rd.reads(selection=ids[i:i + 2000], missing_ok=True):
                rid = str(rec.read_id)
                if rid in seen: continue
                seen.add(rid)
                try:
                    stride, moves, ts, n, seq, hits = plan[rid]
                    if len(moves) != n: bad["moves != bases"] += 1; continue
                    x = (rec.signal.astype(np.float32) + rec.calibration.offset) * rec.calibration.scale
                    body = x[ts:]
                    med = np.median(body); mad = np.median(np.abs(body - med)) * 1.4826
                    if not np.isfinite(mad) or mad < 1e-6: bad["flat signal"] += 1; continue
                    x = np.clip((x - med) / mad, -6, 6)
                    bnd = np.append(moves * stride + ts, len(x))
                    ri = len(read_ids); read_ids.append(rid)
                    for sid, fi in hits:
                        for w, b in enumerate(range(fi - HW, fi + HW + 1)):
                            s, e = int(bnd[b]), int(bnd[b + 1])
                            if e <= s: continue
                            seg = x[s:e]
                            sig[k, w] = seg[0] if e - s == 1 else np.interp(xs, np.linspace(0, 1, e - s), seg)
                            dwell[k, w] = np.log1p(e - s)
                            base[k, w] = BASE.get(seq[b], 4)
                        site[k] = sid; read[k] = ri; k += 1
                except Exception as ex:                      # one bad read must not cost the job
                    bad[type(ex).__name__] += 1
missing = len(ids) - len(seen)
sig, dwell, base, site, read = sig[:k], dwell[:k], base[:k], site[:k], read[:k]
grp = np.array(site_group)[site]; lab = np.array([label_map[g] for g in grp])
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
tmp = a.out + ".tmp.npz"
np.savez_compressed(tmp, sig=sig, dwell=dwell, base=base, label=lab, group=grp, site=site, read=read,
                    read_ids=np.array(read_ids), site_chrom=np.array(site_chrom), site_pos=np.array(site_pos, np.int64),
                    site_strand=np.array(site_strand), sample=np.array(a.sample))
os.replace(tmp, a.out)
print(f"[{a.sample}] wrote {k:,} windows from {len(read_ids):,} reads to {a.out} ({os.path.getsize(a.out)/1e6:.0f} MB, "
      f"{time.time()-t0:.0f}s); reads missing from pod5: {missing:,}; bad: {dict(bad)}")
print(f"[{a.sample}] label counts: {dict(Counter(lab.tolist()))}")
print(f"[{a.sample}] centre base by label: " + "; ".join(
    f"{l}: {dict(Counter('ACGTN'[b] for b in base[lab == l, HW].tolist()))}" for l in sorted(set(lab.tolist()))))
if k == 0: sys.exit("no windows extracted")
