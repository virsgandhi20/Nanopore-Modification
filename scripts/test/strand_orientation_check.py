#!/usr/bin/env python3
"""Which base sequence does the raw signal of a REVERSE-mapped read follow?

rawmod/featurization.py::get_ref_info_from_bam returns, for a reverse read, the
+ strand reference bases in reversed ORDER (no complement); an earlier version
reverse-COMPLEMENTED them. The two differ in every base, so only one of them can
agree with the pore. This measures it directly, the same way the featurizer's
own auto_detect_center_idx does: Pearson r between the observed mean current of
each base and the expected k-mer level, under each convention.

  ref_order   + strand reference bases, walked in the read's 5'->3' order
              (what get_ref_info_from_bam returns today)
  pore_frame  the same, complemented for reverse reads = the strand that
              actually went through the pore (the older _revcomp behaviour)

For forward reads the two are identical (control). Segmentation comes from the
Dorado move table, so no Remora refinement is needed. Login-node safe.
"""
import argparse, statistics as st
import numpy as np, pod5, pysam

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True); ap.add_argument("--pod5", required=True)
ap.add_argument("--level-table", required=True)
ap.add_argument("--ref", default=None, help="indexed FASTA, used only if the BAM has no MD tag")
ap.add_argument("--n", type=int, default=150, help="reads per strand")
ap.add_argument("--min-mapq", type=int, default=60)
a = ap.parse_args()

levels = {}
for line in open(a.level_table):
    p = line.split()
    if len(p) >= 2 and p[0] != "kmer" and not line.startswith("#"):
        try: levels[p[0].upper()] = float(p[1])
        except ValueError: pass
K = len(next(iter(levels)))
print(f"level table: {len(levels):,} {K}-mers")
COMP = str.maketrans("ACGTN", "TGCAN")

fa = pysam.FastaFile(a.ref) if a.ref else None
want = {"+": a.n, "-": a.n}; picked = {}; n_md_missing = 0
with pysam.AlignmentFile(a.bam) as b:
    for r in b:
        s = "-" if r.is_reverse else "+"
        if want[s] <= 0 or r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        if r.mapping_quality < a.min_mapq or not r.has_tag("mv") or r.has_tag("pi"): continue
        if r.cigartuples and (r.cigartuples[0][0] == 5 or r.cigartuples[-1][0] == 5): continue
        mv = r.get_tag("mv"); stride = mv[0]; moves = np.flatnonzero(np.asarray(mv[1:], dtype=np.int8))
        n = r.query_length
        if len(moves) != n: continue
        try:
            pairs = [(q, p, c.upper()) for q, p, c in r.get_aligned_pairs(with_seq=True) if q is not None and p is not None and c]
        except ValueError:                      # no MD tag: take the + strand bases from the reference instead
            n_md_missing += 1
            if fa is None: continue
            span = fa.fetch(r.reference_name, r.reference_start, r.reference_end).upper()
            pairs = [(q, p, span[p - r.reference_start]) for q, p in r.get_aligned_pairs() if q is not None and p is not None]
        picked[r.query_name] = (s, stride, moves, r.get_tag("ts"), n, pairs)
        want[s] -= 1
        if want["+"] <= 0 and want["-"] <= 0: break

print(f"reads picked: {sum(1 for v in picked.values() if v[0]=='+')} forward, {sum(1 for v in picked.values() if v[0]=='-')} reverse"
      + (f"   (BAM has no MD tag on {n_md_missing} reads; reference bases taken from --ref)" if n_md_missing else ""))
if not picked: raise SystemExit("no usable reads (need mv tags, MAPQ >= min, and MD tags or --ref)")
res = {("+", "ref_order"): [], ("+", "pore_frame"): [], ("-", "ref_order"): [], ("-", "pore_frame"): []}
best_c = {k: [] for k in res}
with pod5.Reader(a.pod5) as rd:
    for rec in rd.reads(selection=list(picked), missing_ok=True):
        s, stride, moves, ts, n, pairs = picked[str(rec.read_id)]
        sig = rec.signal.astype(np.float64)
        bnd = np.append(moves * stride + ts, len(sig))
        obs_fwd = np.array([sig[bnd[i]:bnd[i + 1]].mean() if bnd[i + 1] > bnd[i] else np.nan for i in range(n)])
        # reference walk in the read's own 5'->3' order, exactly like get_ref_info_from_bam
        ref = {p: c for _, p, c in pairs}
        lo, hi = min(ref), max(ref)
        rpos = list(range(lo, hi + 1)); rseq = "".join(ref.get(p, "N") for p in rpos)
        obs = {p: obs_fwd[(n - 1 - q) if s == "-" else q] for q, p, _ in pairs}
        if s == "-": rpos, rseq = rpos[::-1], rseq[::-1]
        for name, seq in (("ref_order", rseq), ("pore_frame", rseq.translate(COMP) if s == "-" else rseq)):
            best = (-2, None)
            for c in range(K):
                o, e = [], []
                for i in range(c, len(seq) - K + c + 1):
                    km = seq[i - c:i - c + K]; v = obs.get(rpos[i])
                    if v is None or np.isnan(v) or km not in levels: continue
                    o.append(v); e.append(levels[km])
                if len(o) > 200:
                    rr = float(np.corrcoef(o, e)[0, 1])
                    if rr > best[0]: best = (rr, c)
            if best[1] is not None: res[(s, name)].append(best[0]); best_c[(s, name)].append(best[1])

print(f"\n{'read strand':12s} {'convention':12s} {'reads':>6s} {'median r':>9s} {'min':>7s} {'max':>7s}   best k-mer offset")
for (s, name), v in res.items():
    if v: print(f"{s:12s} {name:12s} {len(v):6d} {st.median(v):9.3f} {min(v):7.3f} {max(v):7.3f}   {st.mode(best_c[(s, name)])}")
print("\nReading: the convention whose r on '-' reads matches the '+' control is the one the signal follows.")
