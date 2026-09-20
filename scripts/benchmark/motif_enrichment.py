#!/usr/bin/env python3
"""Which sequence motifs are actually methylated? (de novo, from site calls)

Takes per-site modification frequencies from any caller, splits covered sites
into methylated (freq >= hi) and unmethylated (freq <= lo), and ranks k-mers
(k = 4..6, modified base at every possible offset) by how many methylated
sites they explain and how specific they are. Also reports how a motif-derived
ground truth agrees with the calls. Used to check a GT preset against the
data with two independent callers before trusting a benchmark row.
"""
import argparse, gzip, sys
from collections import Counter

ap = argparse.ArgumentParser()
ap.add_argument("--sites", required=True); ap.add_argument("--ref", required=True)
ap.add_argument("--format", choices=["unimeth", "bedmethyl"], default="unimeth")
ap.add_argument("--code", default="a", help="bedmethyl mod code to keep")
ap.add_argument("--base", default="A"); ap.add_argument("--gt", default=None)
ap.add_argument("--min-cov", type=int, default=10)
ap.add_argument("--hi", type=float, default=0.7); ap.add_argument("--lo", type=float, default=0.1)
ap.add_argument("--label", default=""); ap.add_argument("--top", type=int, default=12)
a = ap.parse_args()

ref, name = {}, None
with (gzip.open(a.ref, "rt") if a.ref.endswith(("gz", "bgz")) else open(a.ref)) as f:
    for line in f:
        if line.startswith(">"): name = line[1:].split()[0]; ref[name] = []
        else: ref[name].append(line.strip().upper())
ref = {k: "".join(v) for k, v in ref.items()}
comp = str.maketrans("ACGTN", "TGCAN")

gt = set()
if a.gt:
    for line in open(a.gt):
        c = line.split("\t")
        if len(c) >= 2 and c[1].strip().isdigit(): gt.add((c[0], int(c[1])))

F = 6  # flank
meth, unmeth = Counter(), Counter()
n = n_m = n_u = n_wrongbase = 0
gt_cov = gt_m = m_in_gt = 0
with open(a.sites) as f:
    for line in f:
        c = line.split()
        try:
            if a.format == "unimeth":
                chrom, pos, strand, cov, freq = c[0], int(c[1]), c[2], float(c[8]), float(c[9])
            else:
                if c[3] != a.code: continue
                chrom, pos, strand, cov, freq = c[0], int(c[1]), c[5], float(c[9]), float(c[10]) / 100
        except (IndexError, ValueError):
            continue
        if cov < a.min_cov or pos < F or chrom not in ref or pos + F >= len(ref[chrom]): continue
        ctx = ref[chrom][pos - F:pos + F + 1]
        if strand == "-": ctx = ctx.translate(comp)[::-1]
        if ctx[F] != a.base: n_wrongbase += 1; continue
        n += 1
        is_m, is_u = freq >= a.hi, freq <= a.lo
        if (chrom, pos) in gt:
            gt_cov += 1; gt_m += is_m
        if is_m:
            n_m += 1; m_in_gt += (chrom, pos) in gt
        elif is_u: n_u += 1
        else: continue
        tgt = meth if is_m else unmeth
        for k in (4, 5, 6):
            for off in range(k):
                s = F - off
                tgt[(ctx[s:s + k], off)] += 1

print(f"==== {a.label}  ({a.format})")
print(f"  covered {a.base} sites: {n:,}   methylated (freq>={a.hi}): {n_m:,} ({100*n_m/max(n,1):.2f}%)   unmethylated (freq<={a.lo}): {n_u:,}   wrong-base rows skipped: {n_wrongbase:,}")
if a.gt:
    print(f"  GT preset vs calls: {gt_cov:,} GT sites covered, {gt_m:,} of them methylated ({100*gt_m/max(gt_cov,1):.1f}%);"
          f"  of all methylated sites, {m_in_gt:,} ({100*m_in_gt/max(n_m,1):.1f}%) are in the GT")
rows = []
for key, m in meth.items():
    if m < max(30, 0.01 * n_m): continue
    u = unmeth.get(key, 0)
    rows.append((m, m / (m + u), key))
rows = [r for r in rows if r[1] >= 0.8]
rows.sort(key=lambda r: (-r[0], -r[1]))
print(f"  motifs with >=80% of their occurrences methylated (modified base shown in [ ]):")
if not rows: print("    none")
kept = []
for r in rows:   # drop longer k-mers that merely extend a motif already listed
    kmer, off = r[2]
    if any(off - o2 >= 0 and kmer[off - o2:off - o2 + len(k2)] == k2 for _, _, (k2, o2) in kept): continue
    kept.append(r)
for m, prec, (kmer, off) in kept[:a.top]:
    shown = kmer[:off] + "[" + kmer[off] + "]" + kmer[off + 1:]
    print(f"    {shown:10s} methylated occurrences {m:7,}  = {100*m/max(n_m,1):5.1f}% of methylated sites   specificity {100*prec:5.1f}%")
