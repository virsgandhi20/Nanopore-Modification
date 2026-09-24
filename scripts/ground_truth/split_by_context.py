#!/usr/bin/env python3
"""Split a per-position ground-truth BED (contig, 0-based pos [, ...]) into CpG and non-CpG cytosines.

A position is a plus-strand cytosine when the reference base is C (CpG if the next base is G) and a
minus-strand cytosine when the reference base is G (CpG if the previous base is C). Positions on A/T
are reported and dropped. Needs an indexed FASTA (a .fai next to it; copy the reference somewhere
writable and run `samtools faidx` if the original folder is read-only).
"""
import argparse, sys
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--ref", required=True); ap.add_argument("--bed", required=True, nargs="+", help="one or more BEDs (contig, pos, ...)")
ap.add_argument("--out-cpg", required=True); ap.add_argument("--out-noncpg", required=True)
ap.add_argument("--one-based", action="store_true", help="input positions are 1-based (default 0-based)")
a = ap.parse_args()
fa = pysam.FastaFile(a.ref)
seqs = {}
def base(ctg, p):
    s = seqs.get(ctg)
    if s is None:
        s = seqs[ctg] = fa.fetch(ctg).upper()
    return s[p] if 0 <= p < len(s) else "N"
n = {"cpg": 0, "noncpg": 0, "notC": 0, "seen": set()}
with open(a.out_cpg, "w") as fc, open(a.out_noncpg, "w") as fn:
    for path in a.bed:
        for line in open(path):
            f = line.split()
            if len(f) < 2 or f[0].startswith("#"):
                continue
            ctg, p = f[0], int(f[1]) - (1 if a.one_based else 0)
            if (ctg, p) in n["seen"]:
                continue
            n["seen"].add((ctg, p))
            b = base(ctg, p)
            if b == "C":
                ctx = "cpg" if base(ctg, p + 1) == "G" else "noncpg"
            elif b == "G":
                ctx = "cpg" if base(ctg, p - 1) == "C" else "noncpg"
            else:
                n["notC"] += 1; continue
            (fc if ctx == "cpg" else fn).write(f"{ctg}\t{p}\n"); n[ctx] += 1
print(f"split_by_context: {len(n['seen']):,} positions -> CpG {n['cpg']:,}, non-CpG {n['noncpg']:,}, not a C/G {n['notC']:,}", file=sys.stderr)
if n["notC"] > 0.05 * max(1, len(n["seen"])):
    print("WARNING: more than 5 percent of positions are not on a C or G; check the coordinate convention (--one-based?)", file=sys.stderr)
