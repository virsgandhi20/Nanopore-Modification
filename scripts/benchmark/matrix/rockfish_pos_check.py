#!/usr/bin/env python3
"""Which coordinate does Rockfish report? Join its per-read output (read_id, pos, prob) with the BAM's alignments and
count the reference base at pos-1, pos, pos+1 for plus- and minus-strand reads. A 0-based C-of-CpG convention shows C
at pos (plus reads) with G at pos+1; a 1-based one shows G at pos with C at pos-1.
  rockfish_pos_check.py <rockfish out> <bam> <ref.fa>"""
import sys
import pysam
out, bam, ref = sys.argv[1:4]
aln = {}
with pysam.AlignmentFile(bam, "rb", check_sq=False) as b:
    for r in b.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        aln[r.query_name] = (r.reference_name, "-" if r.is_reverse else "+")
fa = pysam.FastaFile(ref); seqs = {}; counts = {"+": {-1: {}, 0: {}, 1: {}}, "-": {-1: {}, 0: {}, 1: {}}}; n = miss = 0; pmin = 10**12; pmax = -1
for line in open(out):
    c = line.split()
    if len(c) < 3 or c[0] == "read_id": continue
    a = aln.get(c[0])
    if a is None: miss += 1; continue
    ctg, strand = a; s = seqs.get(ctg)
    if s is None: s = seqs[ctg] = fa.fetch(ctg).upper()
    p = int(c[1]); n += 1; pmin = min(pmin, p); pmax = max(pmax, p)
    for off in (-1, 0, 1):
        q = p + off; b = s[q] if 0 <= q < len(s) else "?"; d = counts[strand][off]; d[b] = d.get(b, 0) + 1
print(f"{n:,} calls joined to {len(aln):,} primary alignments, {miss:,} calls with no alignment in the BAM; pos range {pmin}-{pmax}")
for strand in ("+", "-"):
    for off in (-1, 0, 1):
        d = counts[strand][off]; tot = sum(d.values()) or 1
        print(f"  {strand} reads, base at pos{off:+d}: " + ", ".join(f"{k} {v / tot:.0%}" for k, v in sorted(d.items(), key=lambda kv: -kv[1])))
