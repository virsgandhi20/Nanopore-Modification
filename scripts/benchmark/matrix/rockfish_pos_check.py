#!/usr/bin/env python3
"""Which coordinate does Rockfish report? Its per-read output (read_id, pos, prob) carries no contig and, as the
M.SssI smoke test showed, the position is an index into the READ, not the reference. This projects each call onto
the reference through the read's primary alignment under two hypotheses and prints the reference base at the
projected position and its neighbours, per strand:
  bam   pos indexes the BAM's SEQ (reverse-complemented for minus reads)
  read  pos indexes the read as sequenced (minus reads: SEQ index = L-1-pos)
The right hypothesis puts ~100% of plus-read calls on a C followed by G; minus-read calls then sit on the C (bam)
or on the G with C before it (read).   rockfish_pos_check.py <rockfish out> <bam> <ref.fa>"""
import sys
import numpy as np
import pysam
out, bam, ref = sys.argv[1:4]
aln = {}
with pysam.AlignmentFile(bam, "rb", check_sq=False) as b:
    for r in b.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        aln[r.query_name] = r
fa = pysam.FastaFile(ref); seqs = {}; q2r = {}
def proj(r):
    m = q2r.get(r.query_name)
    if m is None:
        m = np.full(r.query_length, -1, dtype=np.int64)
        for q, p in r.get_aligned_pairs(matches_only=True): m[q] = p
        q2r[r.query_name] = m
    return m
counts = {h: {s: {o: {} for o in (-1, 0, 1)} for s in "+-"} for h in ("bam", "read")}; n = miss = unal = 0
for line in open(out):
    c = line.split()
    if len(c) < 3 or c[0] == "read_id": continue
    r = aln.get(c[0])
    if r is None: miss += 1; continue
    p = int(c[1]); n += 1; m = proj(r); ctg = r.reference_name; s = seqs.get(ctg)
    if s is None: s = seqs[ctg] = fa.fetch(ctg).upper()
    strand = "-" if r.is_reverse else "+"
    for h, q in (("bam", p), ("read", r.query_length - 1 - p if r.is_reverse else p)):
        if not (0 <= q < len(m)) or m[q] < 0: unal += 1; continue
        for o in (-1, 0, 1):
            rp = int(m[q]) + o; b = s[rp] if 0 <= rp < len(s) else "?"; d = counts[h][strand][o]; d[b] = d.get(b, 0) + 1
print(f"{n:,} calls, {len(aln):,} primary alignments, {miss:,} calls from reads not in the BAM, {unal:,} projections onto unaligned bases")
for h in ("bam", "read"):
    print(f"hypothesis {h}:")
    for strand in "+-":
        for o in (-1, 0, 1):
            d = counts[h][strand][o]; tot = sum(d.values()) or 1
            print(f"  {strand} reads, ref base at proj{o:+d}: " + ", ".join(f"{k} {v / tot:.0%}" for k, v in sorted(d.items(), key=lambda kv: -kv[1])))
