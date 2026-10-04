#!/usr/bin/env python3
"""Does UniMeth's reported reference position agree with the CIGAR projection of its own read position?

  unimeth_pos_check.py <calls.txt> <sub.bam> <ref.fa> [--type '[m6A]'] [--max-lines 3000000]

Per strand it tallies, for every per-read call, pos - proj(read_pos) under two read-position conventions
(bam: read_pos indexes the BAM SEQ; read: read_pos indexes the read as sequenced, i.e. L-1-read_pos for minus
reads), the reference base under pos for the agreeing and the disagreeing calls, and a per-read histogram of the
agreement fraction (bimodal = some reads wholly misplaced; uniform = per-position drift). Reads absent from the
BAM's primary alignments and read positions inside insertions are counted separately.
"""
import argparse, collections, sys
import numpy as np, pysam
ap = argparse.ArgumentParser()
ap.add_argument("calls"); ap.add_argument("bam"); ap.add_argument("ref")
ap.add_argument("--type", default="[m6A]"); ap.add_argument("--max-lines", type=int, default=3000000)
a = ap.parse_args()
aln = {}
with pysam.AlignmentFile(a.bam, "rb", check_sq=False) as b:
    for r in b.fetch(until_eof=True):
        if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
        m = np.full(r.query_length, -1, dtype=np.int64)
        for qq, pp in r.get_aligned_pairs(matches_only=True): m[qq] = pp
        aln[r.query_name] = (r.reference_name, "-" if r.is_reverse else "+", m)
print(f"primary alignments: {len(aln):,}")
fa = pysam.FastaFile(a.ref)
D = {"+": collections.Counter(), "-": collections.Counter()}       # delta tallies, bam convention
R = {"+": collections.Counter(), "-": collections.Counter()}       # delta tallies, read convention
base = {"+": collections.defaultdict(collections.Counter), "-": collections.defaultdict(collections.Counter)}
per_read = collections.defaultdict(lambda: [0, 0]); strand_mismatch = missing = n = neg = 0
basecache = {}
def refbase(ch, p):
    if p < 0 or ch not in fa.references: return "?"
    k = (ch, p // 100000)
    if k not in basecache:
        basecache[k] = fa.fetch(ch, k[1] * 100000, min(fa.get_reference_length(ch), k[1] * 100000 + 100000)).upper()
    s = basecache[k]; i = p - k[1] * 100000
    return s[i] if 0 <= i < len(s) else "?"
with open(a.calls) as f:
    for ln, line in enumerate(f):
        if ln >= a.max_lines: break
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[6] != a.type: continue
        try: pos = int(c[1]); rp = int(c[5])
        except ValueError: continue
        n += 1; ch, s, rid = c[0], c[2], c[4]
        if pos < 0: neg += 1
        al = aln.get(rid)
        if al is None: missing += 1; continue
        rch, rs, m = al; L = len(m)
        if rs != s: strand_mismatch += 1
        for tally, q in ((D[s], rp), (R[s], (L - 1 - rp) if rs == "-" else rp)):
            if not (0 <= q < L): tally["read_pos out of range"] += 1; continue
            pr = int(m[q])
            if pr < 0: tally["in insertion"] += 1; continue
            d = pos - pr
            tally[d if abs(d) <= 2 else ("<-2" if d < 0 else ">2")] += 1
            if tally is R[s]:   # the read-as-sequenced convention is UniMeth's (read_pos counts from the read's first base on both strands)
                base[s]["agree" if d == 0 else "disagree"][refbase(ch, pos)] += 1
                pr_ok = per_read[rid]; pr_ok[1] += 1; pr_ok[0] += (d == 0)
print(f"calls of type {a.type}: {n:,}; read not in BAM primaries: {missing:,}; strand differs from BAM: {strand_mismatch:,}; negative pos: {neg:,}")
for s in "+-":
    tot = sum(D[s].values()) or 1
    print(f"\nstrand {s}: {tot:,} calls")
    print("  delta = pos - CIGAR(read_pos), read_pos indexing the BAM SEQ:   " + ", ".join(f"{k}: {v:,} ({100*v/tot:.1f}%)" for k, v in sorted(D[s].items(), key=lambda kv: -kv[1])))
    print("  same with read_pos indexing the read as sequenced:               " + ", ".join(f"{k}: {v:,} ({100*v/tot:.1f}%)" for k, v in sorted(R[s].items(), key=lambda kv: -kv[1])))
    for cls in ("agree", "disagree"):
        bc = base[s][cls]; t = sum(bc.values()) or 1
        print(f"  reference base under pos, {cls:8} ({t:,}): " + ", ".join(f"{b}={100*v/t:.1f}%" for b, v in sorted(bc.items(), key=lambda kv: -kv[1])[:5]))
h = collections.Counter()
for ok, t in per_read.values():
    if t >= 20: h[min(int(10 * ok / t), 9)] += 1
print("\nper-read agreement fraction under the read-as-sequenced convention (reads with >= 20 calls), decile histogram 0.0-0.1 ... 0.9-1.0:")
print("  " + "  ".join(f"{k/10:.1f}-{(k+1)/10:.1f}: {h.get(k,0):,}" for k in range(10)))
