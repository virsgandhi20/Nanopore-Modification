#!/usr/bin/env python3
"""Standardize a tool's per-site output to:  chrom  pos(0-based)  cov  call_freq  mean_P
One line per position (and per strand where the tool reports strands; score_sites.py collapses them).

  sites_std.py --tool unimeth  --types '[CpG],[CHG],[CHH]'  <calls.txt>  <out.tsv>
      per-read UniMeth TSV (chrom, pos, strand, pos_in_strand, read_id, read_pos, type, prob_0, prob_1, label);
      rows of other types are dropped (a 5mC model asked for CpG/CHG/CHH may still emit nothing else, but a
      patched clone can carry [5hmU] and CpG rows in one file). call_freq = fraction of reads with prob_1 > 0.5.
  sites_std.py --tool deepmod2  <deepmod2 output dir>  <out.tsv>      (via deepmod2_sites.py; mean_P = call_freq)
  sites_std.py --tool rockfish --bam <sub.bam> [--offset 0|-1]  <rockfish out>  <out.tsv>
      per-read TSV (read_id, pos, prob; no contig): the contig and strand come from the read's primary alignment
      in the BAM; --offset -1 if the positions turn out 1-based (rockfish_pos_check.py).
"""
import argparse, os, subprocess, sys
ap = argparse.ArgumentParser()
ap.add_argument("--tool", required=True, choices=["unimeth", "deepmod2", "rockfish"])
ap.add_argument("--types", default="", help="unimeth: comma-separated mod tokens to keep, e.g. '[m6A]'")
ap.add_argument("--bam", default=None, help="rockfish: BAM giving each read's contig and strand"); ap.add_argument("--offset", type=int, default=0)
ap.add_argument("src"); ap.add_argument("out")
a = ap.parse_args()
n_out = 0
if a.tool == "unimeth":
    keep = {t.strip() for t in a.types.split(",") if t.strip()}
    acc = {}; n_in = n_drop = 0; seen_types = {}
    with open(a.src) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            n_in += 1; t = c[6]; seen_types[t] = seen_types.get(t, 0) + 1
            if keep and t not in keep:
                n_drop += 1; continue
            try:
                p = float(c[8]); pos = int(c[1])
            except ValueError:
                continue
            k = (c[0], pos, c[2]); v = acc.get(k)
            if v is None:
                acc[k] = [1, p, 1 if p > 0.5 else 0]
            else:
                v[0] += 1; v[1] += p; v[2] += 1 if p > 0.5 else 0
    with open(a.out, "w") as fo:
        for (chrom, pos, strand), (n, ps, nc) in sorted(acc.items()):
            fo.write(f"{chrom}\t{pos}\t{n}\t{nc / n:.6f}\t{ps / n:.6f}\n"); n_out += 1
    print(f"sites_std unimeth: {n_in:,} per-read calls, types {seen_types}, kept {n_in - n_drop:,} -> {n_out:,} site/strand rows")
elif a.tool == "deepmod2":
    here = os.path.dirname(os.path.abspath(__file__)); tmp = a.out + ".dm2.tmp"
    r = subprocess.run([sys.executable, os.path.join(here, "..", "deepmod2_sites.py"), a.src, tmp], capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit("sites_std deepmod2: " + (r.stderr or r.stdout).strip())
    with open(tmp) as f, open(a.out, "w") as fo:
        for line in f:
            c = line.rstrip("\n").split("\t")
            fo.write(f"{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\t{c[3]}\n"); n_out += 1
    os.remove(tmp); print(f"sites_std deepmod2: {n_out:,} site rows ({r.stdout.strip()})")
else:
    # Rockfish per-read TSV (r10.4.1 branch): read_id, pos, prob, one line per read and CpG, no contig column. The
    # contig and strand are taken from the read's primary alignment in the BAM; --offset shifts the position
    # (-1 for a 1-based file). Values outside [0,1] (the -l logits) go through a sigmoid.
    import math
    if not a.bam: sys.exit("sites_std rockfish: --bam <sub.bam> is required (the output has no contig column)")
    import pysam
    aln = {}
    with pysam.AlignmentFile(a.bam, "rb", check_sq=False) as b:
        for r in b.fetch(until_eof=True):
            if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
            aln[r.query_name] = (r.reference_name, "-" if r.is_reverse else "+")
    acc = {}; n_in = n_miss = 0
    with open(a.src) as f:
        for line in f:
            c = line.split()
            if len(c) < 3 or c[0] == "read_id": continue
            al = aln.get(c[0])
            if al is None: n_miss += 1; continue
            try: pos = int(c[1]) + a.offset; p = float(c[2])
            except ValueError: continue
            if p < 0 or p > 1: p = 1 / (1 + math.exp(-p))
            k = (al[0], pos, al[1]); v = acc.get(k); n_in += 1
            if v is None: acc[k] = [1, p, 1 if p > 0.5 else 0]
            else: v[0] += 1; v[1] += p; v[2] += 1 if p > 0.5 else 0
    with open(a.out, "w") as fo:
        for (chrom, pos, strand), (n, ps, nc) in sorted(acc.items()):
            fo.write(f"{chrom}\t{pos}\t{n}\t{nc / n:.6f}\t{ps / n:.6f}\n"); n_out += 1
    print(f"sites_std rockfish: {n_in:,} per-read calls joined to the BAM ({n_miss:,} reads not in it), {len(aln):,} primary alignments -> {n_out:,} site/strand rows")
