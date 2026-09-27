#!/usr/bin/env python3
"""Standardize a tool's per-site output to:  chrom  pos(0-based)  cov  call_freq  mean_P
One line per position (and per strand where the tool reports strands; score_sites.py collapses them).

  sites_std.py --tool unimeth  --types '[CpG],[CHG],[CHH]'  <calls.txt>  <out.tsv>
      per-read UniMeth TSV (chrom, pos, strand, pos_in_strand, read_id, read_pos, type, prob_0, prob_1, label);
      rows of other types are dropped (a 5mC model asked for CpG/CHG/CHH may still emit nothing else, but a
      patched clone can carry [5hmU] and CpG rows in one file). call_freq = fraction of reads with prob_1 > 0.5.
  sites_std.py --tool deepmod2  <deepmod2 output dir>  <out.tsv>      (via deepmod2_sites.py; mean_P = call_freq)
  sites_std.py --tool rockfish  <rockfish output>      <out.tsv>      (wired once the smoke test shows the format)
"""
import argparse, os, subprocess, sys
ap = argparse.ArgumentParser()
ap.add_argument("--tool", required=True, choices=["unimeth", "deepmod2", "rockfish"])
ap.add_argument("--types", default="", help="unimeth: comma-separated mod tokens to keep, e.g. '[m6A]'")
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
    sys.exit("sites_std rockfish: output format not wired yet (run the smoke test in rockfish_setup.sh and paste its output)")
