#!/usr/bin/env python3
"""Stranded site lists (BED6) for per-read typing on E. coli-like genomes.

  dam   GATC, the A, both strands            (6mA in dam+ strains)
  dcm   CCWGG, the inner C, both strands     (5mC in dcm+ strains)
  cpg   CG, the C, both strands              (5mC after M.SssI)
  bgA / bgC   A or C positions with NO dam/dcm site inside the +/-clear window,
              so a read window around them holds no modified base in the wild
              type or the dam-/dcm- strain (never used from the M.SssI sample)

dam/dcm/cpg sites are also dropped when another motif's site falls inside
their window, so every window carries exactly one candidate modification.
"""
import argparse, gzip, os, re
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--ref", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--clear", type=int, default=12); ap.add_argument("--n-bg", type=int, default=30000)
ap.add_argument("--motifs", default="dam,dcm,cpg"); ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
seqs, name = {}, None
with (gzip.open(a.ref, "rt") if a.ref.endswith("gz") else open(a.ref)) as f:
    for line in f:
        if line.startswith(">"): name = line[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(line.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items()}
# motif -> (regex, plus-strand offset of the modified base, minus-strand offset in + coordinates)
M = {"dam": ("(?=GATC)", 1, 2), "dcm": ("(?=CC[AT]GG)", 1, 3), "cpg": ("(?=CG)", 0, 1)}
use = a.motifs.split(",")
out = {g: [] for g in use + ["bgA", "bgC"]}
for chrom, s in seqs.items():
    n = len(s); occupied = {g: np.zeros(n, bool) for g in use}; sites = {g: [] for g in use}
    for g in use:
        rx, op, om = M[g]
        for m in re.finditer(rx, s):
            sites[g] += [(m.start() + op, "+"), (m.start() + om, "-")]
            occupied[g][m.start() + op] = occupied[g][m.start() + om] = True
    def near(mask):                      # any True within +/-clear
        c = np.concatenate([[0], np.cumsum(mask)]); i = np.arange(n)
        return (c[np.minimum(n, i + a.clear + 1)] - c[np.maximum(0, i - a.clear)]) > 0
    nearg = {g: near(occupied[g]) for g in use}
    # CpG is unmodified in every sample except the M.SssI one, and only cpg sites
    # are taken from that sample, so CpG proximity matters to nobody else.
    native = [h for h in use if h != "cpg"]
    for g in use:
        others = np.zeros(n, bool)
        for h in native:
            if h != g: others |= nearg[h]
        out[g] += [(chrom, p, st) for p, st in sites[g] if not others[p]]
    anymod = np.zeros(n, bool)
    for g in native: anymod |= nearg[g]
    arr = np.frombuffer(s.encode(), dtype="S1")
    for b, comp, grp in (("A", "T", "bgA"), ("C", "G", "bgC")):
        ok = ~anymod; ok[:a.clear] = ok[n - a.clear:] = False
        out[grp] += [(chrom, int(p), "+") for p in np.flatnonzero(ok & (arr == b.encode()))]
        out[grp] += [(chrom, int(p), "-") for p in np.flatnonzero(ok & (arr == comp.encode()))]
for g, rows in out.items():
    if g.startswith("bg") and len(rows) > a.n_bg: rows = [rows[i] for i in sorted(rng.choice(len(rows), a.n_bg, replace=False))]
    with open(os.path.join(a.out, g + ".bed"), "w") as f:
        for chrom, p, st in rows: f.write(f"{chrom}\t{p}\t{p+1}\t{g}\t0\t{st}\n")
    print(f"  {g}: {len(rows):,} sites")
