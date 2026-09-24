#!/usr/bin/env python3
"""Standardize DeepMod2's per_site output to: chrom, 0-based pos, coverage, modified fraction (0-1).

Usage: deepmod2_sites.py <deepmod2 output dir> <out.tsv> [--prefix NEG_]
DeepMod2 names its columns by version; the header is matched by name and the script stops with the
header printed if a needed column cannot be found, so the mapping can be set by hand.
"""
import glob, sys
args = [x for x in sys.argv[1:] if not x.startswith("--")]
prefix = sys.argv[sys.argv.index("--prefix") + 1] if "--prefix" in sys.argv else ""
calls_dir, out = args
files = sorted(glob.glob(f"{calls_dir}/*per_site*"))
if not files:
    sys.exit(f"deepmod2_sites: no per_site file in {calls_dir}")
need = {"chrom": ["chromosome", "chrom", "contig"], "pos": ["position_before", "start", "pos0", "position"],
        "cov": ["coverage", "cov", "total_reads"], "frac": ["mod_fraction", "mod_percentage", "fraction", "modified_fraction", "mod_frac"]}
rows = 0
with open(files[0]) as f, open(out, "w") as fo:
    hdr = f.readline().lstrip("#").rstrip("\n").split("\t"); col = {h.strip().lower(): i for i, h in enumerate(hdr)}
    idx = {k: next((col[c] for c in v if c in col), None) for k, v in need.items()}
    if None in idx.values():
        sys.exit(f"deepmod2_sites: unexpected header {hdr} -> {idx}; set the columns by hand")
    scale = 100.0 if "percentage" in hdr[idx["frac"]].lower() else 1.0
    for line in f:
        c = line.rstrip("\n").split("\t")
        try:
            fo.write(f"{prefix}{c[idx['chrom']]}\t{int(c[idx['pos']])}\t{c[idx['cov']]}\t{float(c[idx['frac']]) / scale}\n"); rows += 1
        except (ValueError, IndexError):
            pass
print(f"deepmod2_sites: {rows:,} site rows from {files[0]} (columns {idx}, scale {scale})")
