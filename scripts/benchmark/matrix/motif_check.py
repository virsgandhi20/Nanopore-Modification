#!/usr/bin/env python3
"""Is each SMRT motif actually methylated in our sample? For every spec, the fraction of its sites that a caller's
per-site output calls modified (>= 70% of reads) among sites covered >= mincov, and the fraction called unmodified
(<= 10%). A motif whose covered sites are mostly called unmethylated is 'exclude (absent in sample)', the T. denticola / J99 preset case;
the methylated fraction a caller reaches is the caller's recall, not the truth (SMRT: 96-99.7%), so it does not exclude a motif.

  --calls  a modkit bedMethyl (code in column 4, valid coverage column 10, percent modified column 11) with --format bedmethyl,
           or a matrix sites.std.tsv (chrom, pos, cov, call_freq, mean_P) with --format std
  motif_check.py --ref ref.fa --list motifs_hpylori_smrt.tsv --strain 26695 --type 6mA --calls dorado_6mA.bed --format bedmethyl --code a
"""
import argparse, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from motif_sites import spec_patterns, specs_from_list
ap = argparse.ArgumentParser()
ap.add_argument("--ref", required=True); ap.add_argument("--list", required=True); ap.add_argument("--strain", required=True); ap.add_argument("--type", required=True)
ap.add_argument("--calls", required=True); ap.add_argument("--format", choices=["bedmethyl", "std"], default="bedmethyl"); ap.add_argument("--code", default=None)
ap.add_argument("--mincov", type=int, default=10); ap.add_argument("--hi", type=float, default=70); ap.add_argument("--lo", type=float, default=10)
a = ap.parse_args()
import pysam
freq = {}
for line in open(a.calls):
    if line.startswith(("#", "chrom", "track")): continue
    c = line.rstrip("\n").split("\t")
    try:
        if a.format == "bedmethyl":
            if a.code and c[3] != a.code: continue
            cov, pct = float(c[9]), float(c[10]); key = (c[0], int(c[1]))
        else:
            cov, pct = float(c[2]), 100 * float(c[3]); key = (c[0], int(c[1]))
    except (IndexError, ValueError): continue
    if cov >= a.mincov: freq[key] = max(pct, freq.get(key, -1))      # both strands of a position: keep the higher call
fa = pysam.FastaFile(a.ref); seqs = {ctg: fa.fetch(ctg).upper() for ctg in fa.references}
print("spec\tsites\tcovered\tpct_methylated\tpct_unmethylated\tverdict")
for s in specs_from_list(a.list, a.strain, a.type):
    sites = set()
    for ctg, seq in seqs.items():
        for pat, o in spec_patterns(s):
            for m in pat.finditer(seq): sites.add((ctg, m.start() + o))
    cov = [freq[k] for k in sites if k in freq]
    if not cov: print(f"{s}\t{len(sites):,}\t0\t-\t-\tno covered sites"); continue
    hi = sum(1 for p in cov if p >= a.hi) / len(cov); lo = sum(1 for p in cov if p <= a.lo) / len(cov)
    # absent = most covered sites called unmethylated (the T. denticola / J99 preset case); partial = a real motif the
    # caller only half detects (Dorado 4mC on GAAGA/TCTTC: 49% at >= 70%, 3% at <= 10%); keep = methylated and not absent
    verdict = "exclude (absent in sample)" if lo >= 0.5 or hi < 0.1 else ("keep" if hi >= 0.5 else "keep (caller detects it only partially)")
    print(f"{s}\t{len(sites):,}\t{len(cov):,}\t{100 * hi:.1f}\t{100 * lo:.1f}\t{verdict}")
