#!/usr/bin/env python3
"""Score fine-tuned UniMeth 5hmU calls: every T in bc07 (native SPO1) is modified,
every T in bc01 (PCR) is not. Per-read AUROC over all calls, then site-level via
the shared scorer with the bc01 sites as NEG_-prefixed negatives."""
import argparse, subprocess, sys
import numpy as np
from sklearn.metrics import roc_auc_score
ap = argparse.ArgumentParser()
for k in ("pos-calls", "neg-calls", "pos-sites", "neg-sites", "scorer", "out"): ap.add_argument("--" + k, required=True)
ap.add_argument("--min-cov", default="10")
a = ap.parse_args()
def probs(path):
    p = []
    for line in open(path):
        c = line.rstrip("\n").split("\t")
        if len(c) >= 10 and c[6] == "[5hmU]":
            try: p.append(float(c[8]))
            except ValueError: pass
    return np.array(p)
pp, pn = probs(a.pos_calls), probs(a.neg_calls)
lines = [f"per-read calls: bc07 (5hmU) {len(pp):,} mean P {pp.mean():.4f} frac>0.5 {(pp>0.5).mean():.4f} | bc01 (none) {len(pn):,} mean P {pn.mean():.4f} frac>0.5 {(pn>0.5).mean():.4f}",
         f"per-read AUROC (bc07 vs bc01): {roc_auc_score(np.r_[np.ones(len(pp)), np.zeros(len(pn))], np.r_[pp, pn]):.4f}"]
# site level: GT = every site reported in bc07; negatives = bc01 sites renamed NEG_<chrom>
comb, gt = a.out + ".sites.tsv", a.out + ".gt.bed"
with open(comb, "w") as fo, open(gt, "w") as fg:
    for line in open(a.pos_sites):
        if line.startswith("chrom"): continue
        c = line.split("\t"); fo.write(line); fg.write(f"{c[0]}\t{c[1]}\t{int(c[1])+1}\n")
    for line in open(a.neg_sites):
        if line.startswith("chrom"): continue
        fo.write("NEG_" + line)
npos = sum(1 for _ in open(gt)); lines.append(f"site-level GT: {npos:,} bc07 T sites; negatives: bc01 sites as NEG_ contigs; coverage floor {a.min_cov}")
for tag, extra in (("call_freq", []), ("mean_P", ["--num-col", "5"])):
    r = subprocess.run([sys.executable, a.scorer, "--calls", comb, "--gt", gt, "--min-cov", a.min_cov, "--chrom-col", "0", "--pos-col", "1", "--cov-col", "8", "--freq-col", "9", "--label", f"SPO1_5hmU/UniMeth_finetuned/{tag}"] + extra, capture_output=True, text=True)
    lines.append("site-level " + tag + ": " + (r.stdout.strip().splitlines()[-1] if r.stdout.strip() else "FAILED " + r.stderr.strip()[-200:]))
open(a.out, "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
