#!/usr/bin/env python3
"""
Base-split analysis of the 5hmU holdout from the pooled LOMO.

The base-specificity hypothesis makes a directional prediction here. SPO1's
hmU sites sit at both A and T reference positions. When 5hmU is held out,
training still contains 6mA (adenine) but no thymine modification, so the
model should transfer partially to SPO1's A positions and fail on its T
positions. This script reads the saved per-fold probabilities and scores the
two subsets separately: positives at A centers vs the shared C/G negatives,
and positives at T centers vs the same negatives.

Alignment: probs_5hmU_seed*.npz rows are barcode06 then barcode07 in CSV
order (see train_lomo_pooled.train_one). We reload the kmers in that order
and hard-verify by recomputing the labels; any mismatch aborts.
"""
import os
import glob
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score

SPO1_WGS = ["barcode06", "barcode07"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lomo-dir", default="/fs/nexus-scratch/vgandhi/lomo_pooled")
    ap.add_argument("--spo1-dir",
                    default="/fs/nexus-scratch/vgandhi/orca_feat_spo1/run1_jan31/single_end")
    ap.add_argument("--mod-bases", default="AT")
    ap.add_argument("--out", default=None,
                    help="optional TSV path for the summary")
    args = ap.parse_args()

    centers = pd.concat([
        pd.read_csv(os.path.join(args.spo1_dir, f"{bc}_sub0.05",
                                 f"{bc}.merged.feature.per.site"),
                    usecols=["kmer"])["kmer"]
        for bc in SPO1_WGS], ignore_index=True).str[4].str.upper()
    y_expect = centers.isin(list(args.mod_bases)).to_numpy().astype(np.int64)

    rows = []
    for f in sorted(glob.glob(os.path.join(args.lomo_dir, "probs_5hmU_seed*.npz"))):
        seed = f.rsplit("seed", 1)[1].split(".")[0]
        d = np.load(f)
        p, y = d["probs"], d["labels"]
        assert len(p) == len(y_expect), \
            f"{f}: {len(p)} probs vs {len(y_expect)} CSV rows -- CSVs changed?"
        assert (y == y_expect).all(), f"{f}: label mismatch, alignment broken"

        neg = ~centers.isin(list(args.mod_bases)).to_numpy()
        for base in args.mod_bases:
            pos = (centers == base).to_numpy()
            keep = pos | neg
            yb, pb = y[keep], p[keep]
            auprc = average_precision_score(yb, pb)
            rate = yb.mean()
            rows.append({"seed": seed, "base": base, "n_pos": int(yb.sum()),
                         "n_neg": int((yb == 0).sum()),
                         "auprc": auprc, "pos_rate": rate,
                         "lift": auprc / rate})
        # overall for reference
        auprc = average_precision_score(y, p)
        rows.append({"seed": seed, "base": "A+T", "n_pos": int(y.sum()),
                     "n_neg": int((y == 0).sum()), "auprc": auprc,
                     "pos_rate": y.mean(), "lift": auprc / y.mean()})

    df = pd.DataFrame(rows)
    summary = (df.groupby("base")
                 .agg(auprc_mean=("auprc", "mean"), auprc_std=("auprc", "std"),
                      lift_mean=("lift", "mean"), n_pos=("n_pos", "first"),
                      n_neg=("n_neg", "first"))
                 .round(4))
    print(summary.to_string())
    print("\nPrediction if base-specificity holds: A (6mA in training) >> T "
          "(no thymine mod in training).")
    if args.out:
        summary.to_csv(args.out, sep="\t")
        print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
