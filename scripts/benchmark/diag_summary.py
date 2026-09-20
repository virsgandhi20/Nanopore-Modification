#!/usr/bin/env python3
"""
Per-read separation of UniMeth calls at ground-truth sites vs background.

For a diagnostic we do not need per-site coverage: if the model sees the
modification at all, per-read P(mod) at known-modified positions must be
higher than at other positions of the same base. Reports mean P(mod) at GT
vs background, fraction called, and the per-read AUROC (0.5 = blind).
"""
import argparse
import sys

import numpy as np
from sklearn.metrics import roc_auc_score


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    gt = set()
    with open(a.gt) as f:
        for line in f:
            c = line.split("\t")
            if len(c) >= 2:
                gt.add((c[0], int(c[1])))

    y, p = [], []
    with open(a.calls) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 10:
                continue
            try:
                pos = int(c[1]); prob = float(c[8])
            except ValueError:
                continue
            if pos < 0:                      # soft-clipped / inserted base
                continue
            y.append(1 if (c[0], pos) in gt else 0); p.append(prob)
    y, p = np.array(y), np.array(p)
    if len(y) == 0 or y.sum() in (0, len(y)):
        row = f"{a.label}\t{len(y)}\t{int(y.sum()) if len(y) else 0}\tnan\tnan\tnan\tnan\tnan"
    else:
        row = (f"{a.label}\t{len(y)}\t{int(y.sum())}\t{p[y == 1].mean():.4f}\t{p[y == 0].mean():.4f}\t"
               f"{(p[y == 1] > 0.5).mean():.4f}\t{(p[y == 0] > 0.5).mean():.4f}\t{roc_auc_score(y, p):.4f}")
    hdr = "config\tn_calls\tn_at_gt\tmeanP_gt\tmeanP_bg\tfrac_called_gt\tfrac_called_bg\tper_read_auroc"
    import os
    new = not os.path.exists(a.out)
    with open(a.out, "a") as fo:
        if new:
            fo.write(hdr + "\n")
        fo.write(row + "\n")
    print(hdr); print(row)


if __name__ == "__main__":
    main()
