#!/usr/bin/env python3
"""
Site-level scoring of an external caller's per-site output against a
ground-truth BED, for the RawMod paper's Table 1 (site-level AUROC/AUPRC).

Inputs
  --calls   TSV/bedMethyl of per-site calls. Column indices are configurable
            so the same scorer works for UniMeth's call_modification_frequency
            output (default: chrom=0, pos=1, cov=8, freq=9) and for modkit /
            Dorado bedMethyl (chrom=0, pos=1, cov=9, freq=10).
  --gt      BED of modified positions (contig, 0-based pos); every other
            called position of the candidate base is a negative.
  --candidates  optional BED restricting which positions are scored (e.g. all
            motif occurrences + background sites); if omitted, every called
            site with coverage >= --min-cov is scored.

The score for a site is its methylation frequency (fraction of reads called
modified). AUROC/AUPRC are computed over sites, positives = GT membership.
"""
import argparse
import gzip
import sys

import numpy as np
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score


def opener(p):
    return gzip.open(p, "rt") if p.endswith(".gz") else open(p)


def load_bed_positions(path):
    s = set()
    with opener(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            c = line.split("\t")
            s.add((c[0], int(c[1])))
    return s


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls", required=True)
    ap.add_argument("--gt", required=True)
    ap.add_argument("--candidates", default=None)
    ap.add_argument("--min-cov", type=int, default=10)
    ap.add_argument("--chrom-col", type=int, default=0)
    ap.add_argument("--pos-col", type=int, default=1)
    ap.add_argument("--cov-col", type=int, default=8)
    ap.add_argument("--freq-col", type=int, default=9)
    ap.add_argument("--num-col", type=int, default=None,
                    help="score = this column / coverage instead of freq-col (UniMeth: 5 = "
                         "prob_1_sum, i.e. mean P(mod); keeps the ranking that a 0.5-thresholded "
                         "call frequency throws away when the model is under-confident)")
    ap.add_argument("--freq-scale", type=float, default=1.0,
                    help="divide freq by this (100 for bedMethyl percent)")
    ap.add_argument("--code-col", type=int, default=None, help="bedMethyl: column holding the mod code (3)")
    ap.add_argument("--code", default=None, help="keep only rows with this mod code (a, m, h, 21839 ...)")
    ap.add_argument("--fill-missing", action="store_true",
                    help="with --candidates: every candidate position the tool did not score counts as called unmodified "
                         "(score 0), so a tool with no model for the row's chemistry lands at 0.5 instead of being N/A; "
                         "positions the tool scored below --min-cov are still excluded")
    ap.add_argument("--label", default="")
    ap.add_argument("--out", default=None, help="append one TSV row here")
    a = ap.parse_args()

    gt = load_bed_positions(a.gt)
    cand = load_bed_positions(a.candidates) if a.candidates else None

    # UniMeth's frequency table keeps + and - strand rows for the same position;
    # the GT is per position, so collapse strands: coverage-weighted frequency.
    acc = {}
    n_skip = 0
    with opener(a.calls) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "chrom", "track")):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 6: c = line.split()              # older modkit: space-separated tail
            if a.code is not None and (len(c) <= a.code_col or c[a.code_col] != a.code):
                continue
            try:
                key = (c[a.chrom_col], int(c[a.pos_col]))
                cov = float(c[a.cov_col])
                if a.num_col is not None:
                    freq = float(c[a.num_col]) / cov if cov > 0 else 0.0
                else:
                    freq = float(c[a.freq_col]) / a.freq_scale
            except (IndexError, ValueError):
                n_skip += 1; continue
            if cand is not None and key not in cand:
                continue
            c0, f0 = acc.get(key, (0.0, 0.0))
            acc[key] = (c0 + cov, f0 + cov * freq)

    keys = [k for k, (c, _) in acc.items() if c >= a.min_cov]
    p_list = [acc[k][1] / acc[k][0] for k in keys]
    n_filled = 0
    if a.fill_missing:
        if cand is None:
            sys.exit("--fill-missing needs --candidates (the set of positions every tool is judged on)")
        scored = set(keys); low = {k for k in acc if k not in scored}          # scored below the floor: excluded, not filled
        for k in cand:
            if k not in scored and k not in low:
                keys.append(k); p_list.append(0.0); n_filled += 1
    y = np.array([1 if k in gt else 0 for k in keys])
    p = np.array(p_list)
    if len(y) == 0 or y.sum() == 0 or y.sum() == len(y):
        print(f"{a.label}\tn={len(y)}\tpos={int(y.sum())}\tAUROC=nan (degenerate)", file=sys.stderr)
        sys.exit(1)
    auroc = roc_auc_score(y, p)
    auprc = average_precision_score(y, p)
    f1 = f1_score(y, (p >= 0.5).astype(int))
    row = (f"{a.label}\t{len(y)}\t{int(y.sum())}\t{y.mean():.4f}\t"
           f"{auroc:.4f}\t{auprc:.4f}\t{f1:.4f}")
    print("label\tn_sites\tn_pos\tpos_rate\tauroc\tauprc\tf1_at_0.5")
    print(row)
    if n_skip:
        print(f"(skipped {n_skip} unparsable lines)", file=sys.stderr)
    if a.fill_missing:
        print(f"(fill-missing: {n_filled} of {len(y)} candidate positions had no call and were scored 0)", file=sys.stderr)
    if a.out:
        new = not __import__("os").path.exists(a.out)
        with open(a.out, "a") as fo:
            if new:
                fo.write("label\tn_sites\tn_pos\tpos_rate\tauroc\tauprc\tf1_at_0.5\n")
            fo.write(row + "\n")


if __name__ == "__main__":
    main()
