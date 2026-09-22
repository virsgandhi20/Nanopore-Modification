#!/usr/bin/env python3
"""Split one pod5 into train / val pod5 files by read id (UniMeth's finetuner pairs
each pod5 with the BAM that holds the same reads, so validation reads must come
from a separate pod5 of the same sample)."""
import argparse, sys
import numpy as np, pod5

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True); ap.add_argument("--train", required=True); ap.add_argument("--val", required=True)
ap.add_argument("--val-frac", type=float, default=0.1); ap.add_argument("--max-reads", type=int, default=0, help="cap total reads (0 = all)")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()
rng = np.random.default_rng(a.seed)
with pod5.Reader(a.inp) as rd:
    ids = [str(x) for x in rd.read_ids]
rng.shuffle(ids)
if a.max_reads: ids = ids[:a.max_reads]
n_val = max(1, int(len(ids) * a.val_frac)); val = set(ids[:n_val]); train = set(ids[n_val:])
counts = {"train": 0, "val": 0}
with pod5.Reader(a.inp) as rd, pod5.Writer(a.train) as wt, pod5.Writer(a.val) as wv:
    for rec in rd.reads(selection=ids):
        rid = str(rec.read_id)
        if rid in val: wv.add_read(rec.to_read()); counts["val"] += 1
        elif rid in train: wt.add_read(rec.to_read()); counts["train"] += 1
print(f"{a.inp}: {counts['train']:,} train reads -> {a.train}, {counts['val']:,} val reads -> {a.val}")
if min(counts.values()) == 0: sys.exit("empty split")
