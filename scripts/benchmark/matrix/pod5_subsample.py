#!/usr/bin/env python3
"""Random subset of a pod5 file: pod5_subsample.py <in.pod5> <out.pod5> <fraction> [seed]
Two-sample rows need matched depth on both samples (a 5x deeper modified sample lets rare miscalls clear the coverage
floor on one side only), so the oligo samples that were not part of a fine-tune split get the same 20% as those that were."""
import random, sys
import pod5
src, dst, frac = sys.argv[1], sys.argv[2], float(sys.argv[3]); seed = int(sys.argv[4]) if len(sys.argv) > 4 else 0
with pod5.Reader(src) as r:
    ids = [str(x) for x in r.read_ids]
random.Random(seed).shuffle(ids); keep = set(ids[: int(len(ids) * frac)])
n = 0
with pod5.Reader(src) as r, pod5.Writer(dst) as w:
    for rec in r.reads(selection=list(keep)):
        w.add_read(rec.to_read()); n += 1
print(f"{dst}: {n:,} of {len(ids):,} reads (fraction {frac}, seed {seed})")
