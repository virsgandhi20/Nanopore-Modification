#!/usr/bin/env python3
"""
How often does ONE site carry MORE THAN ONE modification type?

Motivation: RawMod assigns a single label per site (and chem_array types a
SPO1 site by its DOMINANT modkit code). That is fine for binary detection but
ill-defined for typing if reads at the same base carry different
modifications (e.g. 5mC on some molecules, 5hmC on others). This script puts
a number on that, from two sources:

  --bedmethyl  modkit pileup output(s): per (contig,pos,strand), the percent
               modified for each code (m=5mC, h=5hmC, a=6mA). A site is
               "mixed" when >=2 codes each exceed --mix-frac at coverage
               >= --min-cov.
  --gt         our differential HP beds (gt_4mC/gt_5mC/gt_6mA): a position
               present in >=2 type sets is mixed.

Prints counts and fractions; writes an optional TSV of mixed sites.
"""
import argparse
import collections
import sys

CODES = {"m": "5mC", "h": "5hmC", "a": "6mA"}


def read_bedmethyl(paths, min_cov):
    sites = collections.defaultdict(dict)   # (contig,pos,strand) -> {type: frac}
    for p in paths:
        with open(p) as f:
            for line in f:
                c = line.rstrip("\n").split("\t")
                if len(c) < 11:
                    continue
                code = c[3]
                if code not in CODES:
                    continue
                cov, frac = int(c[9]), float(c[10])
                if cov < min_cov:
                    continue
                key = (c[0], int(c[1]), c[5])
                t = CODES[code]
                sites[key][t] = max(frac, sites[key].get(t, 0.0))
    return sites


def read_gt(paths):
    sites = collections.defaultdict(set)
    for p in paths:
        with open(p) as f:
            for line in f:
                c = line.split("\t")
                if len(c) < 4:
                    continue
                sites[(c[0], int(c[1]), c[5].strip() if len(c) > 5 else ".")].add(c[3])
    return sites


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bedmethyl", nargs="*", default=[])
    ap.add_argument("--gt", nargs="*", default=[])
    ap.add_argument("--min-cov", type=int, default=10)
    ap.add_argument("--call-frac", type=float, default=50.0,
                    help="percent for a code to count as 'present' at a site")
    ap.add_argument("--mix-frac", type=float, default=20.0,
                    help="percent each of >=2 codes must reach to call a site mixed")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.bedmethyl:
        sites = read_bedmethyl(a.bedmethyl, a.min_cov)
        called = {k: v for k, v in sites.items() if max(v.values()) >= a.call_frac}
        mixed = {k: v for k, v in called.items()
                 if sum(1 for x in v.values() if x >= a.mix_frac) >= 2}
        pairs = collections.Counter(
            "+".join(sorted(t for t, x in v.items() if x >= a.mix_frac)) for v in mixed.values())
        print(f"[bedmethyl] covered sites (cov>={a.min_cov}): {len(sites):,}")
        print(f"[bedmethyl] called modified (any code >= {a.call_frac}%): {len(called):,}")
        print(f"[bedmethyl] MIXED (>=2 codes each >= {a.mix_frac}%): {len(mixed):,} "
              f"= {100.0 * len(mixed) / max(1, len(called)):.2f}% of called sites")
        for k, n in pairs.most_common():
            print(f"    {k}: {n:,}")
        if a.out:
            with open(a.out, "w") as fo:
                fo.write("contig\tpos\tstrand\ttypes\n")
                for (c, p, s), v in sorted(mixed.items()):
                    fo.write(f"{c}\t{p}\t{s}\t" + ";".join(f"{t}={x:.0f}" for t, x in sorted(v.items())) + "\n")
            print(f"wrote {a.out}")

    if a.gt:
        gt = read_gt(a.gt)
        multi = {k: v for k, v in gt.items() if len(v) >= 2}
        pairs = collections.Counter("+".join(sorted(v)) for v in multi.values())
        print(f"[gt] labeled positions: {len(gt):,}")
        print(f"[gt] positions in >=2 type sets: {len(multi):,} "
              f"= {100.0 * len(multi) / max(1, len(gt)):.2f}%")
        for k, n in pairs.most_common():
            print(f"    {k}: {n:,}")


if __name__ == "__main__":
    main()
