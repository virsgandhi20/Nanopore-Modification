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
FAMILY = {"5mC": "C", "5hmC": "C", "6mA": "A"}   # canonical base each code modifies


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
    ap.add_argument("--control", nargs="*", default=[],
                    help="bedMethyl from an UNMODIFIED sample (e.g. PCR); sites it also calls are dropped")
    ap.add_argument("--ctrl-max", type=float, default=10.0,
                    help="max total modified percent allowed in the control")
    ap.add_argument("--min-cov", type=int, default=10)
    ap.add_argument("--call-frac", type=float, default=50.0,
                    help="total modified percent for a site to count as modified")
    ap.add_argument("--mix-frac", type=float, default=20.0,
                    help="percent each of >=2 codes must reach to call a site mixed")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.bedmethyl:
        sites = read_bedmethyl(a.bedmethyl, a.min_cov)
        ctrl = read_bedmethyl(a.control, a.min_cov) if a.control else {}
        print(f"[bedmethyl] covered sites (cov>={a.min_cov}): {len(sites):,}")
        # A true mixture can only happen WITHIN one canonical base: 5mC and 5hmC
        # are both on C; 6mA is alone on A. Codes from different families at one
        # position come from reads that disagree on the base (basecalling
        # mismatches), not from co-occurring modifications, so they are counted
        # separately as artifacts and never enter the mixture fraction.
        cross = sum(1 for v in sites.values()
                    if len({FAMILY[t] for t, x in v.items() if x >= a.mix_frac}) >= 2)
        print(f"[bedmethyl] cross-base co-occurrence (artifact, reads disagree on the base): {cross:,}")
        mixed_out = {}
        for fam in sorted(set(FAMILY.values())):
            members = [t for t in FAMILY if FAMILY[t] == fam]
            fam_sites = {}
            for k, v in sites.items():
                fv = {t: x for t, x in v.items() if t in members}
                if not fv:
                    continue
                if ctrl:   # differential: drop sites the control also calls
                    cv = sum(x for t, x in ctrl.get(k, {}).items() if t in members)
                    if cv > a.ctrl_max:
                        continue
                fam_sites[k] = fv
            called = {k: v for k, v in fam_sites.items() if sum(v.values()) >= a.call_frac}
            mixed = {k: v for k, v in called.items()
                     if sum(1 for x in v.values() if x >= a.mix_frac) >= 2}
            tag = " (after control filter)" if ctrl else ""
            print(f"[{fam}-family: {'/'.join(members)}]{tag} called modified: {len(called):,}; "
                  f"MIXED: {len(mixed):,} = {100.0 * len(mixed) / max(1, len(called)):.2f}%")
            mixed_out.update(mixed)
        if a.out:
            with open(a.out, "w") as fo:
                fo.write("contig\tpos\tstrand\ttypes\n")
                for (c, p_, s_), v in sorted(mixed_out.items()):
                    fo.write(f"{c}\t{p_}\t{s_}\t" + ";".join(f"{t}={x:.0f}" for t, x in sorted(v.items())) + "\n")
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
