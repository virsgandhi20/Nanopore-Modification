#!/usr/bin/env python3
"""Positions of methylated bases for a set of motif specs, as a 0-based BED (chrom, pos, pos+1).

A spec is  motif:<IUPAC>:<offset>:<+|both>  (as in rows.tsv): the base at <offset> of every motif match on the plus strand,
and with 'both' the complementary base of every match of the reverse complement (a site on the minus strand, recorded at its
plus-strand coordinate, as the collection's gt_minus files do). --list reads the specs from motifs_hpylori_smrt.tsv.

  motif_sites.py --ref ref.fa --list motifs_hpylori_smrt.tsv --strain 26695 --type 4mC --out smrt_4mC_26695.bed
  motif_sites.py --ref ref.fa --spec motif:GAAGA:3:both --spec motif:CATG:1:both --out sites.bed
"""
import argparse, re, sys
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
         "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}
def revcomp(m): return "".join(COMP[c] for c in reversed(m))
def spec_patterns(spec):
    kind, motif, off, strands = spec.split(":"); off = int(off)
    assert kind == "motif" and strands in ("+", "both"), spec
    pats = [(re.compile("(?=" + "".join(IUPAC[c] for c in motif) + ")"), off)]
    if strands == "both": pats.append((re.compile("(?=" + "".join(IUPAC[c] for c in revcomp(motif)) + ")"), len(motif) - 1 - off))
    return pats
def specs_from_list(path, strain, mtype):
    out = []
    for line in open(path):
        if line.startswith("#") or not line.strip(): continue
        c = line.rstrip("\n").split("\t")
        if c[0] == strain and c[3] == mtype and c[6] == "yes": out.append(f"motif:{c[1]}:{c[2]}:both")
    return out
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True); ap.add_argument("--spec", action="append", default=[])
    ap.add_argument("--list"); ap.add_argument("--strain"); ap.add_argument("--type")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    specs = list(a.spec) + (specs_from_list(a.list, a.strain, a.type) if a.list else [])
    if not specs: sys.exit("no specs")
    import pysam
    fa = pysam.FastaFile(a.ref); pos = set(); per = {}
    for s in specs:
        n = 0
        for ctg in fa.references:
            seq = fa.fetch(ctg).upper()
            for pat, o in spec_patterns(s):
                for m in pat.finditer(seq): pos.add((ctg, m.start() + o)); n += 1
        per[s] = n
    with open(a.out, "w") as f:
        for ctg, p in sorted(pos): f.write(f"{ctg}\t{p}\t{p + 1}\n")
    for s, n in per.items(): print(f"{s}\t{n:,} sites", file=sys.stderr)
    print(f"{len(pos):,} positions -> {a.out}", file=sys.stderr)
if __name__ == "__main__": main()
