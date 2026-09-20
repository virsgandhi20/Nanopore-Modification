#!/usr/bin/env python3
"""What does forward-reads-only featurization (--strand +) see at each GT site?

The GT BEDs are (chrom, pos) with no strand. A site produced by a minus-strand
motif hit is modified on the MINUS strand, so forward reads carry an unmodified
base there. They can still see a modification if a plus-strand site lies inside
the +/- half-window image; otherwise the positive label has no modified signal
behind it. Counts both cases per preset. Needs only the reference.
"""
import argparse, bisect, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import motif_gt as M

ap = argparse.ArgumentParser()
ap.add_argument("--ref", required=True); ap.add_argument("--preset", required=True)
ap.add_argument("--half-window", type=int, default=10)
a = ap.parse_args()

name, chunks, seqs = None, [], {}
with M.open_fasta(a.ref) as f:
    for line in f:
        if line.startswith(">"):
            if name: seqs[name] = "".join(chunks).upper()
            name, chunks = line[1:].split()[0], []
        elif name: chunks.append(line.strip())
if name: seqs[name] = "".join(chunks).upper()

iupac = {"N": "[ACGT]", "W": "[AT]", "S": "[CG]", "R": "[AG]", "Y": "[CT]", "M": "[AC]", "K": "[GT]",
         "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]"}
plus, minus = {}, {}
for motif, off, strand in M._PRESETS[a.preset]:
    rx = "".join(iupac.get(c, c) for c in motif)
    for chrom, seq in seqs.items():
        for c, pos, s in M.find_motif_positions(seq, rx, off, strand, chrom, 0):
            base = seq[pos]
            # a '+' entry that labels T/G is the plus coordinate of a minus-strand A/C
            on_minus = (s == "-") or base in "TG"
            (minus if on_minus else plus).setdefault(c, set()).add(pos)
n_plus = sum(len(v) for v in plus.values()); n_minus = sum(len(v) for v in minus.values())
covered = 0
for c, ps in minus.items():
    srt = sorted(plus.get(c, ()))
    for p in ps:
        i = bisect.bisect_left(srt, p - a.half_window)
        covered += i < len(srt) and srt[i] <= p + a.half_window
print(f"{a.preset}: {n_plus + n_minus:,} GT sites = {n_plus:,} modified on the plus strand + {n_minus:,} on the minus strand")
if n_minus:
    print(f"   minus-strand sites with a plus-strand modified base within +/-{a.half_window}: {covered:,} ({100*covered/n_minus:.1f}%)"
          f"   -> {n_minus - covered:,} positive labels ({100*(n_minus-covered)/(n_plus+n_minus):.1f}% of all) show forward reads NO modified base")
