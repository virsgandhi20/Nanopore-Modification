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
# declared modification per entry, read from the preset comments (as audit_presets.py does)
import re
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "motif_gt.py")).read().splitlines()
def declared(motif, off):
    for ln in src:
        if f"'{motif}'" in ln and re.search(rf",\s*{off}\s*,", ln) and "#" in ln:
            m = re.search(r"(6mA|m6A|4mC|m4C|5mC|m5C)", ln.split("#", 1)[1])
            return m.group(1) if m else None
comp = str.maketrans("ACGTNWSMKRYBDHV", "TGCANWSKMYRVHDB")

plus, minus, slip = {}, {}, 0
self_covered = set()
for motif, off, strand in M._PRESETS[a.preset]:
    rx = "".join(iupac.get(c, c) for c in motif)
    d = declared(motif, off)
    minus_base = {"A": "T", "C": "G"}
    legit_minus = None if d is None else minus_base["A" if "A" in d else "C"]
    palin = motif.translate(comp)[::-1] == motif
    mirror_ok = palin and abs((len(motif) - 1 - off) - off) <= a.half_window
    for chrom, seq in seqs.items():
        for c, pos, st_ in M.find_motif_positions(seq, rx, off, strand, chrom, 0):
            base = seq[pos]
            if st_ == "-" or (base in "TG" and (legit_minus is None or base == legit_minus)):
                minus.setdefault(c, set()).add(pos)
                if mirror_ok or st_ == "-" and palin: self_covered.add((c, pos))
            elif base in "TG":
                slip += 1          # labelled base cannot be the declared mod on either strand
            else:
                plus.setdefault(c, set()).add(pos)
n_plus = sum(len(v) for v in plus.values()); n_minus = sum(len(v) for v in minus.values())
covered = 0
for c, ps in minus.items():
    srt = sorted(plus.get(c, ()))
    for p in ps:
        i = bisect.bisect_left(srt, p - a.half_window)
        covered += ((c, p) in self_covered) or (i < len(srt) and srt[i] <= p + a.half_window)
tot = n_plus + n_minus + slip
print(f"{a.preset}: {tot:,} GT labels = {n_plus:,} on the plus strand + {n_minus:,} on the minus strand + {slip:,} on a base that cannot carry the declared mod (offset slip)")
if n_minus:
    print(f"   minus-strand labels with a plus-strand modified base within +/-{a.half_window}: {covered:,} ({100*covered/n_minus:.1f}%)"
          f"   -> {n_minus - covered:,} labels ({100*(n_minus-covered)/tot:.1f}% of all) show forward reads NO modified base")
