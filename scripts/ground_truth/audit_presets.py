#!/usr/bin/env python3
"""Consistency audit of the motif ground-truth presets in motif_gt.py.

Needs no data. For each (motif, offset, strand) entry it reports the base that
gets labelled and checks it against the modification the preset comment names:
a 6mA site must be labelled at an A (plus strand) or at a T (the plus-strand
coordinate of a minus-strand A); a 4mC/5mC site at a C or a G. Anything else
cannot be the modified base under either convention. Also flags palindromic
motifs that are labelled on one strand only.
"""
import os, re

HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, "motif_gt.py")).read()
i = src.index("_PRESETS = {"); j = src.index("\n}\n", i) + 3
ns = {}; exec(src[i:j], ns)
lines = src[i:j].splitlines()
comp = str.maketrans("ACGTNWSMKRYBDHV", "TGCANWSKMYRVHDB")
rc = lambda s: s.translate(comp)[::-1]

def declared_mod(motif, off):
    for ln in lines:
        if f"'{motif}'" in ln and re.search(rf",\s*{off}\s*,", ln):
            m = re.search(r"(6mA|m6A|4mC|m4C|5mC|m5C)", ln.split("#", 1)[-1] if "#" in ln else "")
            return m.group(1) if m else None

bad = 0
for name, entries in ns["_PRESETS"].items():
    print(name)
    for motif, off, strand in entries:
        base, mod = motif[off], declared_mod(motif, off)
        ok = {"A", "T"} if mod and "A" in mod else {"C", "G"} if mod else {"A", "C", "G", "T"}
        notes = []
        if base not in ok:
            notes.append(f"labelled base {base} cannot be a {mod} site on either strand"); bad += 1
        if rc(motif) == motif and strand == "+":
            notes.append(f"palindromic but '+' only: offset {len(motif)-1-off} ({motif[len(motif)-1-off]}) is never labelled"); bad += 1
        print(f"   {motif:12s} offset={off} base={base} strand={strand:4s} declared={mod or '?':4s} {'  <-- ' + '; '.join(notes) if notes else ''}")
print(f"\n{bad} inconsistent entr{'y' if bad == 1 else 'ies'}")
