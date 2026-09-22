#!/usr/bin/env python3
"""Write per-read 5hmU labels into a BAM as SAM modification tags (T+g, ML 255 or 0).

SPO1 phage DNA carries 5hmU in place of every thymine, and the PCR amplicon
libraries of the same genome carry none, so the label of every T in a read is
known from the sample alone: --state modified writes ML=255 for every T,
--state unmodified writes ML=0. UniMeth's finetuner reads exactly these tags
(pysam modified_bases_forward, key ('T', 0, 'g')). Existing MM/ML tags (e.g.
5mC/6mA calls) are kept and the T entries appended. Secondary, supplementary
and unmapped records are dropped; record order is preserved.
"""
import argparse, sys
from array import array
import pysam

ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--state", choices=["modified", "unmodified"], required=True)
ap.add_argument("--limit", type=int, default=0, help="stop after N kept records (tests)")
a = ap.parse_args()
ml_value = 255 if a.state == "modified" else 0
kept = dropped = 0; n_t = 0
with pysam.AlignmentFile(a.inp, "rb", check_sq=False) as fin, pysam.AlignmentFile(a.out, "wb", template=fin) as fout:
    for r in fin:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            dropped += 1; continue
        fseq = r.get_forward_sequence()
        if not fseq:
            dropped += 1; continue
        nT = fseq.upper().count("T")
        mm = "T+g?," + ",".join(["0"] * nT) + ";" if nT else ""
        ml = array("B", [ml_value] * nT)
        if r.has_tag("MM") and r.get_tag("MM"):
            mm = r.get_tag("MM") + mm
            old = r.get_tag("ML"); ml = array("B", list(old)) + ml
        if mm:
            r.set_tag("MM", mm, "Z"); r.set_tag("ML", ml)
        fout.write(r); kept += 1; n_t += nT
        if a.limit and kept >= a.limit: break
pysam.index(a.out)
print(f"{a.out}: kept {kept:,} records ({dropped:,} dropped), {n_t:,} T positions tagged {a.state} (ML={ml_value})")
if kept == 0: sys.exit("no records written")
