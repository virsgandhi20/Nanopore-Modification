#!/usr/bin/env python3
"""Write per-read modification labels for one base into a BAM as SAM MM/ML tags, from a list of modified reference sites.

Every occurrence of --base in a read (in the read's own orientation) gets an ML value: 255 when the sample is
--state modified AND the base sits on a reference position listed in --sites (plus-strand coordinates; a base on a
minus-strand read is looked up at its own plus coordinate, so list both strands' positions when both carry the
modification), 0 otherwise. So a control sample (--state unmodified) gets 0 everywhere, and in the modified sample the
bases outside the designed sites are labelled unmodified. --code is the SAM modification code (h = 5hmC, 21839 = 4mC).
UniMeth's finetuner reads these tags (pysam modified_bases_forward). Existing MM/ML entries are kept and appended to.
Secondary, supplementary and unmapped records are dropped; the output is coordinate-sorted and indexed.
"""
import argparse, os, sys
from array import array
import pysam
ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="inp", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--base", required=True, choices=["A", "C", "G", "T"]); ap.add_argument("--code", required=True, help="SAM mod code: h, m, a, g, 21839 ...")
ap.add_argument("--sites", required=True, help="BED of modified reference positions (contig, 0-based pos); read only when --state modified")
ap.add_argument("--state", choices=["modified", "unmodified"], required=True)
ap.add_argument("--limit", type=int, default=0, help="stop after N kept records (tests)")
a = ap.parse_args()
sites = set()
if a.state == "modified":
    for line in open(a.sites):
        c = line.split()
        if len(c) >= 2 and not c[0].startswith("#"): sites.add((c[0], int(c[1])))
kept = dropped = n_base = n_mod = 0
tmp = a.out + ".unsorted.bam"                                       # tagged in input order, then coordinate-sorted (the refined oligo BAMs are unsorted)
with pysam.AlignmentFile(a.inp, "rb", check_sq=False) as fin, pysam.AlignmentFile(tmp, "wb", template=fin) as fout:
    for r in fin:
        if r.is_unmapped or r.is_secondary or r.is_supplementary: dropped += 1; continue
        fseq = r.get_forward_sequence()
        if not fseq: dropped += 1; continue
        fseq = fseq.upper(); L = len(fseq)
        idx = [i for i, b in enumerate(fseq) if b == a.base]          # forward-read positions of the base
        if not idx: fout.write(r); kept += 1; continue
        ml = [0] * len(idx)
        if sites:
            q2r = {}
            for q, p in r.get_aligned_pairs(matches_only=True): q2r[q] = p
            ctg = r.reference_name
            for j, i in enumerate(idx):
                q = (L - 1 - i) if r.is_reverse else i                # forward index -> SEQ index
                p = q2r.get(q)
                if p is not None and (ctg, p) in sites: ml[j] = 255; n_mod += 1
        mm = f"{a.base}+{a.code}?," + ",".join(["0"] * len(idx)) + ";"
        mlarr = array("B", ml)
        if r.has_tag("MM") and r.get_tag("MM"):
            mm = r.get_tag("MM") + mm; mlarr = array("B", list(r.get_tag("ML"))) + mlarr
        r.set_tag("MM", mm, "Z"); r.set_tag("ML", mlarr)
        fout.write(r); kept += 1; n_base += len(idx)
        if a.limit and kept >= a.limit: break
if kept == 0: sys.exit("no records written")
pysam.sort("-@", "4", "-o", a.out, tmp); os.remove(tmp); pysam.index(a.out)
print(f"{a.out}: kept {kept:,} records ({dropped:,} dropped), {n_base:,} {a.base} positions tagged {a.base}+{a.code}, {n_mod:,} of them modified (ML=255); sorted and indexed")
