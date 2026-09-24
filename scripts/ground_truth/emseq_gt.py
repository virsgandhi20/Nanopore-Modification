#!/usr/bin/env python3
"""Per-site methylation ground truth from a Bismark (EM-seq / bisulfite) BAM, in one streaming pass.

Bismark stores, per read, the methylation call of every cytosine in the XM tag (Z/z = CpG
methylated/unmethylated, X/x = CHG, H/h = CHH, '.' = not a cytosine) and the strand the read came
from in XG (CT = original top strand, the C sits at the aligned position on '+'; GA = bottom strand,
the C is the complement of the G at that position, i.e. on '-'). The XM string follows the SEQ as
stored in the BAM. Mates are adjacent in Bismark order, so a pair's overlap is counted once (the
second mate's positions inside the first mate's span are skipped). No index is needed and nothing
intermediate is written; memory is 2 x uint16 per strand per base of the contigs kept
(region-limit mouse-sized genomes).

Output (prefix <out>):
  <out>.sites.tsv.gz       contig, pos (0-based plus coordinate of the C), strand, context, n_meth,
                           n_total, frac   for every site with n_total >= --min-cov
  <out>.gt_modified.bed    sites with frac >= --hi   (BED6, score = frac)
  <out>.gt_unmodified.bed  sites with frac <= --lo
  <out>.candidates.bed     the union: pass it to score_sites.py --candidates so the ambiguous band
                           between --lo and --hi is left out of the scoring
  <out>.summary.txt
CpG sites are merged across strands (the C at p on '+' with the C at p+1 on '-') and written as BOTH
coordinates with the merged counts, so a per-position scorer sees both; CHG and CHH stay stranded.
"""
import argparse, gzip, os, sys, time
import numpy as np
import pysam

CTX = {"CpG": (ord("Z"), ord("z")), "CHG": (ord("X"), ord("x")), "CHH": (ord("H"), ord("h"))}

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--bam", required=True, help="Bismark BAM (deduplicated, Bismark order; no index needed)")
ap.add_argument("--context", default="CpG", choices=["CpG", "CHG", "CHH", "nonCpG"], help="nonCpG = CHG and CHH")
ap.add_argument("--region", default=None, help="contig or contig:start-end (1-based, inclusive). Default: every contig")
ap.add_argument("--min-cov", type=int, default=10)
ap.add_argument("--hi", type=float, default=0.9, help="frac >= hi -> modified")
ap.add_argument("--lo", type=float, default=0.1, help="frac <= lo -> unmodified")
ap.add_argument("--min-mapq", type=int, default=20)
ap.add_argument("--max-reads", type=int, default=0, help="stop after this many records (tests)")
ap.add_argument("--out", required=True, help="output prefix")
a = ap.parse_args()
t0 = time.time()
contexts = ["CHG", "CHH"] if a.context == "nonCpG" else [a.context]
os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)

bam = pysam.AlignmentFile(a.bam, "rb", check_sq=False)
lengths = dict(zip(bam.references, bam.lengths))
keep, rstart, rend = None, None, None
if a.region:
    ctg, _, span = a.region.partition(":")
    if ctg not in lengths:
        sys.exit(f"contig {ctg!r} is not in the BAM header; contigs start with {bam.references[:5]}")
    keep = {ctg}
    if span:
        s, _, e = span.partition("-"); rstart, rend = int(s) - 1, int(e)          # to 0-based half-open
    else:
        rstart, rend = 0, lengths[ctg]

# counts[contig][context] = {"meth": [plus, minus], "tot": [plus, minus]}   (uint16 arrays, allocated on first use)
counts = {}
def ensure(ctg):
    c = counts.get(ctg)
    if c is None:
        L = lengths[ctg]
        c = {ctx: {"meth": [np.zeros(L, np.uint16), np.zeros(L, np.uint16)],
                   "tot": [np.zeros(L, np.uint16), np.zeros(L, np.uint16)]} for ctx in contexts}
        counts[ctg] = c
    return c

n_rec = n_used = n_skipped_qc = n_overlap_trim = n_indel = 0
prev_name, prev_span = None, None
for r in bam:
    n_rec += 1
    if a.max_reads and n_rec > a.max_reads:
        break
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.is_qcfail or r.is_duplicate or r.mapping_quality < a.min_mapq:
        n_skipped_qc += 1; continue
    ctg = r.reference_name
    if keep is not None and ctg not in keep:
        continue
    try:
        xm = r.get_tag("XM"); xg = r.get_tag("XG")
    except KeyError:
        n_skipped_qc += 1; continue
    xmb = np.frombuffer(xm.encode("ascii"), dtype=np.uint8)
    cig = r.cigartuples
    if cig is not None and len(cig) == 1 and cig[0][0] == 0 and cig[0][1] == len(xmb):     # the common case: all M
        qidx = np.arange(len(xmb)); rpos = qidx + r.reference_start
    else:                                                                                   # indels or clips: exact pairs
        n_indel += 1
        pairs = r.get_aligned_pairs(matches_only=True)
        if not pairs:
            continue
        qidx = np.fromiter((q for q, _ in pairs), dtype=np.int64, count=len(pairs))
        rpos = np.fromiter((p for _, p in pairs), dtype=np.int64, count=len(pairs))
    if prev_name == r.query_name and prev_span is not None:                                  # second mate: count the overlap once
        m = (rpos < prev_span[0]) | (rpos >= prev_span[1])
        if not m.all():
            n_overlap_trim += 1; qidx, rpos = qidx[m], rpos[m]
    prev_name, prev_span = r.query_name, (r.reference_start, r.reference_end)
    if rstart is not None:
        m = (rpos >= rstart) & (rpos < rend)
        qidx, rpos = qidx[m], rpos[m]
    if len(rpos) == 0:
        continue
    codes = xmb[qidx]; s = 1 if xg == "GA" else 0
    c = ensure(ctg); used = False
    for ctx in contexts:
        up, lo = CTX[ctx]
        pm = rpos[codes == up]; pu = rpos[codes == lo]
        if len(pm):
            np.add.at(c[ctx]["meth"][s], pm, 1); np.add.at(c[ctx]["tot"][s], pm, 1); used = True
        if len(pu):
            np.add.at(c[ctx]["tot"][s], pu, 1); used = True
    n_used += used
    if n_rec % 2_000_000 == 0:
        print(f"  {n_rec:,} records, {n_used:,} used, {(time.time() - t0) / 60:.1f} min", file=sys.stderr, flush=True)

# ------------------------------------------------------------------ sites and thresholds
n_sites = n_mod = n_unmod = 0
with gzip.open(a.out + ".sites.tsv.gz", "wt") as fs, open(a.out + ".gt_modified.bed", "w") as fm, \
     open(a.out + ".gt_unmodified.bed", "w") as fu, open(a.out + ".candidates.bed", "w") as fc:
    fs.write("contig\tpos\tstrand\tcontext\tn_meth\tn_total\tfrac\n")
    def emit(ctg, pos, strand, ctx, nm, nt):
        global n_sites, n_mod, n_unmod
        frac = nm / nt
        fs.write(f"{ctg}\t{pos}\t{strand}\t{ctx}\t{nm}\t{nt}\t{frac:.4f}\n"); n_sites += 1
        line = f"{ctg}\t{pos}\t{pos + 1}\t{ctx}\t{frac:.4f}\t{strand}\n"
        if frac >= a.hi:
            fm.write(line); fc.write(line); n_mod += 1
        elif frac <= a.lo:
            fu.write(line); fc.write(line); n_unmod += 1
    for ctg in sorted(counts, key=lambda k: bam.references.index(k)):
        for ctx in contexts:
            M, T = counts[ctg][ctx]["meth"], counts[ctg][ctx]["tot"]
            if ctx == "CpG":                                   # merge + at p with - at p+1
                nm = M[0][:-1].astype(np.int64) + M[1][1:]; nt = T[0][:-1].astype(np.int64) + T[1][1:]
                for p in np.flatnonzero(nt >= a.min_cov):
                    emit(ctg, int(p), "+", ctx, int(nm[p]), int(nt[p])); emit(ctg, int(p) + 1, "-", ctx, int(nm[p]), int(nt[p]))
            else:
                for s, strand in ((0, "+"), (1, "-")):
                    nt = T[s]; nm = M[s]
                    for p in np.flatnonzero(nt >= a.min_cov):
                        emit(ctg, int(p), strand, ctx, int(nm[p]), int(nt[p]))

summary = (f"emseq_gt {a.bam}\n  context={a.context} region={a.region or 'all'} min_cov={a.min_cov} hi={a.hi} lo={a.lo} min_mapq={a.min_mapq}\n"
           f"  records={n_rec:,} used={n_used:,} qc_skipped={n_skipped_qc:,} mate_overlap_trimmed={n_overlap_trim:,} indel_or_clip_path={n_indel:,}\n"
           f"  sites_with_cov>={a.min_cov}: {n_sites:,}  modified(>= {a.hi}): {n_mod:,}  unmodified(<= {a.lo}): {n_unmod:,}  ambiguous: {n_sites - n_mod - n_unmod:,}\n"
           f"  minutes={(time.time() - t0) / 60:.1f}\n")
open(a.out + ".summary.txt", "w").write(summary); print(summary, end="")
