#!/usr/bin/env python3
"""Score every (row, tool) pair whose standardized sites exist and write the matrix.

For each row: resolve the ground truth (files, or every T for refbase rows), split by context if asked, build the
candidate set (two-sample rows: the same positions on the control sample under NEG_-prefixed contigs), keep only
candidates covered by >= mincov reads of the sample's sub.bam (samtools depth), then run score_sites.py with
--fill-missing twice per tool: score = call frequency, and score = mean P. Positions the tool never scored count as
unmodified, which is how a tool with no model for the row's chemistry lands at 0.5 instead of N/A.

Writes <work>/matrix_long.tsv (one line per row x tool x metric) and <work>/matrix_grid.tsv (rows x tools, the
tool's primary metric: mean P for UniMeth models, call frequency for the others).
"""
import argparse, os, subprocess, sys
from matrix_common import load_datasets, load_rows
ap = argparse.ArgumentParser()
ap.add_argument("--work", required=True); ap.add_argument("--datasets", required=True); ap.add_argument("--rows", required=True)
ap.add_argument("--tools", required=True, help="space-separated tool ids, e.g. 'unimeth_5mC unimeth_6mA deepmod2'")
ap.add_argument("--repo", required=True); ap.add_argument("--samtools", default="samtools")
ap.add_argument("--only", default="", help="comma-separated row ids to (re)score")
ap.add_argument("--force", action="store_true", help="rebuild cached candidate/site files")
a = ap.parse_args()
ds = load_datasets(a.datasets); rows = load_rows(a.rows, ds); tools = a.tools.split()
only = {x for x in a.only.split(",") if x}
scorer = os.path.join(a.repo, "scripts/benchmark/score_sites.py"); splitter = os.path.join(a.repo, "scripts/ground_truth/split_by_context.py")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def log(*x): print(*x, flush=True)
def read_positions(paths, refs=None):
    s = set()
    for p in paths:
        for line in open(p):
            c = line.split()
            if len(c) < 2 or c[0].startswith("#"): continue
            if refs is not None and c[0] not in refs: continue
            s.add((c[0], int(c[1])))
    return s
def write_bed(positions, path):
    with open(path, "w") as f:
        for c, p in sorted(positions): f.write(f"{c}\t{p}\t{p + 1}\n")
IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
         "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]"}
def revcomp(m):
    pairs = {"A": "T", "T": "A", "C": "G", "G": "C", "R": "Y", "Y": "R", "S": "S", "W": "W", "K": "M", "M": "K", "B": "V", "V": "B", "D": "H", "H": "D", "N": "N"}
    return "".join(pairs[c] for c in reversed(m))
def spec_positions(spec, ref):
    """refbase:<B>  every position whose reference base is B or its complement (a B on either strand)
       motif:<IUPAC>:<offset>:<+|both>  the base at <offset> of every motif match; both = the reverse complement too
       (the position then sits on the complementary base, like the collection's gt_minus files)"""
    import pysam, re
    fa = pysam.FastaFile(ref); out = set(); kind, _, rest = spec.partition(":")
    if kind == "refbase":
        want = {rest, COMP[rest]}
        for ctg in fa.references:
            seq = fa.fetch(ctg).upper(); out.update((ctg, i) for i, b in enumerate(seq) if b in want)
        return out
    motif, off, strands = rest.split(":"); off = int(off)
    pats = [(re.compile("(?=" + "".join(IUPAC[c] for c in motif) + ")"), off)]
    if strands == "both":
        pats.append((re.compile("(?=" + "".join(IUPAC[c] for c in revcomp(motif)) + ")"), len(motif) - 1 - off))
    for ctg in fa.references:
        seq = fa.fetch(ctg).upper()
        for pat, o in pats:
            out.update((ctg, m.start() + o) for m in pat.finditer(seq))
    return out
def covered(dsid, bed, mincov, out):
    """positions of `bed` covered by >= mincov reads of the sample's sub.bam (samtools depth -a -b)"""
    d = ds[dsid]; sub = os.path.join(a.work, dsid, "sub.bam")
    if not os.path.exists(sub): return None
    r = subprocess.run([a.samtools, "depth", "-a", "-b", bed, sub], capture_output=True, text=True)
    if r.returncode != 0: sys.exit(f"samtools depth failed on {sub}: {r.stderr.strip()[-300:]}")
    keep = set()
    for line in r.stdout.splitlines():
        c = line.split("\t")
        if len(c) >= 3 and int(c[2]) >= mincov: keep.add((c[0], int(c[1]) - 1))     # depth prints 1-based
    return keep
def parse_score(r):
    """last stdout line of score_sites.py -> dict; stderr carries the fill count"""
    out = r.stdout.strip().splitlines(); res = {"status": "failed", "detail": (r.stderr.strip() or r.stdout.strip())[-200:]}
    if r.returncode == 0 and out:
        c = out[-1].split("\t")
        if len(c) >= 7:
            res = {"status": "ok", "n_sites": c[1], "n_pos": c[2], "pos_rate": c[3], "auroc": c[4], "auprc": c[5], "f1": c[6], "detail": ""}
    for line in r.stderr.splitlines():
        if line.startswith("(fill-missing:"): res["filled"] = line.split(":")[1].split("of")[0].strip()
    return res

long_rows = []; grid = {}
for r in rows:
    if only and r["row"] not in only: continue
    R = os.path.join(a.work, "rows", r["row"]); os.makedirs(R, exist_ok=True)
    pos_ds = r["datasets"][0]; neg_ds = r["datasets"][1] if r["type"] == "two" else None
    ref = os.path.join(a.work, pos_ds, "ref.fa")
    if not os.path.exists(ref + ".fai"):
        log(f"[{r['row']}] prep not done for {pos_ds}, skipped"); continue
    import pysam; refs = set(pysam.FastaFile(ref).references)
    # ---- ground truth and candidates, once per row
    gt_bed, cand_bed = os.path.join(R, "gt.bed"), os.path.join(R, "cand_all.bed")
    if a.force or not os.path.exists(cand_bed):
        gt = spec_positions(r["gt"][0], ref) if r["gt"][0].startswith(("refbase:", "motif:")) else read_positions(r["gt"], refs)
        if r["cand"] == ["same"]: cand = set(gt)
        elif r["cand"][0].startswith(("refbase:", "motif:")): cand = spec_positions(r["cand"][0], ref)
        else: cand = read_positions(r["cand"], refs)
        if r["context"] in ("cpg", "noncpg"):
            for name, S in (("gt", gt), ("cand", cand)):
                tmp = os.path.join(R, f"{name}_in.bed"); write_bed(S, tmp)
                oc, on = os.path.join(R, f"{name}_cpg.bed"), os.path.join(R, f"{name}_noncpg.bed")
                sp = subprocess.run([sys.executable, splitter, "--ref", ref, "--bed", tmp, "--out-cpg", oc, "--out-noncpg", on], capture_output=True, text=True)
                if sp.returncode != 0: sys.exit(f"split_by_context failed for {r['row']}: {sp.stderr[-300:]}")
                S.clear(); S.update(read_positions([oc if r["context"] == "cpg" else on]))
        if r["cand"] != ["same"]:                                     # the scored set is the candidate set: positives = gt AND cand
            outside = len(gt - cand)
            if outside: log(f"[{r['row']}] {outside:,} of {len(gt):,} gt positions are not in the candidate file and are left out (e.g. M.SssI gt carries the GATC adenines)")
            gt &= cand
        write_bed(gt, gt_bed); write_bed(cand, cand_bed)
    gt = read_positions([gt_bed]); cand_pos = read_positions([cand_bed])
    # ---- coverage floor per sample (the tools could only have called what the reads cover)
    cov_bed = os.path.join(R, f"cand_cov{r['mincov']}.bed")
    if a.force or not os.path.exists(cov_bed):
        keep = covered(pos_ds, cand_bed, r["mincov"], cov_bed)
        if keep is None: log(f"[{r['row']}] no sub.bam for {pos_ds}, skipped"); continue
        cand = set(keep)
        if neg_ds:
            keep_n = covered(neg_ds, cand_bed, r["mincov"], cov_bed)
            if keep_n is None: log(f"[{r['row']}] no sub.bam for {neg_ds}, skipped"); continue
            cand |= {("NEG_" + c, p) for c, p in keep_n}
        write_bed(cand, cov_bed)
    n_cand = sum(1 for _ in open(cov_bed)); n_gt_cov = len(gt & read_positions([cov_bed]))
    log(f"[{r['row']}] gt {len(gt):,} positions, candidates covered >= {r['mincov']}: {n_cand:,} (positives among them {n_gt_cov:,})")
    # ---- one score per tool and metric
    for tool in tools:
        sites = os.path.join(R, f"{tool}.sites.tsv")
        src_pos = os.path.join(a.work, pos_ds, tool, "sites.std.tsv")
        src_neg = os.path.join(a.work, neg_ds, tool, "sites.std.tsv") if neg_ds else None
        if not os.path.exists(src_pos) or (src_neg and not os.path.exists(src_neg)):
            long_rows.append((r["row"], tool, "-", "pending", "", "", "", "", "", "", "sites missing")); grid[(r["row"], tool)] = "pending"; continue
        if a.force or not os.path.exists(sites) or os.path.getmtime(sites) < max(os.path.getmtime(src_pos), os.path.getmtime(src_neg) if src_neg else 0):
            with open(sites, "w") as fo:
                for line in open(src_pos): fo.write(line)
                if src_neg:
                    for line in open(src_neg): fo.write("NEG_" + line)
        primary = "mean_P" if tool.startswith("unimeth") else "call_freq"
        for metric, col in (("call_freq", "3"), ("mean_P", "4")):
            label = f"{r['row']}/{tool}/{metric}/cov{r['mincov']}"
            res = parse_score(subprocess.run([sys.executable, scorer, "--calls", sites, "--gt", gt_bed, "--candidates", cov_bed, "--min-cov", str(r["mincov"]),
                                              "--chrom-col", "0", "--pos-col", "1", "--cov-col", "2", "--freq-col", col, "--fill-missing", "--label", label],
                                             capture_output=True, text=True))
            long_rows.append((r["row"], tool, metric, res["status"], res.get("n_sites", ""), res.get("n_pos", ""), res.get("pos_rate", ""),
                              res.get("auroc", ""), res.get("auprc", ""), res.get("filled", ""), res.get("detail", "")))
            if metric == primary: grid[(r["row"], tool)] = res.get("auroc", res["status"])
            log(f"  {tool:14s} {metric:9s} {res['status']:7s} AUROC {res.get('auroc', '-'):7s} AUPRC {res.get('auprc', '-'):7s} n {res.get('n_sites', '-'):>9s} filled {res.get('filled', '-'):>8s} {res.get('detail', '')}")
# a partial run (--only, or a tool whose sites are missing) keeps the earlier results of the rows and cells it did not touch
long_path, grid_path = os.path.join(a.work, "matrix_long.tsv"), os.path.join(a.work, "matrix_grid.tsv")
scored_rows = {t[0] for t in long_rows}
if os.path.exists(long_path):
    for line in open(long_path):
        c = line.rstrip("\n").split("\t")
        if c[0] != "row" and len(c) >= 11 and c[0] not in scored_rows: long_rows.append(tuple(c[:11]))
if os.path.exists(grid_path):
    for line in open(grid_path):
        c = line.rstrip("\n").split("\t")
        if c[0] == "row": old_tools = c[2:]; continue
        for t, v in zip(old_tools, c[2:]):
            if v and (c[0], t) not in grid: grid[(c[0], t)] = v
order = {r["row"]: i for i, r in enumerate(rows)}
with open(long_path, "w") as f:
    f.write("row\ttool\tmetric\tstatus\tn_sites\tn_pos\tpos_rate\tauroc\tauprc\tn_filled\tdetail\n")
    for t in sorted(long_rows, key=lambda t: (order.get(t[0], 999), t[1], t[2])): f.write("\t".join(t) + "\n")
all_tools = list(tools) + sorted({t for (_, t) in grid if t not in tools})
with open(grid_path, "w") as f:
    f.write("row\tchem\t" + "\t".join(all_tools) + "\n")
    for r in rows:
        f.write(r["row"] + "\t" + r["chem"] + "\t" + "\t".join(grid.get((r["row"], t), "") for t in all_tools) + "\n")
log(f"wrote {a.work}/matrix_long.tsv and matrix_grid.tsv")
