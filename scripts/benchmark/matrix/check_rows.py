#!/usr/bin/env python3
"""Login-node check of datasets.tsv and rows.tsv before anything is submitted: every path exists, and the ground-truth
positions sit on the base the row says they do (catches 1-based files and wrong-strand splits). Prints a directory
listing for every missing file so the table can be corrected from one paste."""
import argparse, os, sys
from matrix_common import load_datasets, load_rows, expand
ap = argparse.ArgumentParser(); ap.add_argument("--datasets", required=True); ap.add_argument("--rows", required=True)
ap.add_argument("--sample", type=int, default=5000); a = ap.parse_args()
try:
    import pysam
except ImportError:
    pysam = None; print("pysam not importable here: base check skipped (activate the unimeth env)")
ds = load_datasets(a.datasets); rows = load_rows(a.rows, ds); bad = 0; listed = set()
def missing(p, what):
    global bad
    bad += 1; print(f"  MISSING {what}: {p}")
    d = os.path.dirname(p)
    while d and not os.path.isdir(d): d = os.path.dirname(d)
    if d and d not in listed:
        listed.add(d); print(f"    contents of {d}:")
        try:
            for e in sorted(os.listdir(d))[:60]: print("      " + e + ("/" if os.path.isdir(os.path.join(d, e)) else ""))
        except OSError as e: print("      (cannot list: %s)" % e)
print("== datasets")
for k, d in ds.items():
    print(f"{k}: nreads={d['nreads']} reuse={d['reuse']}")
    for what in ("pod5", "bam", "ref"):
        os.path.exists(d[what]) or missing(d[what], what)
    if d["gtdir"] != "-" and not os.path.isdir(d["gtdir"]): missing(d["gtdir"], "gtdir")
print("== rows")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
for r in rows:
    print(f"{r['row']}: {r['chem']} on {r['base']}, {r['type']}-sample {r['datasets']}, context {r['context']}, mincov {r['mincov']}")
    for dsid in r["datasets"]:
        if dsid not in ds: bad += 1; print(f"  UNKNOWN dataset {dsid}"); continue
    pos_ds = ds.get(r["datasets"][0]); files = []
    for col in ("gt", "cand"):
        for p in r[col]:
            if p in ("same", ) or p.startswith("refbase:"): continue
            if not os.path.exists(p): missing(p, col)
            elif col == "gt": files.append(p)
    if pysam and pos_ds and files and r["base"] != "N" and os.path.exists(pos_ds["ref"]):
        fa = pysam.FastaFile(pos_ds["ref"]); refs = set(fa.references); n = ok = notc = 0; cache = {}
        for p in files:
            for line in open(p):
                c = line.split()
                if len(c) < 2 or c[0].startswith("#"): continue
                if c[0] not in refs: continue
                s = cache.get(c[0])
                if s is None: s = cache[c[0]] = fa.fetch(c[0]).upper()
                q = int(c[1]); n += 1
                if 0 <= q < len(s) and s[q] in (r["base"], COMP[r["base"]]): ok += 1
                if n >= a.sample: break
            if n >= a.sample: break
        frac = ok / n if n else 0
        flag = "" if frac >= 0.95 else "   <-- CHECK: positions are not on the expected base (1-based file? wrong strand?)"
        print(f"  gt base check: {ok:,} of {n:,} sampled positions are {r['base']}/{COMP[r['base']]} ({frac:.1%}){flag}")
        if frac < 0.95: bad += 1
print(f"== {bad} problem(s)" if bad else "== all paths present, ground truth on the expected bases")
sys.exit(1 if bad else 0)
