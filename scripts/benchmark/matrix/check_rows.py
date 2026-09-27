#!/usr/bin/env python3
"""Login-node check of datasets.tsv and rows.tsv before anything is submitted: every path exists, the BAMs carry move
tables, and the ground-truth / candidate positions sit on the base the row says they do (catches 1-based files and
wrong-strand splits). Prints a directory listing for every missing file so the tables can be corrected from one paste.
The reference is indexed in the work dir (the collection's folders are read-only), which prep then reuses."""
import argparse, os, shutil, sys, traceback
from matrix_common import load_datasets, load_rows
ap = argparse.ArgumentParser(); ap.add_argument("--datasets", required=True); ap.add_argument("--rows", required=True)
ap.add_argument("--work", required=True); ap.add_argument("--euk", default=os.environ.get("EUK", "/fs/cbcb-lab/storm/vgandhi/euk"))
ap.add_argument("--sample", type=int, default=5000)
ap.add_argument("--detail", default="", help="row id: per-file base composition and sequence context of off-base positions")
a = ap.parse_args()
try:
    import pysam
except ImportError:
    pysam = None; print("pysam not importable here: base and move-table checks skipped (activate the unimeth env)")
ds = load_datasets(a.datasets); rows = load_rows(a.rows, ds); bad = 0; listed = set()
def problem(msg):
    global bad
    bad += 1; print("  " + msg)
def missing(p, what):
    problem(f"MISSING {what}: {p}")
    d = os.path.dirname(p)
    while d and not os.path.isdir(d): d = os.path.dirname(d)
    if d and d not in listed:
        listed.add(d); print(f"    contents of {d}:")
        try:
            for e in sorted(os.listdir(d))[:60]: print("      " + e + ("/" if os.path.isdir(os.path.join(d, e)) else ""))
        except OSError as e: print("      (cannot list: %s)" % e)
def indexed_ref(d):
    """path of an indexed copy of the dataset's reference: work dir, else the euk folder, else make one now"""
    W = os.path.join(a.work, d["id"]); os.makedirs(W, exist_ok=True); w = os.path.join(W, "ref.fa")
    if os.path.exists(w + ".fai"): return w
    if d["reuse"] != "-" and os.path.exists(os.path.join(a.euk, d["reuse"], "ref.fa.fai")):
        for suf in ("", ".fai"):
            os.path.lexists(w + suf) and os.remove(w + suf); os.symlink(os.path.join(a.euk, d["reuse"], "ref.fa" + suf), w + suf)
        return w
    if os.path.exists(d["ref"] + ".fai"):
        for suf in ("", ".fai"):
            os.path.lexists(w + suf) and os.remove(w + suf); os.symlink(os.path.realpath(d["ref"] + suf), w + suf)
        return w
    print(f"    indexing a copy of the reference under {W} (read-only source folder)")
    shutil.copyfile(os.path.realpath(d["ref"]), w); pysam.faidx(w); return w
print("== datasets")
for k, d in ds.items():
    print(f"{k}: nreads={d['nreads']} reuse={d['reuse']}")
    ok = True
    for what in ("pod5", "bam", "ref"):
        if not os.path.exists(d[what]): missing(d[what], what); ok = False
    if d["gtdir"] != "-" and not os.path.isdir(d["gtdir"]): missing(d["gtdir"], "gtdir")
    if ok and pysam:
        try:
            with pysam.AlignmentFile(d["bam"], "rb", check_sq=False) as bam:
                n = mv = 0
                for r in bam.fetch(until_eof=True):
                    n += 1; mv += r.has_tag("mv")
                    if n >= 200: break
                so = bam.header.get("HD", {}).get("SO", "?")
            print(f"  bam: {mv} of {n} first records carry move tables, sort order {so}" + ("" if mv else "   <-- no mv tags: prep will re-basecall from the pod5 (MODE=basecall)"))
        except Exception as e:
            problem(f"bam unreadable: {e}")
print("== rows")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
for r in rows:
    print(f"{r['row']}: {r['chem']} on {r['base']}, {r['type']}-sample {r['datasets']}, context {r['context']}, mincov {r['mincov']}")
    try:
        for dsid in r["datasets"]:
            if dsid not in ds: problem(f"UNKNOWN dataset {dsid}")
        pos_ds = ds.get(r["datasets"][0]); files = {"gt": [], "cand": []}
        for col in ("gt", "cand"):
            for p in r[col]:
                if p == "same" or p.startswith("refbase:"): continue
                if not os.path.exists(p): missing(p, col)
                else: files[col].append(p)
        if pysam and pos_ds and r["base"] != "N" and os.path.exists(pos_ds["ref"]):
            fa = pysam.FastaFile(indexed_ref(pos_ds)); refs = set(fa.references); cache = {}
            if a.detail == r["row"]:
                for col in ("gt", "cand"):
                    for p in files[col]:
                        comp = {}; off = []; n = 0
                        for line in open(p):
                            c = line.split()
                            if len(c) < 2 or c[0].startswith("#") or c[0] not in refs: continue
                            s_ = cache.get(c[0])
                            if s_ is None: s_ = cache[c[0]] = fa.fetch(c[0]).upper()
                            q = int(c[1]); n += 1; b = s_[q] if 0 <= q < len(s_) else "?"; comp[b] = comp.get(b, 0) + 1
                            if b not in (r["base"], COMP[r["base"]]) and len(off) < 12: off.append((c[0], q, s_[max(0, q - 5):q] + "[" + b + "]" + s_[q + 1:q + 6], line.strip()))
                        print(f"  DETAIL {col} {os.path.basename(p)}: {n:,} positions, bases {dict(sorted(comp.items()))}")
                        for ctg, q, ctx, raw in off: print(f"      {ctg}:{q}  {ctx}   line: {raw[:60]}")
            for col in ("gt", "cand"):
                if not files[col]: continue
                n = ok = ncol = unknown = 0
                for p in files[col]:
                    for line in open(p):
                        c = line.split()
                        if len(c) < 2 or c[0].startswith("#"): continue
                        ncol = max(ncol, len(c))
                        if c[0] not in refs: unknown += 1; continue
                        s = cache.get(c[0])
                        if s is None: s = cache[c[0]] = fa.fetch(c[0]).upper()
                        q = int(c[1]); n += 1
                        if 0 <= q < len(s) and s[q] in (r["base"], COMP[r["base"]]): ok += 1
                        if n >= a.sample: break
                    if n >= a.sample: break
                frac = ok / n if n else 0
                flag = "" if frac >= 0.95 else "   <-- CHECK: positions are not on the expected base (1-based file? wrong strand?)"
                print(f"  {col} base check: {ok:,} of {n:,} sampled positions are {r['base']}/{COMP[r['base']]} ({frac:.1%}), {ncol} columns" + (f", {unknown} on contigs not in the reference" if unknown else "") + flag)
                if frac < 0.95: bad += 1
            if files["gt"] and files["cand"]:                       # positives outside the candidate set are never scored
                def full(paths):
                    S = set()
                    for p in paths:
                        for line in open(p):
                            c = line.split()
                            if len(c) >= 2 and not c[0].startswith("#"): S.add((c[0], int(c[1])))
                    return S
                G = full(files["gt"]); Cd = full(files["cand"]); out = G - Cd
                print(f"  gt vs cand: {len(G):,} positives, {len(Cd):,} candidates, {len(out):,} positives NOT in the candidate set" + ("   <-- CHECK: those positives cannot be scored by any tool" if out else ""))
                if out:
                    comp = {}
                    for ctg, q in list(out)[:20000]:
                        s_ = cache.get(ctg)
                        if s_ is None and ctg in refs: s_ = cache[ctg] = fa.fetch(ctg).upper()
                        b = s_[q] if s_ and 0 <= q < len(s_) else "?"; comp[b] = comp.get(b, 0) + 1
                    print(f"    bases of the positives outside the candidates: {dict(sorted(comp.items()))}"); bad += 1
    except Exception:
        problem("check crashed: " + traceback.format_exc().strip().splitlines()[-1])
print(f"== {bad} problem(s)" if bad else "== all paths present, ground truth on the expected bases")
sys.exit(1 if bad else 0)
