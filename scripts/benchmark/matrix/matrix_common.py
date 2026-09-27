"""Shared parsing of datasets.tsv and rows.tsv (paths expanded with $CAT, $EUK, $UMBC, $ME from the environment)."""
import os, re
DEFAULTS = {"CAT": "/fs/cbcb-lab/storm/shared/data", "EUK": "/fs/cbcb-lab/storm/vgandhi/euk",
            "UMBC": "/fs/cbcb-lab/storm/shared/umbc-ont-data", "ME": "/fs/nexus-scratch/vgandhi"}
def expand(s):
    return re.sub(r"\$\{(\w+)\}", lambda m: os.environ.get(m.group(1), DEFAULTS.get(m.group(1), m.group(0))), s)
def _read(path, ncol):
    for line in open(path):
        if not line.strip() or line.startswith("#"): continue
        c = line.rstrip("\n").split("\t")
        if len(c) < ncol: raise SystemExit(f"{path}: line needs {ncol} tab-separated columns: {line.strip()[:80]}")
        yield [x.strip() for x in c]
def load_datasets(path):
    ds = {}
    for c in _read(path, 7):
        ds[c[0]] = {"id": c[0], "gtdir": expand(c[1]), "pod5": expand(c[2]), "bam": expand(c[3]), "ref": expand(c[4]),
                    "nreads": int(c[5]), "reuse": c[6]}
    return ds
def load_rows(path, datasets=None):
    rows = []
    for c in _read(path, 9):
        r = {"row": c[0], "chem": c[1], "base": c[2], "type": c[3], "datasets": c[4].split(","), "context": c[5],
             "mincov": int(c[8]), "note": c[9] if len(c) > 9 else ""}
        gtdir = datasets[r["datasets"][0]]["gtdir"] if datasets and r["datasets"][0] in datasets else None
        def paths(spec):
            out = []
            for p in spec.split(","):
                p = expand(p.strip())
                if p == "same" or p.startswith("refbase:") or os.path.isabs(p) or gtdir is None: out.append(p)
                else: out.append(os.path.join(gtdir, p))
            return out
        r["gt"] = paths(c[6]); r["cand"] = paths(c[7]); rows.append(r)
    return rows
