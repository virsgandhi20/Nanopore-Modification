#!/usr/bin/env python3
"""Second label patch for the UniMeth clone that already carries the 5hmU patch (patch_unimeth_5hmU.py):
--hmC 1 reads C labels from C+h tags (5hmC), --m4C 1 from C+21839 tags (4mC). The C context tokens ([CpG], [CHG],
[CHH]) are reused unchanged, so no vocabulary change: a model fine-tuned this way is a 5hmC (or 4mC) caller at every
cytosine, and the published all-context 5mC checkpoint is the natural starting point. Idempotent, and it repairs the
indentation of an earlier application (the first version indented the new CLI lines wrongly in one file).
Usage: python patch_labels.py /path/to/Unimeth_5hmU"""
import os, re, sys
root = sys.argv[1]; MARK = "# label patch (5hmC/4mC)"
def edit(rel, pairs):
    p = os.path.join(root, rel); s = open(p).read()
    if MARK in s: print("already patched", rel); return
    for old, new, n in pairs:
        assert s.count(old) == n, f"{rel}: expected {n} of {old!r}, found {s.count(old)}"
        s = s.replace(old, new)
    open(p, "w").write(s); print("patched", rel)
edit("unimeth/data/extract.py", [
    ("        if self.detect_hmu == 1:\n            self.detect_mod = ('T', 0, 'g')\n",
     f"        {MARK}\n        if getattr(args, 'hmC', 0) == 1:\n            self.detect_mod = ('C', 0, 'h')\n        elif getattr(args, 'm4C', 0) == 1:\n            self.detect_mod = ('C', 0, 21839)\n        elif self.detect_hmu == 1:\n            self.detect_mod = ('T', 0, 'g')\n", 1),
])
edit("unimeth/utils/bam_tags.py", [
    ("    'g': ('T', ('T', 0, 'g')),\n}", f"    'g': ('T', ('T', 0, 'g')),\n    {MARK}\n    'h': ('C', ('C', 0, 'h')),\n    '21839': ('C', ('C', 0, 21839)),\n    21839: ('C', ('C', 0, 21839)),\n}}", 1),
])
# CLI flags: inserted right after the --hmU line, with THAT line's indentation (4 spaces in training/__main__.py,
# 8 in args_config.py); an earlier application with the wrong indent is re-indented here.
NEW = ["parser.add_argument('--hmC', type=int, default=0, help='Labels from C+h tags: 5hmC at the C context sites (1=yes)')",
       "parser.add_argument('--m4C', type=int, default=0, help='Labels from C+21839 tags: 4mC at the C context sites (1=yes)')"]
for rel in ("unimeth/config/args_config.py", "unimeth/training/__main__.py"):
    p = os.path.join(root, rel); lines = open(p).read().split("\n")
    i_hmu = next((i for i, l in enumerate(lines) if "add_argument('--hmU'" in l), None)
    assert i_hmu is not None, f"{rel}: --hmU line not found (apply patch_unimeth_5hmU.py first)"
    indent = re.match(r"\s*", lines[i_hmu]).group(0)
    block = [indent + MARK] + [indent + n for n in NEW]
    if any(MARK in l for l in lines):
        j = next(i for i, l in enumerate(lines) if MARK in l)
        if lines[j:j + 3] == block: print("already patched", rel); continue
        lines[j:j + 3] = block; print("re-indented", rel)
    else:
        lines[i_hmu + 1:i_hmu + 1] = block; print("patched", rel)
    open(p, "w").write("\n".join(lines))
    compile(open(p).read(), p, "exec")                      # fail loudly here, not at job time
print("label patch applied")
