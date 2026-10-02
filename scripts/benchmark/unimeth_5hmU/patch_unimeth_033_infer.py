#!/usr/bin/env python3
"""Inference-side port of the 5hmU patch to UniMeth v0.3.3 (the release with the authors' normalization fix, issue 21).

The three fine-tuned checkpoints (5hmU, 5hmC, 4mC) were trained in a v0.3.1 clone whose vocabulary carries an 18th
token, '[5hmU]'. To run them under v0.3.3 the inference code needs the same vocabulary, T as a candidate base when
--hmU 1, and the flag itself. Nothing on the training side is touched: fine-tuning stays on the v0.3.1 clone
(patch_unimeth_5hmU.py + patch_labels.py), inference of every UniMeth model moves to this clone.
Every edit is an exact-string replacement with an assertion. Idempotent.
Usage: python patch_unimeth_033_infer.py /path/to/Unimeth_0.3.3_clone"""
import os, sys
root = sys.argv[1]
def edit(rel, pairs):
    p = os.path.join(root, rel); s = open(p).read()
    if all(new in s for _, new, _ in pairs): print("already patched", rel); return
    for old, new, n in pairs:
        assert s.count(old) == n, f"{rel}: expected {n} occurrence(s) of {old!r}, found {s.count(old)}"
        s = s.replace(old, new)
    open(p, "w").write(s); compile(s, p, "exec"); print("patched", rel)
edit("unimeth/config/model_config.py", [
    ("'[R10]', '[4khz]', '[5khz]']", "'[R10]', '[4khz]', '[5khz]', '[5hmU]']", 1),
    ("    '[m6A]': TOKENIZER['[m6A]'],\n}", "    '[m6A]': TOKENIZER['[m6A]'],\n    '[5hmU]': TOKENIZER['[5hmU]'],\n}", 1),
    ("    m6A: int = 0\n", "    m6A: int = 0\n    hmU: int = 0\n", 1),
])
edit("unimeth/data/sites.py", [
    ("                           detect_chh: int, detect_m6a: int) -> list:",
     "                           detect_chh: int, detect_m6a: int, detect_hmu: int = 0) -> list:", 1),
    ("        elif seq[i] == 'A':\n            if detect_m6a:\n                pred_pos.append(i)\n",
     "        elif seq[i] == 'A':\n            if detect_m6a:\n                pred_pos.append(i)\n        elif seq[i] == 'T':\n            if detect_hmu:\n                pred_pos.append(i)\n", 1),
    ("                   detect_chh: int, detect_m6a: int) -> str | None:",
     "                   detect_chh: int, detect_m6a: int, detect_hmu: int = 0) -> str | None:", 1),
    ("    elif seq[pos] == 'A':\n        return '[m6A]' if detect_m6a else None\n",
     "    elif seq[pos] == 'A':\n        return '[m6A]' if detect_m6a else None\n    elif seq[pos] == 'T':\n        return '[5hmU]' if detect_hmu else None\n", 1),
])
edit("unimeth/data/extract.py", [
    ("        self.detect_m6a = getattr(args, 'm6A', 0)\n",
     "        self.detect_m6a = getattr(args, 'm6A', 0)\n        self.detect_hmu = getattr(args, 'hmU', 0)\n", 1),
    ("        if self.detect_m6a == 1:\n            self.detect_mod = ('A', 0, 'a')",
     "        if self.detect_hmu == 1:\n            self.detect_mod = ('T', 0, 'g')\n        elif self.detect_m6a == 1:\n            self.detect_mod = ('A', 0, 'a')", 1),
    ("seq, self.detect_cpg, self.detect_chg, self.detect_chh, self.detect_m6a\n",
     "seq, self.detect_cpg, self.detect_chg, self.detect_chh, self.detect_m6a, self.detect_hmu\n", 2),
])
edit("unimeth/data/pipeline.py", [
    ("               detect_cpg=0, detect_chg=0, detect_chh=0, detect_m6a=0):",
     "               detect_cpg=0, detect_chg=0, detect_chh=0, detect_m6a=0, detect_hmu=0):", 1),
    ("get_methy_type(bases, pos_p, detect_cpg, detect_chg, detect_chh, detect_m6a)",
     "get_methy_type(bases, pos_p, detect_cpg, detect_chg, detect_chh, detect_m6a, detect_hmu)", 1),
    ("        detect_m6a=args.m6A\n    )", "        detect_m6a=args.m6A,\n        detect_hmu=getattr(args, 'hmU', 0)\n    )", 2),
])
edit("unimeth/config/modification_names.py", [
    ('        ("6mA", "m6A", "m6A", "6mA detection"),\n    ):',
     '        ("6mA", "m6A", "m6A", "6mA detection"),\n        ("5hmU", "hmU", "hmU", "5hmU detection at every T (fine-tuned model)"),\n    ):', 1),
])
edit("unimeth/eval/metrics.py", [
    ("        methy_types = ['[CpG]', '[CHG]', '[CHH]', '[m6A]']", "        methy_types = ['[CpG]', '[CHG]', '[CHH]', '[m6A]', '[5hmU]']", 1),
])
print("v0.3.3 inference patch applied")
