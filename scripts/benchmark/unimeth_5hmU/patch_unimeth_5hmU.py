#!/usr/bin/env python3
"""Teach a UniMeth clone (v0.3.1, commit 11215d4) a fifth modification: 5hmU on T.

Applied to a SEPARATE clone; the benchmark install is untouched. Every edit is an
exact-string replacement with an assertion, so an upstream change fails loudly
instead of silently producing a half-patched package.

What changes
  vocab           '[5hmU]' appended (id 17); existing ids unchanged
  site detection  T positions are candidate sites when --hmU 1
  labels          read from BAM tag T+g (SAM code g = 5-hydroxymethyluracil)
  checkpoint load a 17-token checkpoint is expanded to 18 rows; the new row is
                  initialised from the [m6A] row (closest task: a per-base,
                  non-CpG modification)
  trainer         final weights saved as <out>/final.pt (upstream only saves
                  every 5000 steps, so short runs saved nothing); a few
                  hyperparameters made overridable through environment variables
Usage: python patch_unimeth_5hmU.py /path/to/Unimeth_5hmU
"""
import os, sys

root = sys.argv[1]
def edit(rel, pairs):
    p = os.path.join(root, rel); s = open(p).read()
    for old, new, n in pairs:
        assert s.count(old) == n, f"{rel}: expected {n} occurrence(s) of {old!r}, found {s.count(old)}"
        s = s.replace(old, new)
    open(p, "w").write(s); print("patched", rel)

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
edit("unimeth/config/args_config.py", [
    ("        parser.add_argument('--m6A', type=int, default=0, help='Enable m6A detection (1=yes)')\n",
     "        parser.add_argument('--m6A', type=int, default=0, help='Enable m6A detection (1=yes)')\n"
     "        parser.add_argument('--hmU', type=int, default=0, help='Enable 5hmU detection at every T (1=yes)')\n", 1),
])
edit("unimeth/training/__main__.py", [
    ("    parser.add_argument('--m6A', type=int, default=0, help='Enable m6A detection (1=yes)')\n",
     "    parser.add_argument('--m6A', type=int, default=0, help='Enable m6A detection (1=yes)')\n"
     "    parser.add_argument('--hmU', type=int, default=0, help='Enable 5hmU detection at every T (1=yes)')\n", 1),
])
edit("unimeth/eval/metrics.py", [
    ("        methy_types = ['[CpG]', '[CHG]', '[CHH]', '[m6A]']", "        methy_types = ['[CpG]', '[CHG]', '[CHH]', '[m6A]', '[5hmU]']", 1),
])
edit("unimeth/utils/bam_tags.py", [
    ("    'a': ('A', ('A', 0, 'a')),\n}", "    'a': ('A', ('A', 0, 'a')),\n    '5hmU': ('T', ('T', 0, 'g')),\n    'g': ('T', ('T', 0, 'g')),\n}", 1),
])
edit("unimeth/model/loader.py", [
    ("""    state_dict = torch.load(checkpoint_path, map_location='cpu')
    model.load_state_dict(state_dict, strict=True)
    return model""",
     """    state_dict = torch.load(checkpoint_path, map_location='cpu')
    # 5hmU patch: a checkpoint written with a smaller vocabulary is expanded along
    # its vocabulary dimension; new rows start as a copy of the [m6A] row.
    from unimeth.config import TOKENIZER
    own = model.state_dict(); src_row = TOKENIZER['[m6A]']; grown = []
    for k, v in list(state_dict.items()):
        if k not in own or own[k].shape == v.shape or own[k].ndim != v.ndim:
            continue
        diff = [d for d in range(v.ndim) if own[k].shape[d] != v.shape[d]]
        if len(diff) != 1 or own[k].shape[diff[0]] < v.shape[diff[0]]:
            continue
        d = diff[0]; extra = own[k].shape[d] - v.shape[d]
        fill = v.narrow(d, src_row, 1).repeat_interleave(extra, dim=d) if v.shape[d] > src_row else torch.zeros_like(v.narrow(d, 0, 1)).repeat_interleave(extra, dim=d)
        state_dict[k] = torch.cat([v, fill.to(v.dtype)], dim=d); grown.append(f"{k} {tuple(v.shape)}->{tuple(state_dict[k].shape)}")
    if grown:
        print("[5hmU patch] expanded checkpoint tensors: " + "; ".join(grown), flush=True)
    model.load_state_dict(state_dict, strict=True)
    return model""", 1),
])
edit("unimeth/training/trainer_base.py", [
    ("            bf16=True,", "            bf16=os.environ.get('UNIMETH_BF16', '1') == '1',", 1),
    ("            dataloader_num_workers=4,", "            dataloader_num_workers=int(os.environ.get('UNIMETH_DL_WORKERS', '4')),", 1),
    ("            dataloader_prefetch_factor=2,", "            dataloader_prefetch_factor=2 if int(os.environ.get('UNIMETH_DL_WORKERS', '4')) > 0 else None,", 1),
    ('            output_dir=f"models/{self.mode}/{self.args.run_name}",',
     '            output_dir=os.environ.get("UNIMETH_OUT_DIR") or f"models/{self.mode}/{self.args.run_name}",', 1),
    ("            logging_steps=500,", "            logging_steps=int(os.environ.get('UNIMETH_LOG_STEPS', '500')),", 1),
    ("        return steps.get(self.mode, 5000)\n", "        return int(os.environ.get('UNIMETH_SAVE_STEPS', steps.get(self.mode, 5000)))\n", 1),
    ("        return steps.get(self.mode, 10000)\n", "        return int(os.environ.get('UNIMETH_EVAL_STEPS', steps.get(self.mode, 10000)))\n", 1),
    ("""        if self.mode != 'pretrain':
            local_print(trainer.evaluate())

        trainer.train()""",
     """        if self.mode != 'pretrain' and os.environ.get('UNIMETH_SKIP_INITIAL_EVAL', '0') != '1':
            local_print(trainer.evaluate())
        # 5hmU patch: resume from the last checkpoint in the output dir (preempted jobs are requeued)
        _out = self.training_args.output_dir
        _resume = os.environ.get('UNIMETH_RESUME', '0') == '1' and os.path.isdir(_out) and any(d.startswith('checkpoint-') for d in os.listdir(_out))
        trainer.train(resume_from_checkpoint=True if _resume else None)
        # 5hmU patch: always leave the final weights behind, in the plain state-dict
        # form that load_model() reads (upstream only saves every save_steps).
        import torch as _torch
        final = os.path.join(self.training_args.output_dir, 'final.pt')
        _torch.save(trainer.accelerator.unwrap_model(trainer.model).state_dict(), final)
        local_print(f'saved final weights to {final}')
        if self.mode != 'pretrain':
            local_print(trainer.evaluate())""", 1),
])
print("all patches applied")
