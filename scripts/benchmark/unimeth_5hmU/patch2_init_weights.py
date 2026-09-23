#!/usr/bin/env python3
"""Second, incremental patch for the 5hmU UniMeth clone (apply after patch_unimeth_5hmU.py; idempotent).

Adds UNIMETH_INIT_WEIGHTS=<checkpoint-N/pytorch_model.bin or model.safetensors>: before training
starts, load those weights into the model. This continues a run from a Trainer checkpoint WITHOUT
the Trainer's own resume path, which goes through transformers' torch.load wrapper and is refused
on torch < 2.6 (CVE-2025-32434). Only the weights carry over; the optimizer and the learning-rate
schedule restart, so pass the REMAINING number of steps as --max_steps.
Usage: patch2_init_weights.py <path to the patched Unimeth clone>
"""
import os, sys
root = sys.argv[1] if len(sys.argv) > 1 else "."
p = os.path.join(root, "unimeth", "training", "trainer_base.py")
s = open(p).read()
MARK = "# 5hmU patch 2: init weights from a checkpoint"
OLD_MARK = "# 5hmU patch 2: init weights from safetensors"
anchor = "        trainer.train(resume_from_checkpoint=True if _resume else None)\n"
assert s.count(anchor) == 1, "anchor not found: apply patch_unimeth_5hmU.py first"
if MARK in s:
    print("patch2 already applied:", p); sys.exit(0)
if OLD_MARK in s:                                     # strip the earlier variant
    i = s.index("        " + OLD_MARK); j = s.index(anchor); s = s[:i] + s[j:]
block = '''        %s (continue a run without the Trainer's torch.load-based resume)
        _init = os.environ.get('UNIMETH_INIT_WEIGHTS')
        if _init:
            if _init.endswith('.safetensors'):
                from safetensors.torch import load_file as _load_sf
                _sd = _load_sf(_init)
            else:
                import torch as _torch_init
                _sd = _torch_init.load(_init, map_location='cpu', weights_only=True)
            _model = trainer.accelerator.unwrap_model(trainer.model)
            _res = _model.load_state_dict(_sd, strict=False)
            _cur = _model.state_dict()
            _loaded_ptrs = {_cur[k].data_ptr() for k in _sd if k in _cur}
            _really_missing = [k for k in _res.missing_keys if _cur[k].data_ptr() not in _loaded_ptrs]   # tied weights may be absent from the file but share storage with a loaded one
            local_print(f'[5hmU] init weights from {_init}: {sum(1 for k in _sd if k in _cur)}/{len(_cur)} tensors loaded, tied={len(_res.missing_keys) - len(_really_missing)}, missing={_really_missing[:5]}, unexpected={_res.unexpected_keys[:5]}')
            if _really_missing or _res.unexpected_keys:
                raise RuntimeError('UNIMETH_INIT_WEIGHTS does not match the model')
''' % MARK
s = s.replace(anchor, block + anchor)
open(p, "w").write(s); print("patch2 applied:", p)
