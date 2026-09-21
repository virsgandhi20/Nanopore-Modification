#!/usr/bin/env python3
"""Is modification TYPE decodable from a frozen RawMod encoder?

RawMod (ConvFormerV2) encodes each read of a pileup image, lets reads attend to
each other, mean-pools to a 96-d site representation and outputs one detection
logit. This probe freezes a published checkpoint, reads out
  per-read embeddings BEFORE the cross-read Transformer (read_encoder output),
  per-read embeddings AFTER it (encoder output, the tap point for a typing head),
  the pooled site representation (what the detection head sees),
and fits a linear classifier for modification type on each, on the ONT oligo
samples where every image of a sample carries one known chemistry. A per-read
raw-signal baseline (mean current + dwell per base, 42-d) shows what the
embeddings add. Splits are by site (contig:pos), as in train_read_typing.py.

Images: positives (label 1) of ONT_{5mC,5hmC,6mA}_<strand>.h5, and ONT_control
images at the same positions as 'none'. Never raises on a missing file.
"""
import argparse, hashlib, json, os, sys, time
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import torch
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score
from sklearn.preprocessing import StandardScaler

ap = argparse.ArgumentParser()
ap.add_argument("--features", default="/fs/cbcb-lab/storm/bds062/rawmod_strand_resolved/features")
ap.add_argument("--checkpoints", required=True, help="name=path[,name=path...]")
ap.add_argument("--repo", default=str(Path(__file__).resolve().parents[2]))
ap.add_argument("--strand", default="plus"); ap.add_argument("--max-per-class", type=int, default=6000)
ap.add_argument("--out", required=True); ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--prefix", default="ONT_", help="file prefix (tests use a fake one)")
a = ap.parse_args()
t0 = time.time(); os.makedirs(a.out, exist_ok=True); rng = np.random.default_rng(a.seed)
for p in ("scripts/train", "scripts/test", "rawmod"): sys.path.insert(0, os.path.join(a.repo, p))
from model import PileupDataset                      # noqa: E402  (RawMod's own loader: adds the two delta channels)
from score_genome import load_model                  # noqa: E402  (RawMod's own checkpoint loader)
dev = "cuda" if torch.cuda.is_available() else "cpu"
result = {"args": vars(a), "device": dev, "checkpoints": {}, "warnings": []}
def finish(status):
    result["status"] = status; result["minutes"] = round((time.time() - t0) / 60, 1)
    json.dump(result, open(os.path.join(a.out, "metrics.json"), "w"), indent=1, default=str)
    open(os.path.join(a.out, "status.txt"), "w").write(f"rawmod_embedding_probe\t{status}\t{result['minutes']} min\n"); print("==", status)

try:
    # ------------------------------------------------------------ choose images
    files, pos_sites = {}, set()
    for chem in ("5mC", "5hmC", "6mA", "control"):
        f = os.path.join(a.features, f"{a.prefix}{chem}_{a.strand}.h5")
        if not os.path.exists(f): result["warnings"].append(f"missing {f}"); continue
        with h5py.File(f, "r") as h:
            lab = h["labels"][:]; names = h["ref_names"][:].astype(str); pos = h["ref_pos"][:]
        files[chem] = (f, lab, np.char.add(np.char.add(names, ":"), pos.astype(str)))
        if chem != "control": pos_sites |= set(files[chem][2][lab == 1].tolist())
    chosen = {}
    for chem, (f, lab, key) in files.items():
        idx = np.flatnonzero(lab == 1) if chem != "control" else np.flatnonzero(np.isin(key, list(pos_sites)))
        if len(idx) > a.max_per_class: idx = np.sort(rng.choice(idx, a.max_per_class, replace=False))
        chosen["none" if chem == "control" else chem] = (f, idx, key[idx])
        print(f"  {chem}: {len(idx):,} images at {len(set(key[idx].tolist())):,} sites", flush=True)
    classes = [c for c in ("none", "5mC", "5hmC", "6mA") if c in chosen and len(chosen[c][1])]
    if len(classes) < 2: raise RuntimeError(f"need at least two classes, have {classes}")
    result["classes"] = classes; result["n_images"] = {c: int(len(chosen[c][1])) for c in classes}

    def split_of(keys):
        h = np.array([int(hashlib.md5(k.encode()).hexdigest()[:6], 16) % 100 for k in keys]); return h >= 80     # True = test

    def probe(X, y, test, groups=None, tag=""):
        """Linear probe. Returns per-row metrics and, if groups given, metrics after averaging probabilities per group."""
        sc = StandardScaler().fit(X[~test]); clf = LogisticRegression(max_iter=300, C=1.0)
        clf.fit(sc.transform(X[~test]), y[~test]); P = clf.predict_proba(sc.transform(X[test])); yt = y[test]
        r = {"n_train": int((~test).sum()), "n_test": int(test.sum()), "acc": round(float((P.argmax(1) == yt).mean()), 4),
             "macro_f1": round(float(f1_score(yt, P.argmax(1), average="macro")), 4)}
        cm = np.zeros((len(classes), len(classes)), int)
        for t_, p_ in zip(yt, P.argmax(1)): cm[t_, p_] += 1
        r["confusion_rows_true_cols_pred"] = cm.tolist()
        if groups is not None:
            g = groups[test]; u, inv = np.unique(g, return_inverse=True); S = np.zeros((len(u), P.shape[1])); np.add.at(S, inv, P)
            yg = np.zeros(len(u), int); yg[inv] = yt
            r["pooled_by_image_acc"] = round(float((S.argmax(1) == yg).mean()), 4)
            r["pooled_by_image_macro_f1"] = round(float(f1_score(yg, S.argmax(1), average="macro")), 4)
        print(f"    {tag:42s} acc={r['acc']} macroF1={r['macro_f1']} pooled_acc={r.get('pooled_by_image_acc')}", flush=True)
        return r

    for spec in a.checkpoints.split(","):
        name, ck = spec.split("=", 1)
        if not os.path.exists(ck): result["warnings"].append(f"missing checkpoint {ck}"); continue
        print(f"== checkpoint {name}", flush=True)
        model, _ = load_model(ck, dev); model.eval(); cap = {}
        h1 = model.read_encoder.register_forward_hook(lambda m, i, o: cap.__setitem__("pre", o.detach()))
        h2 = model.encoder.register_forward_hook(lambda m, i, o: cap.__setitem__("post", o.detach()))
        site_rep, pre, post, raw, y_img, y_read, img_of_read, key_img = [], [], [], [], [], [], [], []
        n_img = 0
        for ci, c in enumerate(classes):
            f, idx, keys = chosen[c]
            with h5py.File(f, "r") as h: n = h["tensors"].shape[0]
            ds = PileupDataset([f], idx.astype(np.int64), [n], augment=False, seed=0, signal_noise_std=0.0, delta_channels=True, preload=True)      # one sequential read per file; no worker processes to fail
            dl = torch.utils.data.DataLoader(ds, batch_size=256, shuffle=False, num_workers=0)
            for xb, _ in dl:
                xb = xb.to(dev); B, C, Hh, Ww = xb.shape
                with torch.no_grad(): model(xb)
                pad = xb[:, 0].abs().sum(dim=2) < 1e-6; pad[:, 0] = False
                keep = (~pad).unsqueeze(-1).float(); enc = cap["post"]
                with torch.no_grad(): site_rep.append(model.norm((enc * keep).sum(1) / keep.sum(1).clamp(min=1.0)).cpu().numpy())
                rd = ~pad; rd[:, 0] = False                                    # real reads only, not the reference row
                bi, hi = torch.nonzero(rd, as_tuple=True)
                pre.append(cap["pre"].view(B, Hh, -1)[bi, hi].cpu().numpy()); post.append(enc[bi, hi].cpu().numpy())
                L = Ww // 21
                sigm = xb[:, 0].view(B, Hh, 21, L).mean(3)[bi, hi]; dw = xb[:, 1].view(B, Hh, 21, L).mean(3)[bi, hi]
                raw.append(torch.cat([sigm, dw], 1).cpu().numpy())
                y_read.append(np.full(len(bi), ci)); img_of_read.append(bi.cpu().numpy() + n_img); n_img += B
            y_img.append(np.full(len(idx), ci)); key_img.append(keys)
        h1.remove(); h2.remove()
        site_rep, pre, post, raw = (np.concatenate(v) for v in (site_rep, pre, post, raw))
        y_img, y_read, img_of_read, key_img = (np.concatenate(v) for v in (y_img, y_read, img_of_read, key_img))
        test_img = split_of(key_img); test_read = test_img[img_of_read]
        print(f"  {len(y_img):,} images, {len(y_read):,} reads; test sites hold {int(test_img.sum()):,} images", flush=True)
        R = {"n_reads": int(len(y_read)),
             "site_rep_96d": probe(site_rep, y_img, test_img, tag="pooled site representation"),
             "read_post_transformer": probe(post, y_read, test_read, img_of_read, "per-read, after cross-read Transformer"),
             "read_pre_transformer": probe(pre, y_read, test_read, img_of_read, "per-read, read encoder only"),
             "read_raw_signal_baseline": probe(raw, y_read, test_read, img_of_read, "per-read raw mean signal + dwell (42-d)")}
        if "5mC" in classes and "5hmC" in classes:                          # the hard pair on its own
            pr = np.isin(y_read, [classes.index("5mC"), classes.index("5hmC")]); pi = np.isin(y_img, [classes.index("5mC"), classes.index("5hmC")])
            full = classes; classes = ["5mC", "5hmC"]
            remap = lambda v: (v == full.index("5hmC")).astype(int)
            R["pair_5mC_vs_5hmC"] = {"site_rep_96d": probe(site_rep[pi], remap(y_img[pi]), test_img[pi], tag="5mC vs 5hmC, site representation"),
                                     "read_post_transformer": probe(post[pr], remap(y_read[pr]), test_read[pr], img_of_read[pr], "5mC vs 5hmC, per-read after Transformer"),
                                     "read_raw_signal_baseline": probe(raw[pr], remap(y_read[pr]), test_read[pr], img_of_read[pr], "5mC vs 5hmC, per-read raw baseline")}
            classes = full
        result["checkpoints"][name] = R
    finish("OK" if not result["warnings"] else "OK_WITH_WARNINGS")
except Exception as ex:
    import traceback; traceback.print_exc(); result["error"] = repr(ex); finish(f"FAILED: {ex!r}"); sys.exit(1)
