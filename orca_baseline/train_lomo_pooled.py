#!/usr/bin/env python3
"""
Pooled leave-one-modification-out across datasets (Bhargav's "LOMO pooled").

Pools the ONT synthetic conditions (per-site GT BEDs) with the UMBC SPO1
barcodes and holds out one modification type at a time. SPO1 labels come from
the biology rather than a BED: the phage genome carries 5hmU in place of
thymine genome-wide, so in the WGS barcodes every A/T reference position is
modified (T on the forward strand, A positions via reverse-strand reads) and
C/G positions are not. PCR barcodes are unmodified everywhere.

Modification pool:
  5mC, 5hmC, 6mA   from the synthetic conditions (BED-labeled)
  5hmU             from SPO1 WGS barcodes (center-base labeled)
Controls (always in training, never held out):
  synthetic control + SPO1 PCR barcodes

Uses the ORCA architecture and the diagnostics from train_orca_features.py
(pos_rate, lift, precision at matched recall, embedding rank). --grl-max
defaults to 0.1 per the GRL sweep: ORCA's default 1.0 collapses the embedding
and destabilizes across seeds.
"""
import os
import argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import TensorDataset, DataLoader, Subset

from train_orca_features import (ORCAcls, METRICS, evaluate, embedding_rank,
                                 grl_lambda, load_condition, IDX)

torch.backends.mkldnn.enabled = False

SYN_MODS = ["control", "5mC", "5hmC", "6mA"]
SPO1_WGS = ["barcode06", "barcode07"]
SPO1_PCR = ["barcode01", "barcode02", "barcode03", "barcode04", "barcode05"]


def load_spo1_barcode(spo1_dir, bc, modified, mod_bases="AT"):
    """One SPO1 barcode -> X (N,5,C), y (N,).

    modified=True (WGS): label 1 where the center reference base is in
    mod_bases (hmU replaces T genome-wide; A positions carry the signal via
    reverse-strand reads), else 0. modified=False (PCR): all 0.
    """
    csv = os.path.join(spo1_dir, f"{bc}_sub0.05", f"{bc}.merged.feature.per.site")
    df = pd.read_csv(csv)
    feat_cols = [c for c in df.columns if c not in IDX]
    assert len(feat_cols) % 5 == 0
    C = len(feat_cols) // 5
    X = np.nan_to_num(df[feat_cols].to_numpy(np.float32).reshape(-1, 5, C))
    if modified:
        center = df["kmer"].str[4].str.upper()   # 9-mer center = the site base
        y = center.isin(list(mod_bases)).to_numpy().astype(np.int64)
    else:
        y = np.zeros(len(df), dtype=np.int64)
    return X, y, C


def build_pool(args):
    """Return {group_name: (X, y)} plus the list of holdable modifications."""
    data = {}
    for m in SYN_MODS:
        X, y, C = load_condition(args.syn_feat_dir, m + args.syn_suffix, m,
                                 args.gt_dir)
        data["syn_" + m] = (X, y)
        print(f"  syn_{m}: {len(y)} sites, {int(y.sum())} pos")
    for bc in SPO1_WGS:
        X, y, C = load_spo1_barcode(args.spo1_dir, bc, True, args.mod_bases)
        data["spo1_" + bc] = (X, y)
        print(f"  spo1_{bc} (WGS/5hmU): {len(y)} sites, {int(y.sum())} pos")
    for bc in SPO1_PCR:
        X, y, C = load_spo1_barcode(args.spo1_dir, bc, False)
        data["spo1_" + bc] = (X, y)
        print(f"  spo1_{bc} (PCR/ctrl): {len(y)} sites, 0 pos")
    return data


# which groups carry each modification (held-out test set), and which groups
# are pure control (always kept in training)
MOD_GROUPS = {
    "5mC":  ["syn_5mC"],
    "5hmC": ["syn_5hmC"],
    "6mA":  ["syn_6mA"],
    "5hmU": ["spo1_" + bc for bc in SPO1_WGS],
}
CTRL_GROUPS = ["syn_control"] + ["spo1_" + bc for bc in SPO1_PCR]


def train_one(held, data, args, device, seed=0):
    import copy
    train_groups = CTRL_GROUPS + [g for m, gs in MOD_GROUPS.items()
                                  if m != held for g in gs]
    test_groups = MOD_GROUPS[held]
    dom_ids = {g: i for i, g in enumerate(train_groups)}

    Xtr = np.concatenate([data[g][0] for g in train_groups])
    ytr = np.concatenate([data[g][1] for g in train_groups])
    dtr = np.concatenate([np.full(len(data[g][1]), dom_ids[g])
                          for g in train_groups])
    Xte = np.concatenate([data[g][0] for g in test_groups])
    yte = np.concatenate([data[g][1] for g in test_groups])

    if args.cap and len(Xtr) > args.cap:
        idx = np.random.default_rng(seed).choice(len(Xtr), args.cap, False)
        Xtr, ytr, dtr = Xtr[idx], ytr[idx], dtr[idx]

    rng = np.random.default_rng(seed)
    perm = rng.permutation(len(Xtr))
    nval = max(args.batch, len(Xtr) // 10)
    vi, ti = perm[:nval], perm[nval:]

    def dl(X, y, d, idx=None, shuffle=False):
        ds = TensorDataset(torch.tensor(X, dtype=torch.float32),
                           torch.tensor(y, dtype=torch.long),
                           torch.tensor(d, dtype=torch.long))
        if idx is not None:
            ds = Subset(ds, idx)
        return DataLoader(ds, batch_size=args.batch, shuffle=shuffle,
                          drop_last=shuffle)

    dl_tr = dl(Xtr, ytr, dtr, ti, shuffle=True)
    dl_val = dl(Xtr, ytr, dtr, vi)
    dl_te = dl(Xte, yte, np.zeros(len(yte)))

    npos = int(ytr[ti].sum())
    pos_w = (len(ti) - npos) / max(1, npos)
    torch.manual_seed(seed)
    C = Xtr.shape[2]
    model = ORCAcls(in_ch=C, window=5, n_domains=len(train_groups)).to(device)
    print(f"[hold={held} seed={seed}] train {len(ti)} ({npos} pos, "
          f"{len(train_groups)} domains), test {len(yte)} ({int(yte.sum())} pos)")
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr)
    nll = nn.NLLLoss(weight=torch.tensor([1.0, float(pos_w)]).to(device))
    ce_dom = nn.NLLLoss()

    total = args.epochs * max(1, len(dl_tr)); step = 0
    best_val = -1.0; best_state = None
    for ep in range(1, args.epochs + 1):
        model.train()
        for x, y, d in dl_tr:
            x, y, d = x.to(device), y.to(device), d.to(device)
            lam = args.grl_max * grl_lambda(step, total)
            ce, dm = model(x, lam)
            loss = nll(ce, y) + args.dom_weight * ce_dom(dm, d)
            opt.zero_grad(); loss.backward(); opt.step(); step += 1
        vm = evaluate(model, dl_val, device)
        if vm["auprc"] == vm["auprc"] and vm["auprc"] > best_val:
            best_val = vm["auprc"]; best_state = copy.deepcopy(model.state_dict())
        print(f"  [{held}] ep {ep:2d} val_AUPRC={vm['auprc']:.3f} (best {best_val:.3f})")
    if best_state:
        model.load_state_dict(best_state)
    er = embedding_rank(model, dl_val, device)
    m, p, y = evaluate(model, dl_te, device, return_probs=True)
    m["rank99"], m["eff_rank"] = er["rank99"], er["eff_rank"]
    np.savez_compressed(os.path.join(args.out_dir, f"probs_{held}_seed{seed}.npz"),
                        probs=p, labels=y)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--syn-feat-dir", default="/fs/nexus-scratch/vgandhi/orca_feat")
    ap.add_argument("--syn-suffix", default="_rep1")
    ap.add_argument("--gt-dir",
                    default="/fs/nexus-scratch/bds062/data/ont-os/references")
    ap.add_argument("--spo1-dir",
                    default="/fs/nexus-scratch/vgandhi/orca_feat_spo1/run1_jan31/single_end")
    ap.add_argument("--mod-bases", default="AT",
                    help="reference center bases labeled 5hmU in SPO1 WGS")
    ap.add_argument("--out-dir", default="/fs/nexus-scratch/vgandhi/lomo_pooled")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--lr", type=float, default=5e-4)
    ap.add_argument("--cap", type=int, default=400000,
                    help="max pooled training sites (memory guard)")
    ap.add_argument("--grl-max", type=float, default=0.1)
    ap.add_argument("--dom-weight", type=float, default=1.0)
    ap.add_argument("--holdouts", default="5mC,5hmC,6mA,5hmU")
    ap.add_argument("--seeds", default="0,1,2")
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Device: {device}")

    data = build_pool(args)
    seeds = [int(s) for s in args.seeds.split(",")]
    agg = {}
    for held in args.holdouts.split(","):
        runs = [train_one(held, data, args, device, seed=s) for s in seeds]
        agg[held] = {k: (float(np.nanmean([r[k] for r in runs])),
                         float(np.nanstd([r[k] for r in runs])))
                     for k in METRICS}
        print(f"### {held}: AUPRC={agg[held]['auprc'][0]:.3f} "
              f"lift={agg[held]['lift'][0]:.1f}x "
              f"P@R30={agg[held]['prec_at_r30'][0]:.3f}\n")

    with open(os.path.join(args.out_dir, "lomo_metrics.tsv"), "w") as fh:
        fh.write("held_out\t" + "\t".join(
            f"{k}\t{k}_std" for k in METRICS) + "\n")
        for h, m in agg.items():
            fh.write(h + "\t" + "\t".join(
                f"{m[k][0]:.4f}\t{m[k][1]:.4f}" for k in METRICS) + "\n")
    print("Saved lomo_metrics.tsv")


if __name__ == "__main__":
    main()
