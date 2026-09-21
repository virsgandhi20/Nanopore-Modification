#!/usr/bin/env python3
"""Per-read modification typing on signal windows (extract_read_windows.py output).

One script, many experiments: which samples train, which classes exist, which
input channels the model may see, and which held-out sets are scored are all
arguments. Splits are BY SITE (chrom:pos:strand, shared across samples so a
modified site and its matched unmodified control always land in the same split).

Reports per test set: per-read accuracy / macro-F1 / confusion / one-vs-rest
AUROC, modified-vs-none AUROC, site-level metrics (mean probability over the
reads of a site), open-set detection of labels the model was never trained on,
and a stoichiometry simulation (mix modified and control reads of the same
site at known fractions, estimate the fraction from per-read calls).

Built to finish: every evaluation is isolated, and metrics.json + status.txt
are always written, with the reason if something went wrong.
"""
import argparse, hashlib, json, os, sys, time, traceback
from collections import Counter, defaultdict

import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import f1_score, roc_auc_score

ap = argparse.ArgumentParser()
ap.add_argument("--name", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--train", required=True, help="comma list of npz")
ap.add_argument("--eval", default="", help="name=npz+npz[,name=npz...] scored in full, never trained on")
ap.add_argument("--classes", required=True, help="comma list, 'none' first")
ap.add_argument("--inputs", default="sig,dwell,seq")
ap.add_argument("--groups", default="", help="keep only these site groups (train and test)")
ap.add_argument("--train-samples", default="", help="substring filter: only these samples train; the rest of --train is test-only")
ap.add_argument("--epochs", type=int, default=8); ap.add_argument("--batch", type=int, default=512)
ap.add_argument("--lr", type=float, default=2e-3); ap.add_argument("--max-per-class", type=int, default=150000)
ap.add_argument("--merge", default="", help="new=a+b+c[,new2=...]: relabel before anything else (e.g. mod=5mC+5hmC+6mA+4mC)")
ap.add_argument("--renorm", choices=["none", "window"], default="none",
                help="window: re-centre and re-scale each 21-base window by its own median / MAD, removing "
                     "read-level normalization differences (139-base oligo reads vs multi-kb genomic reads)")
ap.add_argument("--crop", type=int, default=0, help="keep only the centre +/- CROP bases of each window (0 = all)")
ap.add_argument("--seed", type=int, default=0); ap.add_argument("--note", default="")
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
t0 = time.time()
LOG = open(os.path.join(a.out, "log.txt"), "w")
def say(*x):
    s = " ".join(str(v) for v in x); print(s, flush=True); LOG.write(s + "\n"); LOG.flush()
result = {"name": a.name, "note": a.note, "args": vars(a), "sets": {}, "warnings": []}
def finish(status):
    result["status"] = status; result["minutes"] = round((time.time() - t0) / 60, 1)
    json.dump(result, open(os.path.join(a.out, "metrics.json"), "w"), indent=1, default=str)
    open(os.path.join(a.out, "status.txt"), "w").write(f"{a.name}\t{status}\t{result['minutes']} min\n")
    say(f"== {a.name}: {status} ({result['minutes']} min)")

CLASSES = a.classes.split(","); CI = {c: i for i, c in enumerate(CLASSES)}
INPUTS = set(a.inputs.split(",")); GROUPS = set(g for g in a.groups.split(",") if g)
MERGE = {old: new for spec in a.merge.split(",") if spec for new, olds in [spec.split("=", 1)] for old in olds.split("+")}
torch.manual_seed(a.seed); np.random.seed(a.seed); rng = np.random.default_rng(a.seed)
dev = "cuda" if torch.cuda.is_available() else "cpu"


def load(paths):
    parts = []
    for p in paths:
        if not os.path.exists(p):
            result["warnings"].append(f"missing input {p}"); say("WARNING missing", p); continue
        z = np.load(p, allow_pickle=False)
        sample = str(z["sample"]); keep = np.ones(len(z["label"]), bool)
        if GROUPS: keep &= np.isin(z["group"], list(GROUPS))
        sid = z["site"][keep]
        key = np.char.add(np.char.add(np.char.add(z["site_chrom"][sid], ":"), z["site_pos"][sid].astype(str)),
                          np.char.add(":", z["site_strand"][sid]))
        lab_all = z["label"]
        if MERGE: lab_all = np.array([MERGE.get(l, l) for l in lab_all.tolist()])
        sig, dwl, bas = z["sig"][keep], z["dwell"][keep], z["base"][keep]
        if a.crop:
            c = sig.shape[1] // 2; sig, dwl, bas = sig[:, c - a.crop:c + a.crop + 1], dwl[:, c - a.crop:c + a.crop + 1], bas[:, c - a.crop:c + a.crop + 1]
        if a.renorm == "window":
            x = sig.astype(np.float32); flat = x.reshape(len(x), -1)
            med = np.median(flat, 1)[:, None, None]; mad = (np.median(np.abs(flat - med[:, :, 0]), 1) * 1.4826 + 1e-3)[:, None, None]
            sig = np.clip((x - med) / mad, -6, 6).astype(np.float16)
        parts.append(dict(sig=sig, dwell=dwl, base=bas, label=lab_all[keep],
                          group=z["group"][keep], key=key, sample=np.full(keep.sum(), sample)))
        say(f"  loaded {sample}: {int(keep.sum()):,} windows {dict(Counter(lab_all[keep].tolist()))}")
    if not parts: return None
    return {k: np.concatenate([p[k] for p in parts]) for k in parts[0]}


def split_of(keys):
    u, inv = np.unique(keys, return_inverse=True)
    h = np.array([int(hashlib.md5(k.encode()).hexdigest()[:6], 16) % 100 for k in u])
    return np.where(h < 70, 0, np.where(h < 85, 1, 2))[inv]          # 0 train, 1 val, 2 test


class Net(nn.Module):
    def __init__(self, cin, ncls):
        super().__init__()
        def blk(i, o, k, s=1): return nn.Sequential(nn.Conv1d(i, o, k, s, k // 2), nn.BatchNorm1d(o), nn.GELU())
        self.f = nn.Sequential(nn.BatchNorm1d(cin), blk(cin, 64, 7), blk(64, 64, 5), blk(64, 128, 5, 2),
                               blk(128, 128, 3), blk(128, 192, 3, 2), blk(192, 192, 3))
        self.h = nn.Sequential(nn.Linear(384, 128), nn.GELU(), nn.Dropout(0.2), nn.Linear(128, ncls))
    def forward(self, x):
        z = self.f(x); return self.h(torch.cat([z.mean(-1), z.amax(-1)], 1))


def batch_x(d, idx):
    ch = []
    L = d["sig"].shape[2]
    if "sig" in INPUTS: ch.append(torch.from_numpy(d["sig"][idx].astype(np.float32)).flatten(1).unsqueeze(1))
    if "dwell" in INPUTS: ch.append(torch.from_numpy(d["dwell"][idx].astype(np.float32)).repeat_interleave(L, 1).unsqueeze(1))
    if "seq" in INPUTS:
        b = torch.from_numpy(d["base"][idx].astype(np.int64))
        oh = torch.nn.functional.one_hot(b, 5)[..., :4].float().permute(0, 2, 1).repeat_interleave(L, 2)
        ch.append(oh)
    return torch.cat(ch, 1)


@torch.no_grad()
def predict(model, d, idx):
    model.eval(); out = []
    for i in range(0, len(idx), 4096):
        out.append(torch.softmax(model(batch_x(d, idx[i:i + 4096]).to(dev)), 1).cpu().numpy())
    return np.concatenate(out) if out else np.zeros((0, len(CLASSES)))


def safe_auc(y, s):
    try: return round(float(roc_auc_score(y, s)), 4) if 0 < y.sum() < len(y) else None
    except Exception: return None


def score(d, idx, P, tag):
    """All metrics for one test set. idx indexes d; P are class probabilities for idx."""
    lab = d["label"][idx]; known = np.isin(lab, CLASSES)
    r = {"n_windows": int(len(idx)), "labels": dict(Counter(lab.tolist()))}
    if known.sum():
        y = np.array([CI[l] for l in lab[known]]); pk = P[known]; yhat = pk.argmax(1)
        present = sorted(set(y.tolist()))
        r["read_acc"] = round(float((y == yhat).mean()), 4)
        r["read_macro_f1"] = round(float(f1_score(y, yhat, labels=present, average="macro")), 4)
        r["read_recall_by_class"] = {CLASSES[c]: round(float((yhat[y == c] == c).mean()), 4) for c in present}
        cm = np.zeros((len(CLASSES), len(CLASSES)), int)
        for t, p_ in zip(y, yhat): cm[t, p_] += 1
        r["confusion_rows_true_cols_pred"] = {"classes": CLASSES, "matrix": cm.tolist()}
        r["read_auroc_ovr"] = {CLASSES[c]: safe_auc((y == c).astype(int), pk[:, c]) for c in present}
        if "none" in CI: r["read_auroc_mod_vs_none"] = safe_auc((y != CI["none"]).astype(int), 1 - pk[:, CI["none"]])
        # site level: mean probability over the reads of one site in one sample
        sk = np.char.add(np.char.add(d["sample"][idx][known], "|"), d["key"][idx][known])
        u, inv = np.unique(sk, return_inverse=True); cnt = np.bincount(inv)
        good = cnt >= 5
        if good.sum() >= 20:
            ps = np.zeros((len(u), pk.shape[1])); np.add.at(ps, inv, pk); ps /= cnt[:, None]
            ys = np.zeros(len(u), int); ys[inv] = y
            ps, ys = ps[good], ys[good]
            r["site_n"] = int(good.sum()); r["site_acc"] = round(float((ps.argmax(1) == ys).mean()), 4)
            r["site_macro_f1"] = round(float(f1_score(ys, ps.argmax(1), labels=sorted(set(ys.tolist())), average="macro")), 4)
            if "none" in CI: r["site_auroc_mod_vs_none"] = safe_auc((ys != CI["none"]).astype(int), 1 - ps[:, CI["none"]])
        # stoichiometry: same site seen modified (one sample) and unmodified (another)
        if "none" in CI:
            keyk = d["key"][idx][known]; pm = 1 - pk[:, CI["none"]]; ismod = y != CI["none"]
            bysite = defaultdict(lambda: ([], []))
            for kk, m, v in zip(keyk, ismod, pm): bysite[kk][0 if m else 1].append(v)
            pairs = [(np.array(m_), np.array(n_)) for m_, n_ in bysite.values() if len(m_) >= 6 and len(n_) >= 6]
            if len(pairs) >= 8:
                fr = [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]; est = {f: [] for f in fr}
                for m_, n_ in pairs[:3000]:
                    for f in fr:
                        km = int(round(16 * f))
                        x = np.concatenate([rng.choice(m_, km), rng.choice(n_, 16 - km)]) if 0 < km < 16 else (rng.choice(m_, 16) if km else rng.choice(n_, 16))
                        est[f].append(float((x > 0.5).mean()))
                tf = np.repeat(fr, [len(est[f]) for f in fr]); ef = np.concatenate([est[f] for f in fr])
                r["stoichiometry"] = {"sites": len(pairs), "reads_per_mix": 16,
                                      "mean_estimate_by_true_fraction": {str(f): round(float(np.mean(est[f])), 3) for f in fr},
                                      "mae": round(float(np.abs(tf - ef).mean()), 4), "pearson_r": round(float(np.corrcoef(tf, ef)[0, 1]), 4)}
    unk = ~known
    if unk.sum() >= 50 and known.sum() >= 50:              # open set: labels never trained on
        conf = P.max(1)
        r["openset"] = {"unknown_labels": dict(Counter(lab[unk].tolist())),
                        "auroc_unknown_by_low_confidence": safe_auc(unk.astype(int), -conf),
                        "unknown_called_as": {CLASSES[c]: round(float((P[unk].argmax(1) == c).mean()), 4) for c in range(len(CLASSES))},
                        "mean_confidence_known": round(float(conf[known].mean()), 4), "mean_confidence_unknown": round(float(conf[unk].mean()), 4)}
        knownmod = known & (lab != "none")
        if knownmod.sum() >= 50:
            m = unk | knownmod
            r["openset"]["auroc_unknown_vs_known_modified"] = safe_auc(unk[m].astype(int), -conf[m])
    say(f"  [{tag}] n={r['n_windows']:,} read_acc={r.get('read_acc')} macroF1={r.get('read_macro_f1')} "
        f"modAUROC={r.get('read_auroc_mod_vs_none')} site_acc={r.get('site_acc')} recall={r.get('read_recall_by_class')}")
    return r


try:
    say(f"== {a.name} | classes {CLASSES} | inputs {sorted(INPUTS)} | device {dev} | {a.note}")
    D = load(a.train.split(","))
    if D is None: raise RuntimeError("none of the training inputs exist")
    sp = split_of(D["key"])
    trainable = np.isin(D["label"], CLASSES)
    if a.train_samples:
        subs = a.train_samples.split(",")
        insamp = np.array([any(s in x for s in subs) for x in D["sample"]])
        sp = np.where(insamp, sp, 2)                       # other samples are test-only, in full
        trainable &= insamp
    tr = np.flatnonzero((sp == 0) & trainable); va = np.flatnonzero((sp == 1) & trainable); te = np.flatnonzero(sp == 2)
    # cap and report class balance
    keep = []
    for c in CLASSES:
        ii = tr[D["label"][tr] == c]
        if len(ii) > a.max_per_class: ii = rng.choice(ii, a.max_per_class, replace=False)
        keep.append(ii)
    counts = {c: int(len(k)) for c, k in zip(CLASSES, keep)}
    say(f"  train windows by class: {counts}; val {len(va):,}; test {len(te):,}")
    if sum(v > 0 for v in counts.values()) < 2: raise RuntimeError(f"fewer than two classes have training data: {counts}")
    for c, v in counts.items():
        if v == 0: result["warnings"].append(f"class {c} has no training windows")
    tr = np.concatenate(keep); result["train_counts"] = counts
    w = torch.tensor([0.0 if counts[c] == 0 else len(tr) / (len(CLASSES) * counts[c]) for c in CLASSES], dtype=torch.float32).to(dev)
    ytr = torch.tensor([CI[l] for l in D["label"][tr]])
    cin = batch_x(D, tr[:2]).shape[1]
    model = Net(cin, len(CLASSES)).to(dev)
    opt = torch.optim.AdamW(model.parameters(), lr=a.lr, weight_decay=1e-2)
    steps = a.epochs * max(1, len(tr) // a.batch); sched = torch.optim.lr_scheduler.OneCycleLR(opt, a.lr, total_steps=steps + 5)
    lossf = nn.CrossEntropyLoss(weight=w); best, best_state = -1, None
    for ep in range(a.epochs):
        model.train(); perm = rng.permutation(len(tr)); tot = 0.0; nb = 0
        for i in range(0, len(perm) - a.batch + 1, a.batch) if len(perm) >= a.batch else [0]:
            j = perm[i:i + a.batch]
            loss = lossf(model(batch_x(D, tr[j]).to(dev)), ytr[j].to(dev))
            opt.zero_grad(set_to_none=True); loss.backward(); opt.step(); sched.step(); tot += float(loss); nb += 1
        if len(va):
            pv = predict(model, D, va); yv = np.array([CI[l] for l in D["label"][va]])
            f1 = float(f1_score(yv, pv.argmax(1), labels=sorted(set(yv.tolist())), average="macro"))
        else: f1 = -tot
        say(f"  epoch {ep+1}/{a.epochs} loss {tot/max(nb,1):.4f} val_macroF1 {f1:.4f} ({(time.time()-t0)/60:.1f} min)")
        if f1 > best: best, best_state = f1, {k: v.detach().cpu().clone() for k, v in model.state_dict().items()}
    model.load_state_dict(best_state); result["best_val_macro_f1"] = round(best, 4)
    torch.save({"state": best_state, "classes": CLASSES, "inputs": sorted(INPUTS), "cin": cin}, os.path.join(a.out, "model.pt"))

    def run_set(tag, d, idx):
        try:
            if len(idx) == 0: result["sets"][tag] = {"skipped": "no windows"}; return
            result["sets"][tag] = score(d, idx, predict(model, d, idx), tag)
        except Exception as ex:
            result["sets"][tag] = {"error": repr(ex)}; result["warnings"].append(f"{tag}: {ex!r}"); say(traceback.format_exc())
    run_set("test", D, te)
    for smp in sorted(set(D["sample"][te].tolist())):
        run_set("test/" + smp, D, te[D["sample"][te] == smp])
    for spec in [s for s in a.eval.split(",") if s]:
        nm, paths = spec.split("=", 1)
        E = load(paths.split("+"))
        if E is None: result["sets"][nm] = {"skipped": "inputs missing"}; continue
        run_set(nm, E, np.arange(len(E["label"])))
    finish("OK" if not result["warnings"] else "OK_WITH_WARNINGS")
except Exception as ex:
    say(traceback.format_exc()); result["error"] = repr(ex); finish(f"FAILED: {ex!r}")
    sys.exit(1)
