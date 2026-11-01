#!/usr/bin/env python3
"""Per-read modification TYPING on RawMod's own encoder, trained on RawMod's own features.

The Sep 21 probe (rawmod_embedding_probe.py) showed that a frozen RawMod checkpoint already carries
modification type in its per-read embeddings (linear probe macro-F1 0.71 to 0.77). This script is the
next step the journal asked for: the like-for-like comparison with the per-read model of
train_read_typing.py, on RawMod's inputs.

  * A typing head (LayerNorm, Linear, GELU, Linear) sits on the per-read rows of RawMod's cross-read
    Transformer output `enc` (B, H, d), i.e. the tap point the probe called "after the Transformer".
  * --mode frozen   : encoder fixed, only the head trains (a non-linear probe).
  * --mode finetune : encoder and head train together; the original detection logit keeps a BCE term
    (--detect-weight) so detection is not traded away, and detection AUROC is reported before and after.
  * Labels are per READ: every read of a modified oligo image carries the sample's chemistry, every read
    of a control image at the same positions is 'none'. Splits are by site (contig:pos), 80/20 by hash,
    exactly as the Sep 21 scripts, so a site never appears on both sides.
  * Open set: chemistries present in --files but absent from --classes are never trained on; their reads
    are scored with 1 - max softmax as an "unknown" score (AUROC known-vs-unknown, and the fraction flagged
    at the threshold that keeps 95% of known test reads). Start from a leave-out checkpoint
    (checkpoints/leaveout/ont5hmC.pt, region_holdout/loco_5hmC.pt) to make the encoder itself naive.
  * --other: chemistries whose reads train an explicit 'other' class (appended to --classes), so the head can
    say "none of the known types". Chemistries in neither list stay unknown and are only tested, which asks
    whether an 'other' trained on one foreign chemistry (say 6mA) catches a different one (say 5hmC).
  * Three unknown scores are reported for every run: 1 - max softmax over the known classes, the head's own
    P(other) when --other is given, and the Mahalanobis distance of a read's hidden-layer feature to the nearest
    training-class centroid (shared covariance), which needs no 'other' class at all.
  * --eval-files: a second set of files (e.g. the replicate-2 libraries) scored in full, never trained on.

Inputs are RawMod feature files (featurization.py HDF5: tensors (N,H,W,9), labels, ref_names, ref_pos),
either Bhargav's strand-resolved set (ONT_<chem>_<strand>.h5) or the benchmark pipeline's
work/<sample>/rawmod_features/<tag>/features_plus.h5. The released package supplies the dataset class
(rawmod.model.PileupDataset adds the two delta channels) and the checkpoint loader
(scripts/evaluate/score_genome.load_model). Never raises on a missing file: it is recorded and skipped.

  python rawmod_typing_head.py --selftest --repo <RawMod clone> --checkpoint <ckpt> --out /tmp/t   # synthetic data, CPU
"""
import argparse, hashlib, json, os, sys, time, tempfile
from pathlib import Path

import h5py
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import f1_score, roc_auc_score

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--repo', default=os.environ.get('RAWMOD_REPO', '/fs/cbcb-lab/storm/vgandhi/rawmod-test/RawMod'),
                help='clone of the released RawMod repository (rawmod package + scripts/evaluate/score_genome.py)')
ap.add_argument('--checkpoint', default='', help='RawMod checkpoint (.pt); default: <repo>/checkpoints/region_holdout/mixed.pt')
ap.add_argument('--files', default='', help='chem=path[,chem=path...]; "none" = the unmodified control; chems not in --classes are unknown (open-set eval only)')
ap.add_argument('--eval-files', default='', help='same format; scored in full, never trained on (e.g. replicate-2 libraries)')
ap.add_argument('--classes', default='none,5mC,5hmC,6mA', help="typing classes, 'none' first")
ap.add_argument('--other', default='', help="chemistries (in --files) trained as an explicit 'other' class; chems in neither --classes nor --other stay unknown and are only tested")
ap.add_argument('--mode', choices=['frozen', 'finetune'], default='frozen')
ap.add_argument('--epochs', type=int, default=6)
ap.add_argument('--batch', type=int, default=128)
ap.add_argument('--lr', type=float, default=1e-4, help='encoder learning rate (finetune)')
ap.add_argument('--head-lr', type=float, default=1e-3)
ap.add_argument('--detect-weight', type=float, default=1.0, help='weight of the detection BCE term in finetune mode')
ap.add_argument('--max-per-class', type=int, default=3000, help='images per class (training classes and unknowns alike)')
ap.add_argument('--test-frac', type=int, default=20, help='percent of sites held out (by hash of contig:pos)')
ap.add_argument('--seed', type=int, default=0)
ap.add_argument('--device', default='')
ap.add_argument('--out', required=True)
ap.add_argument('--name', default='')
ap.add_argument('--note', default='')
ap.add_argument('--selftest', action='store_true', help='build synthetic feature files and run end to end (CPU, ~1 min)')
a = ap.parse_args()

t0 = time.time()
os.makedirs(a.out, exist_ok=True)
dev = a.device or ('cuda' if torch.cuda.is_available() else 'cpu')
torch.manual_seed(a.seed); rng = np.random.default_rng(a.seed)
result = {'name': a.name, 'args': vars(a), 'device': dev, 'warnings': []}


def finish(status):
    result['status'] = status; result['minutes'] = round((time.time() - t0) / 60, 1)
    json.dump(result, open(os.path.join(a.out, 'metrics.json'), 'w'), indent=1, default=str)
    open(os.path.join(a.out, 'status.txt'), 'w').write(f'{a.name or "rawmod_typing_head"}\t{status}\t{result["minutes"]} min\n')
    print('==', status, flush=True)


def warn(msg):
    result['warnings'].append(msg); print('  [warn]', msg, flush=True)


# ------------------------------------------------------------------------------------------- synthetic data
def make_synthetic(d, chems=('none', '5mC', '5hmC', '6mA', '4mC'), n=96, H=16, W=21, L=10, C=9, reads=10):
    """Feature files in featurization.py's layout where the chemistries differ by a signal offset at the centre."""
    files = {}
    for ci, chem in enumerate(chems):
        X = np.zeros((n, H, W * L, C), np.float32)
        ref = rng.normal(0, 1, (n, 1, W * L)).astype(np.float32)        # expected level, row 0
        X[:, 0, :, 0] = ref[:, 0]
        X[:, 1:reads + 1, :, 0] = ref + rng.normal(0, 0.6, (n, reads, W * L))
        if chem != 'none':
            c = W // 2
            X[:, 1:reads + 1, c * L:(c + 1) * L, 0] += 1.2 * np.cos(ci) + 0.8 * np.sin(2 * ci)   # class-specific shift
        X[:, 1:reads + 1, :, 1] = rng.normal(1.5, 0.3, (n, reads, W * L))                  # dwell
        base = rng.integers(0, 4, (n, 1, W))
        oh = np.eye(4, dtype=np.float32)[base][:, :, :, :].repeat(L, axis=2)              # (n,1,W*L,4)
        X[:, :reads + 1, :, 2:6] = oh
        X[:, 1:reads + 1, :, 6] = 1.0; X[:, 1:reads + 1, :, 7] = 0.9; X[:, 1:reads + 1, :, 8] = 1.0
        p = str(Path(d) / f'SYN_{chem}_plus.h5')
        with h5py.File(p, 'w') as hf:
            hf.create_dataset('tensors', data=X, compression='gzip', chunks=(min(64, n), H, W * L, C))
            hf.create_dataset('labels', data=np.full(n, 0 if chem == 'none' else 1, np.int8))
            hf.create_dataset('ref_names', data=np.array(['oligo1'] * n, dtype='S16'))
            hf.create_dataset('ref_pos', data=np.arange(n, dtype=np.int64) * 7 + 100)
            hf.create_dataset('n_reads', data=np.full(n, reads, np.int32))
            hf.create_dataset('image_idx', data=np.arange(n, dtype=np.int64))
            for k, v in dict(half_window=W // 2, L=L, max_reads=H - 1, W=W, n_channels=C, height=H, n_positions=n,
                             n_images=n, max_images_per_base=1, center_idx=6, normalization='MAD', ref_row='row_0',
                             label_semantics='binary_modified_vs_unmodified', multi_image=True).items():
                hf.attrs[k] = v
            hf.attrs['channel_names'] = ['raw_signal', 'dwell_log1p', 'is_A', 'is_C', 'is_G', 'is_T', 'strand', 'mapq_norm', 'matches_ref']
        files[chem] = p
    return files


# ------------------------------------------------------------------------------------------- data
def parse_files(spec):
    out = {}
    for item in filter(None, spec.split(',')):
        chem, path = item.split('=', 1)
        if not os.path.exists(path):
            warn(f'missing {chem}: {path}'); continue
        out[chem] = path
    return out


def choose_images(files, max_per_class):
    """chem -> (path, image indices, site keys). Chemistry files: their positives. Control: images at any
    position that is positive in some chemistry file (falls back to all control images if none overlap)."""
    info = {}
    for chem, path in files.items():
        with h5py.File(path, 'r') as h:
            lab = h['labels'][:].astype(int); names = h['ref_names'][:].astype(str); pos = h['ref_pos'][:]
        info[chem] = (path, lab, np.char.add(np.char.add(names, ':'), pos.astype(str)))
    pos_sites = set()
    for chem, (p, lab, keys) in info.items():
        if chem != 'none':
            pos_sites |= set(keys[lab == 1].tolist())
    chosen = {}
    for chem, (p, lab, keys) in info.items():
        if chem == 'none':
            idx = np.flatnonzero(np.isin(keys, list(pos_sites)))
            if len(idx) == 0:
                warn(f'{p}: no control image shares a position with a positive; using all {len(keys)} control images'); idx = np.arange(len(keys))
        else:
            idx = np.flatnonzero(lab == 1)
        if len(idx) == 0:
            warn(f'{chem}: no images in {p}'); continue
        if len(idx) > max_per_class:
            idx = np.sort(rng.choice(idx, max_per_class, replace=False))
        chosen[chem] = (p, idx, keys[idx])
        print(f'  {chem:6s} {len(idx):6,} images at {len(set(keys[idx].tolist())):6,} sites  {os.path.basename(p)}', flush=True)
    return chosen


def load_tensors(path, idx):
    """(N, 11, H, W) float32 through RawMod's own dataset class (adds the two delta channels)."""
    from rawmod.model import PileupDataset
    with h5py.File(path, 'r') as h:
        n = h['tensors'].shape[0]
    ds = PileupDataset([path], idx.astype(np.int64), [n], augment=False, seed=0, signal_noise_std=0.0,
                       delta_channels=True, preload=True)
    dl = torch.utils.data.DataLoader(ds, batch_size=256, shuffle=False, num_workers=0)
    return torch.cat([xb.float() for xb, _ in dl])


def site_split(keys, frac):
    h = np.array([int(hashlib.md5(k.encode()).hexdigest()[:6], 16) % 100 for k in keys])
    return h >= (100 - frac)            # True = test


# ------------------------------------------------------------------------------------------- model
def read_mask(xb):
    pad = xb[:, 0].abs().sum(dim=2) < 1e-6; pad[:, 0] = False     # as ConvFormerV2.forward
    rd = ~pad; rd[:, 0] = False                                   # real reads only, never the reference row
    return rd


try:
    sys.path.insert(0, a.repo); sys.path.insert(0, os.path.join(a.repo, 'scripts', 'evaluate'))
    from score_genome import load_model   # noqa: E402  (RawMod's own checkpoint loader)
    ckpt = a.checkpoint or os.path.join(a.repo, 'checkpoints', 'region_holdout', 'mixed.pt')
    if not os.path.exists(ckpt):
        raise FileNotFoundError(ckpt)
    result['checkpoint'] = ckpt
    if a.selftest:
        tmp = tempfile.mkdtemp(prefix='typing_selftest_')
        a.files = ','.join(f'{c}={p}' for c, p in make_synthetic(tmp).items())
        a.eval_files = ''; a.max_per_class = min(a.max_per_class, 96); a.epochs = min(a.epochs, 2); a.batch = 32
        a.other = a.other or '4mC'                            # trained 'other'; 5hmC stays unknown when --classes omits it
        print(f'selftest: synthetic files in {tmp}', flush=True)

    classes = a.classes.split(',')
    files = parse_files(a.files)
    if not files:
        raise RuntimeError('no feature files found')
    print('images:', flush=True)
    chosen = choose_images(files, a.max_per_class)
    known = [c for c in classes if c in chosen]
    other_chems = [c for c in filter(None, a.other.split(',')) if c in chosen and c not in classes]
    if other_chems:
        known = known + ['other']                          # one extra output, fed by every --other chemistry
    unknown = [c for c in chosen if c not in classes and c not in other_chems]
    if len(known) < 2:
        raise RuntimeError(f'need at least two training classes, have {known}')
    result['classes'] = known; result['unknown_chems'] = unknown; result['other_chems'] = other_chems
    result['n_images'] = {c: int(len(chosen[c][1])) for c in chosen}
    OTHER = known.index('other') if 'other' in known else -1

    def label_of(c):
        return known.index(c) if c in classes else (OTHER if c in other_chems else -1)

    X, y_img, keys, is_unknown = [], [], [], []
    for c in [c for c in classes if c in chosen] + other_chems + unknown:
        p, idx, k = chosen[c]
        X.append(load_tensors(p, idx)); y_img.append(np.full(len(idx), label_of(c))); keys.append(k)
        is_unknown.append(np.full(len(idx), c in unknown))
    X = torch.cat(X); y_img = np.concatenate(y_img); keys = np.concatenate(keys); is_unknown = np.concatenate(is_unknown)
    test = site_split(keys, a.test_frac) | is_unknown          # unknown chemistries are never trained on
    train_idx = np.flatnonzero(~test); test_idx = np.flatnonzero(test & ~is_unknown); unk_idx = np.flatnonzero(is_unknown)
    print(f'  tensors {tuple(X.shape)}  train images {len(train_idx):,}  test images {len(test_idx):,}  unknown images {len(unk_idx):,}', flush=True)

    ev = {}
    for c, p in parse_files(a.eval_files).items():
        with h5py.File(p, 'r') as h:
            lab = h['labels'][:].astype(int)
        idx = np.flatnonzero(lab == 1) if c != 'none' else np.arange(len(lab))
        if len(idx) > a.max_per_class:
            idx = np.sort(rng.choice(idx, a.max_per_class, replace=False))
        ev[c] = (load_tensors(p, idx), label_of(c))
        print(f'  eval {c:6s} {len(idx):6,} images  {os.path.basename(p)}', flush=True)

    model, arch = load_model(ckpt, dev)
    model.eval()                                           # dropout is 0 at load; eval also keeps forward() to one output
    for prm in model.parameters():
        prm.requires_grad_(a.mode == 'finetune')
    d = int(model.norm.weight.shape[0])
    class Head(nn.Module):
        def __init__(self, d, k):
            super().__init__(); self.norm = nn.LayerNorm(d); self.l1 = nn.Linear(d, 128); self.drop = nn.Dropout(0.1); self.l2 = nn.Linear(128, k)
        def feat(self, x):                                 # the 128-d hidden feature the Mahalanobis score lives in
            return F.gelu(self.l1(self.norm(x)))
        def forward(self, x):
            return self.l2(self.drop(self.feat(x)))
    head = Head(d, len(known)).to(dev)
    cap = {}
    hook = model.encoder.register_forward_hook(lambda m, i, o: cap.__setitem__('enc', o))

    def forward(xb, grad):
        with torch.set_grad_enabled(grad):
            out = model(xb)
            logit = (out[0] if isinstance(out, (tuple, list)) else out).reshape(len(xb))
        return logit, cap['enc']

    # per-read class weights from the training images (reads per image are near constant, so image counts do)
    cnt = np.bincount(y_img[train_idx], minlength=len(known)).astype(float)
    cw = torch.tensor((cnt.sum() / np.maximum(cnt, 1)) / len(known), dtype=torch.float32, device=dev)
    params = [{'params': head.parameters(), 'lr': a.head_lr}]
    if a.mode == 'finetune':
        params.append({'params': [p for p in model.parameters() if p.requires_grad], 'lr': a.lr})
    opt = torch.optim.AdamW(params, weight_decay=1e-4)

    def evaluate(idx_or_tensor, y_true=None):
        """Per-read and per-image probabilities plus detection logits for a set of images."""
        head.eval(); P_read, img_of_read, P_img, det, feats = [], [], [], [], []
        src = X[idx_or_tensor] if isinstance(idx_or_tensor, np.ndarray) else idx_or_tensor
        with torch.no_grad():
            for s in range(0, len(src), 256):
                xb = src[s:s + 256].to(dev); logit, enc = forward(xb, False); rd = read_mask(xb)
                fe = head.feat(enc[rd]); pr = F.softmax(head.l2(fe), dim=1)
                bi = torch.nonzero(rd, as_tuple=True)[0]
                S = torch.zeros(len(xb), len(known), device=dev).index_add_(0, bi, pr)
                S = S / torch.bincount(bi, minlength=len(xb)).clamp(min=1).unsqueeze(1).float()
                P_read.append(pr.cpu().numpy()); img_of_read.append(bi.cpu().numpy() + s); P_img.append(S.cpu().numpy()); det.append(logit.cpu().numpy())
                feats.append(fe.cpu().numpy())
        head.train()
        return np.concatenate(P_read), np.concatenate(img_of_read), np.concatenate(P_img), np.concatenate(det), np.concatenate(feats)

    def metrics(P_read, img_of_read, P_img, y):
        yr = y[img_of_read]; pr = P_read.argmax(1); pi = P_img.argmax(1)
        cm = np.zeros((len(known), len(known)), int)
        for t_, p_ in zip(yr, pr):
            cm[t_, p_] += 1
        r = {'n_images': int(len(y)), 'n_reads': int(len(yr)),
             'read_acc': round(float((pr == yr).mean()), 4), 'read_macro_f1': round(float(f1_score(yr, pr, average='macro')), 4),
             'read_recall': {known[c]: round(float((pr[yr == c] == c).mean()), 4) for c in range(len(known)) if (yr == c).any()},
             'image_acc': round(float((pi == y).mean()), 4), 'image_macro_f1': round(float(f1_score(y, pi, average='macro')), 4),
             'confusion_reads_rows_true_cols_pred': cm.tolist()}
        if '5mC' in known and '5hmC' in known:
            i5, ih = known.index('5mC'), known.index('5hmC'); sel = np.isin(yr, [i5, ih])
            if sel.any():
                pp = P_read[sel][:, ih] / np.maximum(P_read[sel][:, i5] + P_read[sel][:, ih], 1e-9)
                r['pair_5mC_vs_5hmC'] = {'n_reads': int(sel.sum()), 'read_auroc': round(float(roc_auc_score(yr[sel] == ih, pp)), 4),
                                         'read_acc': round(float(((pp > 0.5) == (yr[sel] == ih)).mean()), 4)}
        return r

    def detection(det, y):
        t = (y != known.index('none')).astype(int) if 'none' in known else None
        if t is None or t.min() == t.max():
            return None
        return round(float(roc_auc_score(t, det)), 4)

    # before training
    Pr, ir, Pi, det0, _ = evaluate(test_idx)
    result['detection_auroc_test_before'] = detection(det0, y_img[test_idx])
    print(f'  detection AUROC on test images, {os.path.basename(ckpt)} as released: {result["detection_auroc_test_before"]}', flush=True)

    # training
    head.train(); hist = []
    for ep in range(a.epochs):
        perm = rng.permutation(train_idx); tot, nb = 0.0, 0
        for s in range(0, len(perm), a.batch):
            b = perm[s:s + a.batch]; xb = X[b].to(dev); yb = torch.as_tensor(y_img[b], device=dev)
            logit, enc = forward(xb, a.mode == 'finetune'); rd = read_mask(xb)
            loss = F.cross_entropy(head(enc[rd]), yb.unsqueeze(1).expand(-1, xb.shape[2])[rd], weight=cw)
            if a.mode == 'finetune' and a.detect_weight > 0 and 'none' in known:
                loss = loss + a.detect_weight * F.binary_cross_entropy_with_logits(logit, (yb != known.index('none')).float())
            opt.zero_grad(); loss.backward(); opt.step(); tot += loss.item(); nb += 1
        Pr, ir, Pi, det, _ = evaluate(test_idx); m = metrics(Pr, ir, Pi, y_img[test_idx])
        hist.append({'epoch': ep + 1, 'train_loss': round(tot / max(nb, 1), 4), 'test_read_macro_f1': m['read_macro_f1'], 'test_image_acc': m['image_acc']})
        print(f'  epoch {ep + 1}/{a.epochs}  loss {tot / max(nb, 1):.4f}  test read macro-F1 {m["read_macro_f1"]}  image acc {m["image_acc"]}', flush=True)
    result['history'] = hist

    # final evaluation
    Pr, ir, Pi, det, Fr = evaluate(test_idx)
    result['test'] = metrics(Pr, ir, Pi, y_img[test_idx])
    result['detection_auroc_test_after'] = detection(det, y_img[test_idx])
    print(f'  final: read macro-F1 {result["test"]["read_macro_f1"]}  per-class read recall {result["test"]["read_recall"]}  '
          f'image acc {result["test"]["image_acc"]}  detection AUROC {result["detection_auroc_test_before"]} -> {result["detection_auroc_test_after"]}', flush=True)

    # Mahalanobis distance to the nearest training-class centroid in the head's hidden feature (shared covariance)
    _, it_, _, _, Ft = evaluate(train_idx); yt = y_img[train_idx][it_]
    mu = {c: Ft[yt == c].mean(0) for c in range(len(known)) if (yt == c).any()}
    R = np.concatenate([Ft[yt == c] - mu[c] for c in mu]); icov = np.linalg.inv(R.T @ R / len(R) + 1e-3 * np.eye(R.shape[1]))

    def mahal(Fx):
        return np.min([np.einsum('ij,jk,ik->i', Fx - mu[c], icov, Fx - mu[c]) for c in mu], axis=0)

    def auroc_unknown(score_known, score_unknown):   # higher score = more unknown
        return round(float(roc_auc_score(np.r_[np.zeros(len(score_known)), np.ones(len(score_unknown))], np.r_[score_known, score_unknown])), 4)

    # open set: unknown chemistries vs known test reads (reads of the trained 'other' class are left out of the known side)
    if len(unk_idx):
        Pu, iu, Piu, _, Fu = evaluate(unk_idx)
        yr_test = y_img[test_idx][ir]; sel_k = yr_test != OTHER
        known_conf, unk_conf = Pr[sel_k].max(1), Pu.max(1)
        thr = np.percentile(known_conf, 5)                      # keeps 95% of known test reads
        dk, du = mahal(Fr[sel_k]), mahal(Fu); dthr = np.percentile(dk, 95)
        o = {'n_unknown_reads': int(len(unk_conf)), 'chems': unknown,
             'auroc_known_vs_unknown_reads': auroc_unknown(-known_conf, -unk_conf),
             'unknown_flagged_at_tpr95': round(float((unk_conf < thr).mean()), 4),
             'unknown_assigned_to': {known[c]: round(float((Pu.argmax(1) == c).mean()), 4) for c in range(len(known))},
             'mahalanobis': {'auroc_known_vs_unknown_reads': auroc_unknown(dk, du),
                             'unknown_flagged_at_tpr95': round(float((du > dthr).mean()), 4)},
             'image_level': {'auroc_known_vs_unknown': auroc_unknown(-Pi[y_img[test_idx] != OTHER].max(1), -Piu.max(1))}}
        if OTHER >= 0:
            po_k, po_u = Pr[sel_k][:, OTHER], Pu[:, OTHER]
            o['p_other'] = {'auroc_known_vs_unknown_reads': auroc_unknown(po_k, po_u),
                            'unknown_called_other': round(float((Pu.argmax(1) == OTHER).mean()), 4),
                            'known_called_other': round(float((Pr[sel_k].argmax(1) == OTHER).mean()), 4),
                            'trained_other_recall_on_test': result['test']['read_recall'].get('other')}
        result['open_set'] = o
        print(f'  open set ({",".join(unknown)}): AUROC known-vs-unknown  max-softmax {o["auroc_known_vs_unknown_reads"]}  '
              f'mahalanobis {o["mahalanobis"]["auroc_known_vs_unknown_reads"]}' + (f'  P(other) {o["p_other"]["auroc_known_vs_unknown_reads"]}' if OTHER >= 0 else '')
              + f'  assigned to {o["unknown_assigned_to"]}', flush=True)

    # held-out files (never trained on)
    if ev:
        result['eval_files'] = {}
        for c, (Xe, yi) in ev.items():
            Pe, ie, Pie, dete, Fe = evaluate(Xe)
            e = {'n_images': int(len(Xe)), 'n_reads': int(len(Pe))}
            if yi >= 0:
                e['read_recall'] = round(float((Pe.argmax(1) == yi).mean()), 4); e['image_recall'] = round(float((Pie.argmax(1) == yi).mean()), 4)
            e['assigned_to'] = {known[k]: round(float((Pe.argmax(1) == k).mean()), 4) for k in range(len(known))}
            e['mean_max_softmax'] = round(float(Pe.max(1).mean()), 4)
            e['mean_mahalanobis'] = round(float(mahal(Fe).mean()), 4)
            if OTHER >= 0:
                e['called_other'] = round(float((Pe.argmax(1) == OTHER).mean()), 4)
            result['eval_files'][c] = e
            print(f'  eval {c}: {e}', flush=True)

    hook.remove()
    torch.save({'head': head.state_dict(), 'classes': known, 'other_chems': other_chems, 'checkpoint': ckpt, 'mode': a.mode,
                'mahalanobis': {'mu': {known[c]: mu[c] for c in mu}, 'icov': icov},
                'model': model.state_dict() if a.mode == 'finetune' else None}, os.path.join(a.out, 'typing_head.pt'))
    finish('OK' if not result['warnings'] else 'OK_WITH_WARNINGS')
except Exception as ex:
    import traceback; traceback.print_exc(); result['error'] = repr(ex); finish(f'FAILED: {ex!r}'); sys.exit(1)
