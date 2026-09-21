#!/usr/bin/env python3
"""Summary figure of the first per-read typing results (reads the metric JSONs)."""
import glob, json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "../../benchmark_results/typing_overnight")
M = {os.path.basename(f)[:-5]: json.load(open(f)) for f in glob.glob(D + "/E*.json")}
T = lambda e, s="test": M[e]["sets"][s]
BLUE, ORANGE, AQUA, GREY, INK, INK2, SURF = "#2a78d6", "#eb6834", "#1baf7a", "#a9a8a1", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#c9c8c2", "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left", "figure.facecolor": SURF, "axes.facecolor": SURF})
fig = plt.figure(figsize=(15.5, 10.4))
ax = np.array([[fig.add_axes([0.06, 0.56, 0.40, 0.30]), fig.add_axes([0.62, 0.56, 0.22, 0.30])],
               [fig.add_axes([0.285, 0.08, 0.215, 0.33]), fig.add_axes([0.645, 0.08, 0.33, 0.33])]])
fig.suptitle("Per-read modification typing, first results (1D CNN on 21-base signal windows; test sites never seen in training; one seed)",
             x=0.06, y=0.965, ha="left", fontsize=12.5, color=INK2)

# A ------------------------------------------------------------------ signal vs sequence
a = ax[0, 0]; groups = ["ONT oligos\n(none / 5mC / 5hmC / 6mA)", "E. coli, matched controls\n(none / 6mA / 5mC)"]
series = [("signal + dwell + sequence", BLUE, ["E01_syn_type4_all", "E07_ecoli_type3_all"]), ("signal + dwell", AQUA, ["E02_syn_type4_signal", "E08_ecoli_type3_signal"]),
          ("sequence only (control)", GREY, ["E03_syn_type4_seqonly", "E09_ecoli_type3_seqonly"])]
w = 0.24
for k, (lab, col, exps) in enumerate(series):
    v = [T(e)["read_macro_f1"] for e in exps]; x = np.arange(2) + (k - 1) * (w + 0.02)
    a.bar(x, v, w, color=col, label=lab)
    for xi, vi in zip(x, v): a.text(xi, vi + 0.015, f"{vi:.2f}", ha="center", fontsize=10)
for gi, ch in enumerate((0.25, 1 / 3)):
    a.plot([gi - 0.42, gi + 0.42], [ch, ch], color=INK2, lw=1, ls=(0, (4, 3))); a.text(gi + 0.43, ch, "chance", va="center", fontsize=9, color=INK2)
a.set_xticks(range(2)); a.set_xticklabels(groups); a.set_ylim(0, 1.0); a.set_ylabel("per-read macro-F1"); a.set_xlim(-0.55, 1.75)
a.yaxis.grid(True, color="#e6e5e0", lw=0.8); a.set_axisbelow(True); a.legend(frameon=False, loc="upper left", fontsize=9.5)
a.set_title("A. Typing is read from the signal: a sequence-only model is at chance\n    on the oligos and far behind on E. coli")

# B ------------------------------------------------------------------ 5mC vs 5hmC confusion
b = ax[0, 1]; cm = np.array(T("E04_syn_5mC_vs_5hmC")["confusion_rows_true_cols_pred"]["matrix"], float); cm = cm / cm.sum(1, keepdims=True)
cls = T("E04_syn_5mC_vs_5hmC")["confusion_rows_true_cols_pred"]["classes"]
b.imshow(cm, cmap=matplotlib.colors.LinearSegmentedColormap.from_list("b", ["#eef4fd", "#86b6ef", "#2a78d6", "#0d366b"]), vmin=0, vmax=1)
for i in range(3):
    for j in range(3): b.text(j, i, f"{100*cm[i, j]:.0f}%", ha="center", va="center", fontsize=13, color="white" if cm[i, j] > 0.45 else INK)
b.set_xticks(range(3)); b.set_yticks(range(3)); b.set_xticklabels(cls); b.set_yticklabels(cls); b.set_xlabel("called as"); b.set_ylabel("true label (per read)")
for s in b.spines.values(): s.set_visible(False)
b.set_title(f"B. 5mC vs 5hmC: separable per read\n    (site accuracy after pooling {T('E04_syn_5mC_vs_5hmC')['site_acc']:.2f}, chance 0.33)")

# C ------------------------------------------------------------------ context generalization
c = ax[1, 0]; e5 = M["E05_syn_cross_replicate"]["sets"]; e10 = M["E10_bact_type4_all"]["sets"]
seen = np.mean([e5[f"test/syn_{x}_rep2"]["read_recall_by_class"][x] for x in ("5mC", "5hmC", "6mA")])
unseen = np.mean([e5[f"test/syn_{x}_rep1"]["read_recall_by_class"][x] for x in ("5mC", "5hmC", "6mA")])
rows = [("Oligos: other flow cell,\ncontexts seen in training", seen, BLUE), ("Oligos: same flow cell,\n5-mer contexts NOT seen", unseen, BLUE),
        ("E. coli 6mA: held-out sites,\norganism in training", e10["test/ecoli_wt"]["read_recall_by_class"]["6mA"], ORANGE),
        ("Anabaena 6mA: organism never\nseen, same GATC motif", e10["anabaena"]["read_recall_by_class"]["6mA"], ORANGE),
        ("H. pylori J99 6mA: never seen,\nmotifs partly shared", e10["j99"]["read_recall_by_class"]["6mA"], ORANGE),
        ("T. denticola 6mA: never seen,\nmotifs never seen", e10["tdent"]["read_recall_by_class"]["6mA"], ORANGE)]
y = np.arange(len(rows))[::-1]
c.barh(y, [r[1] for r in rows], 0.55, color=[r[2] for r in rows])
for yi, r in zip(y, rows): c.text(r[1] + 0.012, yi, f"{r[1]:.2f}", va="center", fontsize=10)
c.set_yticks(y); c.set_yticklabels([r[0] for r in rows], fontsize=9.5); c.set_xlim(0, 1.08); c.set_xlabel("per-read recall of the modified class")
c.xaxis.grid(True, color="#e6e5e0", lw=0.8); c.set_axisbelow(True)
c.legend(handles=[matplotlib.patches.Patch(color=BLUE, label="ONT oligos (mean of 5mC, 5hmC, 6mA)"), matplotlib.patches.Patch(color=ORANGE, label="bacteria, four-class model")],
         frameon=False, loc="upper center", bbox_to_anchor=(0.2, -0.13), ncol=2, fontsize=9.5)
c.set_title("C. What transfers is the sequence context,\n    not the chemistry", x=-0.95)

# D ------------------------------------------------------------------ open set
d = ax[1, 1]
o5 = T("E16_openset_5hmC")["openset"]; o4 = T("E15_openset_4mC")["openset"]
bars = [("5hmC reads, model never\ntrained on 5hmC", o5, ["none", "5mC", "6mA"]), ("4mC reads, model never\ntrained on 4mC", o4, ["none", "5mC", "6mA"])]
colr = {"none": GREY, "5mC": BLUE, "6mA": ORANGE}
for yi, (lab, o, order) in zip((1, 0), bars):
    left = 0
    for k in order:
        v = o["unknown_called_as"][k]; d.barh(yi, v, 0.5, left=left, color=colr[k], edgecolor=SURF, linewidth=2)
        if v > 0.12: d.text(left + v / 2, yi, f"{k}\n{100*v:.0f}%", ha="center", va="center", fontsize=10, color="white" if k != "none" else INK)
        left += v
    d.text(1.02, yi, f"flagged as unknown by\nlow confidence: AUROC {o['auroc_unknown_by_low_confidence']:.2f}", va="center", fontsize=9.5, color=INK2)
d.set_yticks([1, 0]); d.set_yticklabels([b_[0] for b_ in bars]); d.set_xlim(0, 1.48); d.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); d.set_xticklabels(["0", "25%", "50%", "75%", "100%"])
d.set_xlabel("what the closed-set model calls them"); d.set_ylim(-0.6, 1.6)
d.legend(handles=[matplotlib.patches.Patch(color=colr[k], label=f"called {k}") for k in ("none", "5mC", "6mA")], frameon=False, loc="upper center", bbox_to_anchor=(0.36, 1.0), ncol=3, fontsize=9.5)
d.set_title("D. A chemistry the model never saw is silently mislabelled\n    (flagging AUROC of 0.5 = cannot tell it is new)", x=-0.30)

out = os.path.join(HERE, "../../figures/typing", "typing_first_results_2026-09-21.png"); fig.savefig(out, dpi=170); print(os.path.normpath(out))
