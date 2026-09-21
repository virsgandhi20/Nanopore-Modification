#!/usr/bin/env python3
"""Figure: is modification type decodable from frozen RawMod checkpoints? (reads P01 metrics)"""
import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "../../benchmark_results/typing_followup")
P = json.load(open(D + "/P01_rawmod_probe.json"))["checkpoints"]
BLUE, ORANGE, AQUA, GREY, INK, INK2, SURF = "#2a78d6", "#eb6834", "#1baf7a", "#a9a8a1", "#0b0b0b", "#52514e", "#fcfcfb"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10.5, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#c9c8c2",
                     "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK, "axes.titlesize": 12, "axes.titleweight": "bold",
                     "axes.titlelocation": "left", "figure.facecolor": SURF, "axes.facecolor": SURF})
fig, (a, b) = plt.subplots(1, 2, figsize=(14.5, 5.6), gridspec_kw={"width_ratios": [1.55, 1]}); fig.subplots_adjust(left=0.06, right=0.98, top=0.80, bottom=0.17, wspace=0.22)
fig.suptitle("Frozen RawMod checkpoints + a linear classifier: modification type on the ONT oligos (none / 5mC / 5hmC / 6mA), test sites unseen",
             x=0.06, y=0.96, ha="left", fontsize=12, color=INK2)
cks = [("mixed", "mixed\n(trained on all)"), ("loco_5mC", "never saw\n5mC"), ("loco_5hmC", "never saw\n5hmC"), ("loco_6mA", "never saw\n6mA")]
series = [("raw signal + dwell per read (no model)", GREY, "read_raw_signal_baseline"), ("per read, read encoder only", AQUA, "read_pre_transformer"),
          ("per read, after cross-read Transformer", BLUE, "read_post_transformer"), ("pooled site representation (96-d)", ORANGE, "site_rep_96d")]
w = 0.19
for k, (lab, col, key) in enumerate(series):
    v = [P[c][key]["macro_f1"] for c, _ in cks]; x = np.arange(len(cks)) + (k - 1.5) * (w + 0.015)
    a.bar(x, v, w, color=col, label=lab)
    for xi, vi in zip(x, v): a.text(xi, vi + 0.012, f"{vi:.2f}", ha="center", fontsize=9.5)
a.axhline(0.25, color=INK2, lw=1, ls=(0, (4, 3))); a.text(len(cks) - 0.42, 0.262, "chance", fontsize=9, color=INK2)
a.set_xticks(range(len(cks))); a.set_xticklabels([l for _, l in cks]); a.set_ylim(0, 1.0); a.set_ylabel("macro-F1, four classes")
a.yaxis.grid(True, color="#e6e5e0", lw=0.8); a.set_axisbelow(True); a.legend(frameon=False, loc="upper left", ncol=2, fontsize=9.5)
a.set_title("A. A detection-only model already encodes type, and most of it\n    is added where reads attend to each other")
cm = np.array(P["loco_5hmC"]["site_rep_96d"]["confusion_rows_true_cols_pred"], float); cm /= cm.sum(1, keepdims=True); cls = ["none", "5mC", "5hmC", "6mA"]
b.imshow(cm, cmap=matplotlib.colors.LinearSegmentedColormap.from_list("b", ["#eef4fd", "#86b6ef", "#2a78d6", "#0d366b"]), vmin=0, vmax=1)
for i in range(4):
    for j in range(4): b.text(j, i, f"{100*cm[i, j]:.0f}%", ha="center", va="center", fontsize=12, color="white" if cm[i, j] > 0.45 else INK)
b.set_xticks(range(4)); b.set_yticks(range(4)); b.set_xticklabels(cls); b.set_yticklabels(cls); b.set_xlabel("called as (linear probe)"); b.set_ylabel("true chemistry of the site")
for s in b.spines.values(): s.set_visible(False)
b.set_title("B. Checkpoint that never saw 5hmC:\n    errors stay within the base", x=-0.12)
out = os.path.join(HERE, "../../figures/typing"); os.makedirs(out, exist_ok=True); out = os.path.join(out, "rawmod_embedding_probe_2026-09-21.png"); fig.savefig(out, dpi=170); print(os.path.normpath(out))
