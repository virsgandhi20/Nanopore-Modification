#!/usr/bin/env python3
"""Summary figure for the UniMeth benchmarking update (2026-09-20).

All numbers are copied from the run logs recorded in
benchmark_results/unimeth/ (diag_summary.tsv, TABLE1_UNIMETH_2026-09-20.md,
gt_check_*.txt), so the figure can be rebuilt without cluster access.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../../benchmark_results/unimeth/figures")
os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, GREY, RED = "#2C6FB7", "#E8892B", "#9AA0A6", "#C0392B"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.titlesize": 12.5, "axes.titleweight": "bold"})

fig = plt.figure(figsize=(15, 9.6))
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.25], hspace=0.62, wspace=0.42)
ax = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[:, 1]), fig.add_subplot(gs[1, 0])]

# ---- A. the normalization bug: same reads, same checkpoint
a = ax[0]
x = np.arange(2); w = 0.36
gt, bg = [0.0009, 0.8618], [0.0010, 0.0112]
b1 = a.bar(x - w/2, gt, w, color=BLUE, label="Dam GATC sites (methylated)")
b2 = a.bar(x + w/2, bg, w, color=GREY, label="all other adenines")
for bars, vals in ((b1, gt), (b2, bg)):
    for r, v in zip(bars, vals):
        a.text(r.get_x() + r.get_width()/2, v + 0.02, f"{v:.3f}" if v < 0.1 else f"{v:.2f}", ha="center", fontsize=10)
a.set_xticks(x); a.set_xticklabels(["UniMeth 5 kHz path\n(pA tags on raw DAC)", "calibrated input\n(--frequency 4khz)"])
a.set_ylabel("mean per-read P(6mA)"); a.set_ylim(0, 1.08)
a.text(0, 0.10, "per-read AUROC 0.48", ha="center", color=RED, fontweight="bold")
a.text(1, 0.95, "per-read AUROC 0.985", ha="center", color=BLUE, fontweight="bold")
a.set_title("A. UniMeth was blind on our data; one flag fixes it\n(E. coli, same reads, same checkpoint)", loc="left")
a.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.84), fontsize=9.5)
a.set_xlim(-0.6, 1.6)

# ---- B. Table 1 rows
b = ax[1]
rows = [("E. coli 6mA (Dam)", 0.9998, 0.9998), ("Anabaena 6mA", 0.9998, 0.9999), ("H. pylori 26695 6mA", 0.9988, 0.9989),
        ("E. coli M.SssI 5mC CpG", 1.0000, 1.0000), ("E. coli Dcm 5mC non-CpG", 0.7244, 0.9992),
        ("H. pylori 26695 5mC", 0.8460, 0.9290), ("T. denticola 6mA  (preset GT)", 0.4905, 0.5540),
        ("H. pylori J99 6mA  (preset GT)", 0.4993, 0.5224)]
y = np.arange(len(rows))[::-1]; h = 0.38
f = [r[1] for r in rows]; m = [r[2] for r in rows]
b.barh(y + h/2, f, h, color=ORANGE, label="UniMeth call frequency (P > 0.5)")
b.barh(y - h/2, m, h, color=BLUE, label="mean P(mod) per site")
for yi, fv, mv in zip(y, f, m):
    b.text(fv + 0.006, yi + h/2, f"{fv:.4f}", va="center", fontsize=9)
    b.text(mv + 0.006, yi - h/2, f"{mv:.4f}", va="center", fontsize=9)
b.axvline(0.5, color=RED, ls="--", lw=1); b.text(0.503, y[0] + 0.62, "chance", color=RED, fontsize=9)
b.axhspan(-0.6, 1.5, color=RED, alpha=0.06)
b.text(0.63, 0.45, "these two rows score the ground truth,\nnot the caller (panel C)", color=RED, fontsize=10, va="center")
b.set_yticks(y); b.set_yticklabels([r[0] for r in rows]); b.set_xlim(0.4, 1.09); b.set_ylim(-0.6, len(rows) - 0.3)
b.set_xlabel("site-level AUROC"); b.set_title("B. UniMeth rows for Table 1 (bacteria)", loc="left")
b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.45, -0.07), ncol=2, fontsize=10)

# ---- C. preset ground truth vs two independent callers
c = ax[2]
orgs = ["E. coli (control)\npreset: GATC\nfound: G[A]TC", "T. denticola\npreset: GATC, TATAC\nfound: RA[A]TTY,\nCTA[A]T, GAAG[A]G",
        "H. pylori J99\npreset: GTNNNNNNAC\nfound: C[A]TG, G[A]TC,\nG[A]GG, G[A]NTC, ..."]
uni, dor = [93.2, 0.01, 2.4], [98.8, 0.02, 2.7]
lab_u, lab_d = ["93.2%", "1 of 9,766", "2.4%"], ["98.8%", "3 of 15,368", "2.7%"]
x = np.arange(3)
c1 = c.bar(x - w/2, uni, w, color=BLUE, label="UniMeth")
c2 = c.bar(x + w/2, dor, w, color=ORANGE, label="Dorado 6mA (modkit)")
for r, t in zip(c1, lab_u): c.text(r.get_x() + r.get_width()/2, r.get_height() + 2, t, ha="center", fontsize=9.5, rotation=0 if "of" not in t else 90, va="bottom")
for r, t in zip(c2, lab_d): c.text(r.get_x() + r.get_width()/2, r.get_height() + 2, t, ha="center", fontsize=9.5, rotation=0 if "of" not in t else 90, va="bottom")
c.set_xticks(x); c.set_xticklabels(orgs, fontsize=9.5); c.set_ylim(0, 118); c.set_ylabel("motif-preset sites called methylated (%)")
c.set_title("C. Two motif presets are not methylated in the data\n(two unrelated callers agree)", loc="left")
c.legend(frameon=False, loc="upper right", fontsize=9.5)

fig.subplots_adjust(left=0.06, right=0.97, top=0.92, bottom=0.13)
p = os.path.join(OUT, "unimeth_update_2026-09-20.png")
fig.savefig(p, dpi=200, facecolor="white"); print(os.path.normpath(p))
