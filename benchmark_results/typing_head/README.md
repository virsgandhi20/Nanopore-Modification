# Typing head on RawMod's encoder: results (Oct 28, 2026)

One metrics file per run (T01.json ...), written by scripts/typing/rawmod_typing_head.py and copied here by
`typing_head.sh status`. Grid and notes: scripts/typing/typing_head.sh. Data: ONT synthetic oligos, strand-resolved
features (Bhargav's Sep 21 set), up to 3,000 images per class, split 80/20 by site, six epochs, one GPU job each.

| Run | Start checkpoint | Encoder | Head classes | Read macro-F1 | 5mC vs 5hmC AUROC | Image acc | Detection AUROC before -> after | Unknown: AUROC, flagged at TPR95, called as |
|---|---|---|---|---|---|---|---|---|
| T01 | region_holdout/mixed.pt | frozen | none,5mC,5hmC,6mA | 0.695 | 0.651 | 0.742 | 0.959 -> 0.959 | - |
| T02 | region_holdout/mixed.pt | fine-tuned (+BCE guard) | none,5mC,5hmC,6mA | 0.712 | 0.635 | 0.755 | 0.959 -> 0.982 | - |
| T03 | region_holdout/mixed.pt | fine-tuned | none,5mC,5hmC,6mA | not run (pipeline features pending) | | | | |
| T04 | leaveout/ont5hmC.pt | fine-tuned | none,5mC,6mA | 0.896 | - | 0.902 | 0.934 -> 0.974 | 0.323, 3.7%, 5mC 89% |
| T05 | leaveout/ont6mA.pt | fine-tuned | none,5mC,5hmC | 0.633 | 0.669 | 0.740 | 0.977 -> 0.993 | 0.301, 1.8%, none 94% |
| T06 | region_holdout/loco_5hmC.pt | fine-tuned | none,5mC,6mA | 0.930 | - | 0.935 | 0.960 -> 0.984 | 0.386, 6.4%, 5mC 89% |
| T07 | rawmod.pt (all data) | frozen | none,5mC,5hmC,6mA | 0.660 | 0.617 | 0.736 | 0.988 -> 0.988 | - |
| T08 | region_holdout/mixed.pt | fine-tuned, no BCE guard | none,5mC,5hmC,6mA | 0.728 | 0.672 | 0.760 | 0.959 -> 0.978 | - |

Reference points on the same oligos: approach 1 (window model, Sep 21) 0.61 four-class macro-F1, 0.91 pair AUROC
(dedicated two-class model), 0.76 on replicate 2; linear probe on frozen mixed.pt 0.67 (0.75 / 0.72 on the loco
checkpoints).

Reading:
- The head beats approach 1 overall (0.69 to 0.73 vs 0.61) and matches or beats the linear probe. Per-class recall:
  none and 6mA 0.85 to 0.94 in every run; 5mC and 5hmC cross into each other (T02: 0.39 / 0.70; T01: 0.55 / 0.50).
- The head loses the hard pair (AUROC 0.62 to 0.67 vs 0.91): the detection representation keeps "modified C" and
  discards "which C". Caveat: approach 1's 0.91 is a dedicated binary model; the head's number is the four-class head
  restricted to the pair's reads.
- Detection never fell; fine-tuning for typing raised it with or without the BCE guard (T02 0.982, T08 0.978), so the
  guard was not needed on this data.
- Open set fails: unknown-vs-known AUROC 0.30 to 0.39 (below chance), 2 to 6% of unknown reads flagged at TPR95.
  Unseen 5hmC is called 5mC 89% of the time (T04, T06); unseen 6mA is called none 94% of the time (T05), even from a
  checkpoint that never saw 6mA. Confidence thresholding does not work; an explicit "other" class has to be trained.
- All-data rawmod.pt (T07) holds slightly less type than the region-holdout model (0.66 vs 0.69) and collapses 5hmC
  into 5mC (5hmC recall 0.13).
