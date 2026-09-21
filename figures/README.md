# Figures

Every figure here is rebuilt by a script in this repository from recorded numbers or metric files, so none of
them needs cluster access to regenerate. Two folders, matching the two lines of work.

## `paper_benchmarking/` - joint work for the RawMod paper (Bhargav Srinivasan, Ernest Zhang, Vir Gandhi)

| figure | what it shows | built by | numbers from |
|---|---|---|---|
| `unimeth_update_2026-09-20.png` | (A) UniMeth's 5 kHz normalization bug, before and after the `--frequency 4khz` workaround, same reads and checkpoint; (B) the eight UniMeth bacteria rows under both site scores; (C) the T. denticola and H. pylori J99 motif presets are not methylated according to two independent callers | `scripts/benchmark/make_unimeth_figures.py` | `benchmark_results/unimeth/` |

The Table 1 columns themselves (UniMeth, Dorado, DeepMod2) are tables, not figures, because every value is 0.99 or above:
`benchmark_results/unimeth/TABLE1_UNIMETH_2026-09-20.md`, `benchmark_results/dorado/TABLE1_DORADO_2026-09-21.md`,
`benchmark_results/deepmod2/TABLE1_DEEPMOD2_2026-09-21.md`. The live shared tracker is the Google Sheet
"RawMod Table 1 benchmark tracker".

## `typing/` - Vir's project: from detecting a modification to identifying it

| figure | what it shows | built by | numbers from |
|---|---|---|---|
| `typing_first_results_2026-09-21.png` | per-read typing with a small CNN: (A) signal versus sequence-only controls; (B) 5mC / 5hmC per-read confusion; (C) recall by whether the sequence context was seen, oligos and four organisms; (D) what a closed-set model calls a chemistry it never saw | `scripts/typing/make_typing_figure.py` | `benchmark_results/typing_overnight/*.json` |
| `rawmod_embedding_probe_2026-09-21.png` | modification type decoded by a linear classifier from frozen RawMod checkpoints, by checkpoint and by tap point, with the confusion matrix of the checkpoint that never saw 5hmC | `scripts/typing/make_probe_figure.py` | `benchmark_results/typing_followup/P01_rawmod_probe.json` |

Naming: `<topic>_<YYYY-MM-DD>.png`, the date being the day the underlying results were produced. A figure that is
superseded keeps its file; the new one gets a new date.
