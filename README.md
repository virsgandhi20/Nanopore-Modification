# RawMod

RawMod detects DNA base modifications from raw nanopore signal without being
told in advance which modification to look for. Reads are stacked into a
per-position pileup image and classified by a per-read convolutional trunk
followed by a cross-read transformer (`ConvFormerV2`, ~209K parameters),
trained with a paired supervised-contrastive objective on top of binary
cross-entropy.

## Pretrained checkpoint

`checkpoints/results20_sad_dim16/mixed.pt` is the general-purpose checkpoint
(trained on all five chemistries) — use it to score new data. The
`loco_<CHEM>.pt` checkpoints hold out `CHEM` from training and exist only to
reproduce the zero-shot numbers below; don't use them to score `CHEM` sites
in new data.

| Fold | AUROC | AUPRC | macro F1 | n_pos / n_test |
|---|---:|---:|---:|---:|
| `mixed` (in-distribution) | 0.995 | 0.999 | 0.958 | 39,028 / 51,423 |
| `loco_5mC` (zero-shot) | 0.868 | 0.687 | 0.780 | 1,174 / 4,049 |
| `loco_5hmC` (zero-shot) | 0.875 | 0.582 | 0.750 | 621 / 3,496 |
| `loco_6mA` (zero-shot) | 0.818 | 0.857 | 0.716 | 11,716 / 19,025 |
| `loco_4mC` (zero-shot) | 0.648 | 0.808 | 0.400 | 4,780 / 7,181 |
| `loco_5hmU` (zero-shot) | 0.761 | 0.929 | 0.576 | 4,658 / 5,707 |

![Zero-shot AUROC by held-out chemistry](docs/figures/loco_results20_sad_dim16.png)

## Overview

A modification classifier trained on motif-derived labels alone can simply
memorize the recognition motif instead of reading the signal. RawMod trains
only on samples with a matched, genuinely unmodified counterpart at the same
coordinate (synthetic control, amplicon-stripped, or whole-genome-amplified
DNA), and evaluates with leave-one-chemistry-out / leave-one-organism-out
splits, so a reported score reflects generalization rather than memorization.

Data artifacts (POD5, BAM, `features.h5`, checkpoints) are not stored in this
repo; scripts reference them by absolute path on shared scratch storage (see
"Paths").

## Repository layout

| Path | Contents |
|---|---|
| `rawmod/` | `featurization.py` builds pileup tensors from POD5/BAM/peaks; `model.py` provides `PileupDataset`, data splits, and evaluation. |
| `scripts/ground_truth/` | Ground-truth extraction (motif-based, bisulfite/EM-seq, pileup-derived) and background-site generation. |
| `scripts/featurize/` | `refeaturize_strand15.py` (matched pool), `refeaturize_benchmark.py` (benchmark organisms + human), `featurize_background.py` (background negatives). |
| `scripts/train/` | `run_pipeline.py` (training loop), `run_convformer_v2.py` (model), `run_matched_loco.py` / `.sh` (training + evaluation entry point). |
| `scripts/test/` | `score_genome.py` and `test_external_sites.py` — score a checkpoint against data outside the built-in folds. |
| `analysis/chem_diversity_sweep/` | Ablation: how many training chemistries are needed for generalization. |
| `analysis/orca_remake/` | Embedding/clustering diagnostics on a trained model. |
| `analysis/visualize_h5_pileup.py` | Render a `features.h5` pileup as a figure. |

## Paths

```
/fs/cbcb-lab/storm/bds062/data/benchmark/                       POD5 and references, benchmark organisms
/fs/cbcb-scratch/bds062/data/human/{hg001,hg002}/pod5/           POD5, human
/fs/cbcb-scratch/bds062/data/gt/                                 ground-truth BED files, all organisms
/fs/cbcb-scratch/bds062/results/benchmark_results/               reads_refined.bam / peaks_refined.tsv
/fs/cbcb-scratch/bds062/results/rawmod_full_pipeline4/features/  features.h5 (matched pool, benchmark organisms, background sites)
/fs/cbcb-scratch/bds062/results/rawmod_matched_loco/             run_matched_loco.py output: models/, metrics/
```

## Installation

```bash
git clone <this repository>
cd RawMod
pip install -e .
```

Requires Python 3.10+; dependencies (`numpy`, `h5py`, `pod5`, `pysam`,
`torch`, `scikit-learn`, `matplotlib`, `tqdm`) are in `pyproject.toml`.

## Reproducing the pipeline

### 1. Ground truth

```bash
python scripts/ground_truth/motif_gt.py --ref REF.fa.gz --preset ecoli_dam --outdir data/gt/Ecoli_DM
python scripts/ground_truth/extract_gt_bismark.py ...     # EM-seq / WGBS
python scripts/ground_truth/extract_gt_from_pileup.py ... # pre-computed bedMethyl
python scripts/ground_truth/generate_background_sites.py  # non-motif negatives, used by logo_bacteria
```

### 2. Basecalling and refinement

```bash
bash scripts/ground_truth/submit_all.sh   # per-dataset table; invokes pipeline.sh
```

Runs Dorado basecalling and Remora move refinement into `reads_refined.bam` /
`peaks_refined.tsv`. Runs once; every featurization script below reuses it.

### 3. Featurization

```bash
python scripts/featurize/refeaturize_strand15.py     # matched pool: ONT, SPO1/UMCES, HP26695
python scripts/featurize/refeaturize_benchmark.py     # benchmark organisms + hg001/hg002 (curriculum only)
python scripts/featurize/featurize_background.py      # background negatives for logo_bacteria
```

Each `--dry-run` first. Output: `(16, 210, 9)` tensors — a reference row plus
15 read rows, 21 window positions at 10 samples/base, 9 channels
(`raw_signal, dwell_log1p, is_A, is_C, is_G, is_T, strand, mapq_norm,
matches_ref`), forward-strand-only (`--strand +`).

### 4. Training and evaluation

Each fold trains and evaluates in one job, writing `metrics/<fold>.tsv`.

```bash
FOLDS="mixed loco_5hmU loco_4mC loco_6mA loco_5mC loco_5hmC logo_bacteria logo_plant logo_mammal" \
OUTDIR=/fs/cbcb-scratch/bds062/results/rawmod_matched_loco/<results_dir> \
RAWMOD_DATA_GEN=strand15 EXTRA_ORGANISMS=1 INCLUDE_HUMAN=1 \
SUPCON_DIM=128 SUPCON_WEIGHT=1.0 SUPCON_TEMP=0.20 \
CURRICULUM=1 CURRICULUM_EPOCHS=15 \
SAD_DIM=16 SAD_WEIGHT=1.0 SAD_ETA=1.0 BCE_WEIGHT=1.0 \
bash scripts/train/run_matched_loco.sh
```

This is the recipe behind `checkpoints/results20_sad_dim16/`. Each fold runs
as an independent SLURM job (~7h/fold, RTX A5000); `TIME_LIMIT`, `PARTITION`,
`GPU_TYPE` override the defaults — see the script header. `--dry-run`
inspects the generated `sbatch` commands without submitting.

**Folds:**

- `mixed` — position-grouped 85/15 split over the whole matched pool.
- `loco_<CHEM>`, `CHEM` in `5hmU, 4mC, 6mA, 5mC, 5hmC` — leave-one-chemistry-out.
- `logo_<group>`, `group` in `bacteria, plant, mammal` — leave-one-organism-group-out.
- `subset_<c1>+<c2>[+c3]` — training-diversity sweep; see `analysis/chem_diversity_sweep/`.

Metrics columns (`metrics/<fold>.tsv`): `micro_f1, mod_f1, unmod_f1,
macro_f1, mod_prec, mod_rec, auprc, auroc, auroc_sad, threshold, n_pos,
n_test`. `auroc` is threshold-free and comparable across folds; `mod_f1` /
`macro_f1` use a per-fold optimal threshold, so compare those only within a
fold, not across folds.

### 5. Scoring external data

**From a site list** (contig/position pairs — a candidate BED, or another
tool's calls):

```bash
python scripts/test/test_external_sites.py \
  --sites <tsv with contig, pos columns> \
  --pod5 <pod5 dir> --bam <reads_refined.bam> --peaks <peaks_refined.tsv> \
  --gt <ground-truth BED> \
  --checkpoint checkpoints/results20_sad_dim16/mixed.pt \
  --out-dir <output dir>
```

Ground truth is always looked up from `--gt`, never from a label column the
input file may already carry.

**From an existing `features.h5`:**

```bash
python scripts/test/score_genome.py \
  --h5 <features.h5> --dataset <name> \
  --checkpoint checkpoints/results20_sad_dim16/mixed.pt \
  --out-dir <output dir>
```

`score_genome.py`'s default `--checkpoint` picks the leave-one-dataset-out
fold for `--dataset`, for non-circularity when that organism was in the
training pool; pass `--checkpoint` explicitly to use `mixed` instead.

## Analysis

```bash
python analysis/visualize_h5_pileup.py --h5 features.h5 --cartoon
bash analysis/chem_diversity_sweep/run_chem_diversity_sweep.sh
python analysis/orca_remake/loco_embedding_cluster.py --help
```

## Notes

**Strand.** Pileups are forward-strand-only (`--strand +`) by construction —
pooling both strands requires the reference row and read rows to agree on
orientation, which an earlier pipeline got wrong. All shipped checkpoints
were trained with `--strand +` and are unaffected.

**Curriculum chemistry overlap.** Benchmark/human organisms unioned into
every fold's curriculum data (`EXTRA_ORGANISMS=1`, `INCLUDE_HUMAN=1`) can
carry the same chemistry as a `loco_<CHEM>` holdout under an unrecognized
label. `BENCH_ORG_CHEMS` in `run_matched_loco.py` excludes any such organism
from the affected fold's training. Not applied to `subset_<...>` folds.
