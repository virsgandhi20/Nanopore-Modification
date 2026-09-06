# RawMod

Detects DNA base modifications from raw nanopore signal. Per-read pileup
images classified by `ConvFormerV2`.

## Pretrained checkpoint

`checkpoints/results20_sad_dim16/mixed.pt` — general-purpose, use this to
score new data. `loco_<CHEM>.pt` hold out `CHEM` from training (zero-shot
reproduction only; don't use to score `CHEM` in new data).

| Fold | AUROC | AUPRC | macro F1 | n_pos / n_test |
|---|---:|---:|---:|---:|
| `mixed` | 0.995 | 0.999 | 0.958 | 39,028 / 51,423 |
| `loco_5mC` | 0.868 | 0.687 | 0.780 | 1,174 / 4,049 |
| `loco_5hmC` | 0.875 | 0.582 | 0.750 | 621 / 3,496 |
| `loco_6mA` | 0.818 | 0.857 | 0.716 | 11,716 / 19,025 |
| `loco_4mC` | 0.648 | 0.808 | 0.400 | 4,780 / 7,181 |
| `loco_5hmU` | 0.761 | 0.929 | 0.576 | 4,658 / 5,707 |

![Zero-shot AUROC by held-out chemistry](docs/figures/loco_results20_sad_dim16.png)

## Paths

```
/fs/cbcb-lab/storm/bds062/data/benchmark/                       POD5 and references, benchmark organisms
/fs/cbcb-scratch/bds062/data/human/{hg001,hg002}/pod5/           POD5, human
/fs/cbcb-scratch/bds062/data/gt/                                 ground-truth BED files, all organisms
/fs/cbcb-scratch/bds062/results/benchmark_results/               reads_refined.bam / peaks_refined.tsv
/fs/cbcb-scratch/bds062/results/rawmod_full_pipeline4/features/  features.h5 (matched pool, benchmark organisms, background sites)
/fs/cbcb-scratch/bds062/results/rawmod_matched_loco/             run_matched_loco.py output: models/, metrics/
```

## Install

```bash
git clone <this repository>
cd RawMod
pip install -e .
```

## 1. Ground truth

```bash
python scripts/ground_truth/motif_gt.py --ref REF.fa.gz --preset ecoli_dam --outdir data/gt/Ecoli_DM
python scripts/ground_truth/extract_gt_bismark.py ...     # EM-seq / WGBS
python scripts/ground_truth/extract_gt_from_pileup.py ... # pre-computed bedMethyl
python scripts/ground_truth/generate_background_sites.py  # non-motif negatives, for logo_bacteria
```

## 2. Basecalling and refinement

```bash
bash scripts/ground_truth/submit_all.sh   # per-dataset table; invokes pipeline.sh
```

## 3. Featurization

```bash
python scripts/featurize/refeaturize_strand15.py --dry-run    # matched pool: ONT, SPO1/UMCES, HP26695
python scripts/featurize/refeaturize_strand15.py

python scripts/featurize/refeaturize_benchmark.py --dry-run   # benchmark organisms + hg001/hg002
python scripts/featurize/refeaturize_benchmark.py

python scripts/featurize/featurize_background.py --dry-run    # background negatives for logo_bacteria
python scripts/featurize/featurize_background.py
```

## 4. Training and evaluation

```bash
FOLDS="mixed loco_5hmU loco_4mC loco_6mA loco_5mC loco_5hmC logo_bacteria logo_plant logo_mammal" \
OUTDIR=/fs/cbcb-scratch/bds062/results/rawmod_matched_loco/<results_dir> \
RAWMOD_DATA_GEN=strand15 EXTRA_ORGANISMS=1 INCLUDE_HUMAN=1 \
SUPCON_DIM=128 SUPCON_WEIGHT=1.0 SUPCON_TEMP=0.20 \
CURRICULUM=1 CURRICULUM_EPOCHS=15 \
SAD_DIM=16 SAD_WEIGHT=1.0 SAD_ETA=1.0 BCE_WEIGHT=1.0 \
bash scripts/train/run_matched_loco.sh
```

Each fold trains and evaluates in one job, writing `metrics/<fold>.tsv`.
`--dry-run` inspects the generated `sbatch` commands without submitting.
`TIME_LIMIT`, `PARTITION`, `GPU_TYPE` override the defaults (see script
header).

## 5. Scoring external data

```bash
# from a site list (contig/position pairs)
python scripts/test/test_external_sites.py \
  --sites <tsv with contig, pos columns> \
  --pod5 <pod5 dir> --bam <reads_refined.bam> --peaks <peaks_refined.tsv> \
  --gt <ground-truth BED> \
  --checkpoint checkpoints/results20_sad_dim16/mixed.pt \
  --out-dir <output dir>

# from an existing features.h5
python scripts/test/score_genome.py \
  --h5 <features.h5> --dataset <name> \
  --checkpoint checkpoints/results20_sad_dim16/mixed.pt \
  --out-dir <output dir>
```

## Analysis

```bash
python analysis/visualize_h5_pileup.py --h5 features.h5 --cartoon
bash analysis/chem_diversity_sweep/run_chem_diversity_sweep.sh
python analysis/orca_remake/loco_embedding_cluster.py --help
```
