#!/bin/bash
# Two questions about the finished UniMeth Table 1 rows (CPU only):
#
# 1. E. coli Dcm scored 0.72 at site level although per-read separation is 0.98.
#    UniMeth's native site score is a call FREQUENCY at P > 0.5, and its 5mC
#    model is under-confident on Dcm (mean P 0.24 at true sites), so most true
#    sites get frequency 0 and tie with the negatives. Rescore every row with
#    mean P(mod) per site, which keeps the ranking.
#
# 2. H. pylori J99 and T. denticola scored 0.50, with HIGHER P(6mA) at
#    background than at the preset's motif sites, from the same pipeline that
#    gives 0.9998 on E. coli and Anabaena. Is the caller wrong, or the motif
#    ground truth? Find the methylated motifs de novo from UniMeth's calls, then
#    repeat with an independent caller (Dorado's 6mA modbam via modkit).
#    E. coli is the control: it must come back as G[A]TC.
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
OUT=/fs/nexus-scratch/vgandhi/unimeth_bench
BENCH=/fs/cbcb-lab/storm/bds062/data/benchmark
HPGT=/fs/nexus-scratch/vgandhi/hp_labels
MODKIT=/fs/nexus-scratch/vgandhi/dist_modkit_v0.6.4_cd85862/modkit
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
RES=$REPO/benchmark_results/unimeth
T=${SLURM_CPUS_PER_TASK:-8}
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null
mkdir -p $RES $OUT/gtcheck

declare -A REF=( [Ecoli_WT_5kHz]=ecoli.fa.gz [Ecoli_DM_MSssI_5kHz]=ecoli.fa.gz [Anabaena_WT_5kHz]=anabaena_sp_PCC7120_ATCC27893.fa.gz
                 [Tdenticola_WT_5kHz]=treponema_denticola_ATCC35405.fa.gz [HPJ99_WT_5kHz]=hpylori_J99_ATCC700824.fa.gz [HP26695_WT_5kHz]=hpylori_26695.fa.gz )
gtof() { case "$1" in HP26695_WT_5kHz.6mA) echo $HPGT/gt_6mA.bed;; HP26695_WT_5kHz.5mC) echo $HPGT/gt_5mC.bed;; *) echo $OUT/$1/gt/gt_modified.bed;; esac; }

# ---- 1. rescore with mean P(mod)
M=$OUT/table1_unimeth_meanprob.tsv; rm -f $M
for W in $OUT/*.6mA $OUT/*.5mC; do
    D=$(basename $W); [ -s $W/sites.tsv ] || continue
    python $REPO/scripts/benchmark/score_sites.py --calls $W/sites.tsv --gt $(gtof $D) --min-cov 10 \
        --chrom-col 0 --pos-col 1 --cov-col 8 --num-col 5 --label "${D/./\/}" --out $M > /dev/null 2>&1 || echo "rescore failed: $D"
done
cp $M $RES/ 2>/dev/null
echo "=== site-level, native call frequency ==="; column -t $OUT/table1_unimeth.tsv
echo; echo "=== site-level, mean P(mod) ==="; column -t $M

# ---- 2a. de novo motifs from UniMeth's own calls
{
for S in Ecoli_WT_5kHz Anabaena_WT_5kHz HP26695_WT_5kHz Tdenticola_WT_5kHz HPJ99_WT_5kHz; do
    python $REPO/scripts/benchmark/motif_enrichment.py --sites $OUT/$S.6mA/sites.tsv --ref $BENCH/references/${REF[$S]} \
        --gt $(gtof $S.6mA) --label "$S  UniMeth 6mA" --base A
    echo
done
} 2>&1 | tee $RES/gt_check_unimeth.txt

# ---- 2b. independent caller: Dorado 6mA modbam -> modkit pileup
{
for S in Ecoli_WT_5kHz Tdenticola_WT_5kHz HPJ99_WT_5kHz; do
    SRC=$BENCH/bacteria/$S/modbam/${S}_sup_v5r3_6mA.bam
    BED=$OUT/gtcheck/${S}_6mA.bed
    if [ ! -s $BED ]; then
        [ -s $SRC ] || { echo "==== $S: no Dorado 6mA modbam at $SRC"; continue; }
        ln -sf $SRC $OUT/gtcheck/${S}_6mA.bam      # share is read-only: index via a symlink
        [ -s $OUT/gtcheck/${S}_6mA.bam.bai ] || $SAM index -@ $T $OUT/gtcheck/${S}_6mA.bam
        $MODKIT pileup $OUT/gtcheck/${S}_6mA.bam $BED --threads $T 2> $OUT/gtcheck/${S}_6mA.log || { echo "==== $S: modkit failed: $(tail -1 $OUT/gtcheck/${S}_6mA.log)"; continue; }
    fi
    python $REPO/scripts/benchmark/motif_enrichment.py --sites $BED --format bedmethyl --code a --ref $BENCH/references/${REF[$S]} \
        --gt $(gtof $S.6mA) --label "$S  Dorado 6mA (modkit)" --base A
    echo
done
} 2>&1 | tee $RES/gt_check_dorado.txt
echo "=== gt check finished $(date) ==="
