#!/bin/bash
# Dorado column of Table 1 for the bacteria rows, from the Dorado modbams that
# ship with the benchmark (sup v5, one modbam per modification model): modkit
# pileup -> the same site-level scorer and the same ground truth as the UniMeth
# column, so the two columns are directly comparable. CPU only.
#
# The benchmark ships four modbams per sample: 4mC (4mC_5mC model), 5mC (all-context),
# 5mCG (CpG-context model) and 6mA. The CpG row is scored with both 5mC models.
#
# DATASETS entries: <sample>[+<negative-control sample>]:<modbam tag>:<mod code>:<gt preset | bed | motif:...>:<row label>
#
# A CpG-context model only reports CpG positions, and in the M.SssI sample every
# CpG is methylated, so that sample alone has no negatives. With +<control> the
# same CpG positions of the unmethylated control sample are appended as
# negatives (contig names prefixed NEG_, so they can never match the ground truth).
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
BENCH=/fs/cbcb-lab/storm/bds062/data/benchmark
OUT=${OUT:-/fs/nexus-scratch/vgandhi/dorado_bench}
UB=/fs/nexus-scratch/vgandhi/unimeth_bench
MODKIT=/fs/nexus-scratch/vgandhi/dist_modkit_v0.6.4_cd85862/modkit
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
T=${SLURM_CPUS_PER_TASK:-8}
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null
mkdir -p $OUT $REPO/benchmark_results/dorado
declare -A REF=( [Ecoli_WT_5kHz]=ecoli.fa.gz [Ecoli_DM_MSssI_5kHz]=ecoli.fa.gz [Anabaena_WT_5kHz]=anabaena_sp_PCC7120_ATCC27893.fa.gz
                 [Tdenticola_WT_5kHz]=treponema_denticola_ATCC35405.fa.gz [HPJ99_WT_5kHz]=hpylori_J99_ATCC700824.fa.gz [HP26695_WT_5kHz]=hpylori_26695.fa.gz )
DATASETS=${DATASETS:-"Ecoli_DM_MSssI_5kHz+Ecoli_DM_5kHz:5mCG:m:ecoli_msssi:Ecoli_MSssI_vs_DM/5mC_CpG_cpgModel Ecoli_DM_MSssI_5kHz+Ecoli_DM_5kHz:5mC:m:ecoli_msssi:Ecoli_MSssI_vs_DM/5mC_CpG_allContextModel Ecoli_DM_MSssI_5kHz:5mC:m:ecoli_msssi:Ecoli_MSssI/5mC_CpG_allContextModel_withinSample Ecoli_WT_5kHz:5mC:m:ecoli_dcm:Ecoli_WT/5mC_nonCpG Ecoli_WT_5kHz:6mA:a:ecoli_dam:Ecoli_WT/6mA Tdenticola_WT_5kHz:6mA:a:tdenticola:Tdenticola/6mA_preset Anabaena_WT_5kHz:6mA:a:anabaena:Anabaena/6mA HPJ99_WT_5kHz:6mA:a:motif:GTNNNNNNAC:8:both:HPJ99/6mA_preset HPJ99_WT_5kHz:4mC:21839:motif:TCNNNNNNNGC:1:both:HPJ99/4mC_preset"}
TABLE=$OUT/table1_dorado.tsv; STATUS=$OUT/status.txt
[ "${APPEND:-0}" = 1 ] || { rm -f $TABLE; : > $STATUS; }
for entry in $DATASETS; do
    S=$(echo $entry | cut -d: -f1); TAG=$(echo $entry | cut -d: -f2); CODE=$(echo $entry | cut -d: -f3)
    NEG=""; [[ "$S" == *+* ]] && { NEG=${S#*+}; S=${S%%+*}; }
    LABEL=${entry##*:}; GT=$(echo $entry | cut -d: -f4- ); GT=${GT%:*}
    echo "==================== $LABEL  ($S, modbam $TAG, code $CODE, gt $GT)  $(date +%H:%M:%S)"
    pile() {  # sample -> sets PBED, returns 1 on failure
        local X=$1; local SRC=$BENCH/bacteria/$X/modbam/${X}_sup_v5r3_${TAG}.bam       # exact: 5mC and 5mCG are different models
        PBED=$OUT/${X}_${TAG}.bed
        [ -s "$SRC" ] || { echo "$LABEL: SKIPPED, no modbam ${X}_sup_v5r3_${TAG}.bam" | tee -a $STATUS; return 1; }
        for OLD in $UB/gtcheck/${X}_${TAG}.bed /fs/nexus-scratch/vgandhi/hp_labels/${X}_${TAG}.bed; do [ ! -s $PBED ] && [ -s $OLD ] && ln -sf $OLD $PBED; done   # reuse pileups from earlier work
        [ -s $PBED ] && return 0
        ln -sf $SRC $OUT/${X}_${TAG}.bam                      # share is read-only: index through a symlink
        [ -s $OUT/${X}_${TAG}.bam.bai ] || $SAM index -@ $T $OUT/${X}_${TAG}.bam
        $MODKIT pileup $OUT/${X}_${TAG}.bam $PBED.tmp --threads $T 2> $OUT/${X}_${TAG}.modkit.log && mv $PBED.tmp $PBED && return 0
        echo "$LABEL: FAILED modkit on $X: $(tail -1 $OUT/${X}_${TAG}.modkit.log)" | tee -a $STATUS; rm -f $PBED.tmp; return 1
    }
    pile $S || continue; BED=$PBED
    if [ -n "$NEG" ]; then
        pile $NEG || continue
        CALLS=$OUT/${S}_vs_${NEG}_${TAG}.bed
        { cat $BED; awk 'BEGIN{OFS="\t"} {$1="NEG_"$1; print}' $PBED; } > $CALLS; BED=$CALLS
    fi
    REFGZ=$BENCH/references/${REF[$S]}
    if [ -f "$GT" ]; then GTBED=$GT; else
        GTBED=$OUT/gt/$(echo "$GT" | tr ':' '_')_$S/gt_modified.bed
        if [ ! -s $GTBED ]; then
            if [[ "$GT" == motif:* ]]; then IFS=: read -r _ MOTIF OFFS STR <<< "$GT"
                python $REPO/scripts/ground_truth/motif_gt.py --ref $REFGZ --motif $MOTIF --mod-base ${MOTIF:$OFFS:1} --mod-offset $OFFS --strand $STR --outdir $(dirname $GTBED) > /dev/null 2>&1
            else python $REPO/scripts/ground_truth/motif_gt.py --ref $REFGZ --preset $GT --outdir $(dirname $GTBED) > /dev/null 2>&1; fi
        fi
    fi
    [ -s "$GTBED" ] || { echo "$LABEL: FAILED, no ground truth at $GTBED" | tee -a $STATUS; continue; }
    if python $REPO/scripts/benchmark/score_sites.py --calls $BED --gt $GTBED --min-cov 10 --chrom-col 0 --pos-col 1 --cov-col 9 \
        --freq-col 10 --freq-scale 100 --code-col 3 --code $CODE --label "$LABEL" --out $TABLE > $OUT/score_$(echo $LABEL | tr '/' '_').log 2>&1
    then echo "$LABEL: OK  $(tail -1 $TABLE | cut -f2-6)" | tee -a $STATUS
    else echo "$LABEL: FAILED scoring: $(tail -1 $OUT/score_$(echo $LABEL | tr '/' '_').log)" | tee -a $STATUS; fi
done
cp $TABLE $STATUS $REPO/benchmark_results/dorado/ 2>/dev/null
echo; echo "=== Dorado column ==="; column -t -s$'\t' $TABLE 2>/dev/null; echo; cat $STATUS
