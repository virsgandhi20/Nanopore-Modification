#!/bin/bash
# ORCA exact featurization for the Kulkarni et al. ONT basemod benchmark
# bacteria samples (/fs/cbcb-lab/storm/bds062/data/benchmark/bacteria).
#
# Layout per sample: pod5/$S.pod5 (single file) + modbam/${S}_sup_v5r3_*.bam
# (aligned, coordinate-sorted). The modbams all contain the same reads with
# different mod tags, so any one of them serves as the eventalign alignment;
# we use the 6mA one. References ship gzipped, f5c needs them uncompressed.
#
# Same single-contig cautions as SPO1: subsample (deep coverage on a small
# genome OOMs ORCA's per-contig concat) and the pileup sentinel (ORCA's
# index_pileup drops the last contig, which is the whole file here).
set -euo pipefail

F5C=${F5C:-/fs/nexus-scratch/vgandhi/f5c-v1.6/f5c_x86_64_linux}
SAM=${SAM:-/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools}
BENCH=${BENCH:-/fs/cbcb-lab/storm/bds062/data/benchmark}
OUT=${OUT:-/fs/nexus-scratch/vgandhi/orca_feat_bench}

SAMPLES=${SAMPLES:-"HP26695_WT_5kHz HP26695_WGA_5kHz"}
REFGZ=${REFGZ:-$BENCH/references/hpylori_26695.fa.gz}
ALNTAG=${ALNTAG:-6mA}           # which modbam to use as the alignment
FRAC=${FRAC:-0.2}               # subsample fraction; empty disables
MINRECALIB=${MINRECALIB:-200}
PORE=${PORE:-r10}
T=${T:-4}

mkdir -p $OUT
# uncompress the reference once
REF=$OUT/$(basename ${REFGZ%.gz})
[ -s $REF ] || zcat $REFGZ > $REF
[ -s $REF.fai ] || $SAM faidx $REF

for S in $SAMPLES; do
    W=$OUT/$S${FRAC:+_sub$FRAC}
    mkdir -p $W
    echo "==================== $S ===================="
    date

    POD5=$BENCH/bacteria/$S/pod5/$S.pod5
    SRCBAM=$BENCH/bacteria/$S/modbam/${S}_sup_v5r3_${ALNTAG}.bam
    for f in $POD5 $SRCBAM; do
        [ -f "$f" ] || { echo "MISSING: $f" >&2; exit 1; }
    done

    [ -f $W/$S.blow5 ] || blue-crab p2s $POD5 -o $W/$S.blow5

    # the share is read-only, so subsample (or just index) into the workspace
    if [ -n "$FRAC" ]; then
        if [ ! -s $W/$S.bam ]; then
            $SAM view -s $FRAC -b -@ $T $SRCBAM > $W/$S.bam
            $SAM index $W/$S.bam
        fi
    else
        [ -s $W/$S.bam ] || { cp $SRCBAM $W/$S.bam; $SAM index $W/$S.bam; }
    fi
    BAM=$W/$S.bam
    echo "$S: alignment reads = $($SAM view -c -@ $T $BAM)"

    [ -f $W/$S.fastq ] || $SAM fastq -F 0x900 $BAM > $W/$S.fastq

    if [ ! -s $W/$S.eventalign ]; then
        $F5C index --slow5 $W/$S.blow5 $W/$S.fastq
        $F5C eventalign --pore $PORE --min-recalib-events $MINRECALIB \
            --signal-index --scale-events --collapse-events --secondary=no -t $T \
            --slow5 $W/$S.blow5 --reads $W/$S.fastq \
            --bam $BAM --genome $REF \
            --summary $W/$S.summary > $W/$S.eventalign 2> $W/eventalign.log
    else
        echo "$S: reusing existing eventalign ($(du -h $W/$S.eventalign | cut -f1))"
    fi

    TOT=$(grep -c "^@" $W/$S.fastq || true)
    BAD=$(grep -c "could not calibrate" $W/eventalign.log || true)
    echo "$S reads=$TOT could-not-calibrate=$BAD"

    [ -s $W/$S.pileup ] || $SAM mpileup -f $REF $BAM > $W/$S.pileup 2>/dev/null

    # ORCA index_pileup last-contig workaround (see featurize_orca_spo1.sh)
    grep -q '^ZZZ_SENTINEL' $W/$S.pileup || printf 'ZZZ_SENTINEL\t1\tN\t0\t*\t*\n' >> $W/$S.pileup

    orca-pred_signal_feature_ext --eventalign $W/$S.eventalign --work_dir $W --prefix $S --n_processes $T
    orca-pred_bascal_feature_ext  --pileup     $W/$S.pileup     --work_dir $W --prefix $S --n_processes $T
    orca-pred_feature_merge       --work_dir   $W --prefix $S --n_processes $T

    echo "$S DONE: merged sites = $(wc -l < $W/$S.merged.feature.per.site)"
done
echo "=== BENCHMARK FEATURIZATION COMPLETE: $(date) ==="
