#!/bin/bash
# Re-basecall SPO1 native barcodes WITH modification models, aligned to the
# SPO1 reference, so that (a) modkit can measure per-site m/h/a fractions and
# (b) per-read modification calls exist as weak labels for per-read typing.
#
# The BAMs on the lab share carry no MM/ML tags (basecalled without mod
# models), which is why modkit pileup found nothing to sample.
#
# Model versions match the RawMod paper (sup@v5.0.0). 4mC_5mC is not included
# because dorado cannot run two models that both call C; 5mC_5hmC + 6mA gives
# the m/h/a codes that mod_types.py expects.
#
# Step 1 (login node, needs internet):  MODE=download bash basecall_spo1_mods.sh
# Step 2 (GPU job):                      BARCODES="barcode06 barcode07" bash basecall_spo1_mods.sh
set -euo pipefail

DOR=${DOR:-/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-1.4.0-linux-x64/bin/dorado}
SAM=${SAM:-/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools}
MODKIT=${MODKIT:-/fs/nexus-scratch/vgandhi/dist_modkit_v0.6.4_cd85862/modkit}
DATA=${DATA:-/fs/cbcb-lab/storm/shared/umbc-ont-data}
REF=${REF:-$DATA/ref/SPO1_FJ230960.1.fasta}
MODELS=${MODELS:-/fs/nexus-scratch/vgandhi/dorado_models}
OUT=${OUT:-/fs/nexus-scratch/vgandhi/spo1_typing}
RUN=${RUN:-run1_jan31}
MODE_DIR=${MODE_DIR:-single_end}
BARCODES=${BARCODES:-"barcode06 barcode07"}
T=${T:-8}

BASE=dna_r10.4.1_e8.2_400bps_sup@v5.0.0
MODS="${BASE}_5mC_5hmC@v3 ${BASE}_6mA@v3"

mkdir -p $MODELS $OUT

if [ "${MODE:-}" = "download" ]; then
    for m in $BASE $MODS; do
        [ -d $MODELS/$m ] && { echo "have $m"; continue; }
        $DOR download --model $m --models-directory $MODELS
    done
    ls $MODELS; exit 0
fi

for m in $BASE $MODS; do
    [ -d $MODELS/$m ] || { echo "MISSING model $m: run with MODE=download on the login node first" >&2; exit 1; }
done
MODARGS=$(echo $MODS | sed "s#\([^ ]*\)#$MODELS/\1#g; s/ /,/g")

for BC in $BARCODES; do
    POD5=$DATA/pod5_by_barcode/$RUN/$MODE_DIR/$BC.pod5
    [ -f $POD5 ] || { echo "MISSING $POD5" >&2; exit 1; }
    echo "==================== $BC ===================="; date
    if [ ! -s $OUT/$BC.mod.sorted.bam ]; then
        $DOR basecaller $MODELS/$BASE $POD5 \
            --modified-bases-models $MODARGS \
            --reference $REF --emit-moves \
            > $OUT/$BC.mod.bam
        $SAM sort -@ $T -o $OUT/$BC.mod.sorted.bam $OUT/$BC.mod.bam
        $SAM index $OUT/$BC.mod.sorted.bam
        rm -f $OUT/$BC.mod.bam
    fi
    echo "$BC: MM-tagged reads in first 2000: $($SAM view $OUT/$BC.mod.sorted.bam | head -2000 | grep -c 'MM:Z:')"
    $MODKIT pileup $OUT/$BC.mod.sorted.bam $OUT/$BC.bed --threads $T 2> $OUT/$BC.modkit.log
    echo "$BC pileup rows: $(wc -l < $OUT/$BC.bed)"
done
echo "=== DONE $(date) ==="
