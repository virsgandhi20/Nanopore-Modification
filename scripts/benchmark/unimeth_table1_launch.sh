#!/bin/bash
# Launch the UniMeth bacteria rows of Table 1 as parallel GPU jobs (login node).
# Clears everything produced under the old, blind normalization first, so no
# stale call file or table row can leak into the new table.
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
OUT=/fs/nexus-scratch/vgandhi/unimeth_bench
HPGT=/fs/nexus-scratch/vgandhi/hp_labels
declare -A LAST

for W in $OUT/*.6mA $OUT/*.5mC; do
    [ -d "$W" ] || continue
    if [ "$(cat $W/.norm 2>/dev/null)" != pa ]; then echo "clearing stale calls in $W"; rm -f $W/calls.txt $W/sites.tsv; fi
done
for T in table1_unimeth.tsv table1_sanity.tsv; do
    [ -s $OUT/$T ] && mv $OUT/$T $OUT/${T%.tsv}.blind_$(date +%m%d_%H%M).tsv
done
rm -rf $OUT/HPJ99_WT_5kHz.6mA/gt    # GT spec for this row changed, regenerate

# J99: GTNNNNNNAC is palindromic, the methylated A is at offset 8 on each strand
# (the hpylori_j99 preset's offset 1 is a T, i.e. only the minus-strand A).
for D in "Ecoli_WT_5kHz:6mA:ecoli_dam" "Ecoli_WT_5kHz:5mC:ecoli_dcm" "Ecoli_DM_MSssI_5kHz:5mC:ecoli_msssi" \
         "Anabaena_WT_5kHz:6mA:anabaena" "Tdenticola_WT_5kHz:6mA:tdenticola" "HPJ99_WT_5kHz:6mA:motif:GTNNNNNNAC:8:both" \
         "HP26695_WT_5kHz:6mA:$HPGT/gt_6mA.bed" "HP26695_WT_5kHz:5mC:$HPGT/gt_5mC.bed"; do
    S=${D%%:*}; N=$(echo $D | cut -d: -f1,2 | tr ':' '_')
    # two rows of one sample share a BAM: if it still has to be basecalled, the
    # second job must wait for the first or both would write the same file
    DEP=""; [ -n "${LAST[$S]:-}" ] && [ ! -s $OUT/bam/$S.moves.bam ] && DEP="--dependency=afterany:${LAST[$S]}"
    J=$(sbatch --parsable $DEP --account=cbcb --partition=cbcb --qos=high --gres=gpu:1 --exclude=cbcb25 --job-name=um_$N --mem=64G \
        --cpus-per-task=12 --time=06:00:00 --output=/fs/nexus-scratch/vgandhi/um_${N}_%j.log \
        --wrap="DATASETS='$D' LIMIT=${LIMIT:-15000} NUM_WORKERS=8 bash $REPO/scripts/benchmark/run_unimeth.sh")
    echo "submitted $J  $D  $DEP"; LAST[$S]=$J
done
squeue -u $USER
