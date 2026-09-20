#!/bin/bash
# Why is UniMeth blind on our E. coli data (mean P(6mA) = 0.0009 at Dam GATC
# sites, where methylation is ~100%)? Controlled 2x2 plus a positive control:
#
#   basecaller   x  normalization branch          model
#   A  dorado 1.4.0   auto (new, >0.7.1)           6mA    <- the failing config
#   B  dorado 1.4.0   forced legacy (0.7.1)        6mA
#   C  dorado 0.9.2   auto (new)                   6mA    <- the paper's basecaller
#   D  dorado 0.9.2   forced legacy (0.7.1)        6mA
#   E  dorado 1.4.0   auto                         5mC vs Dcm CCWGG  (is it 6mA-specific?)
#   F  dorado 0.9.2   auto                         5mC vs Dcm CCWGG
#
# UniMeth's normalization boundary is 0.7.1, so 0.9.2 and 1.4.0 share a branch;
# if C works and A does not, dorado 1.x changed something UniMeth relies on.
# If only legacy works, its version switch is wrong for modern dorado. If
# E/F work but A-D do not, the 6mA model does not detect bacterial Dam 6mA.
set -euo pipefail

BENCH=${BENCH:-/fs/cbcb-lab/storm/bds062/data/benchmark}
BASE=${BASE:-/fs/nexus-scratch/vgandhi/unimeth_bench}
OUT=$BASE/diag
S=Ecoli_WT_5kHz
POD5=$BENCH/bacteria/$S/pod5/$S.pod5
REFGZ=$BENCH/references/ecoli.fa.gz
REFFA=$BASE/ref/ecoli.fa
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
DOR092=/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-0.9.2-linux-x64/bin/dorado
DMODEL=/fs/nexus-scratch/vgandhi/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
CKPT=/fs/nexus-scratch/vgandhi/unimeth_models/checkpoints
BAM140=$BASE/bam/$S.moves.bam
NREADS=${NREADS:-2500}      # reads to basecall with 0.9.2
LIMIT=${LIMIT:-600}         # reads per UniMeth config
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
SUMMARY=$REPO/benchmark_results/unimeth/diag_summary.tsv

mkdir -p $OUT $BASE/ref $(dirname $SUMMARY); rm -f $SUMMARY
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth
[ -s $REFFA ] || zcat $REFGZ > $REFFA

# ground truth: Dam (6mA at GATC) and Dcm (5mC at CCWGG)
for P in ecoli_dam ecoli_dcm; do
    [ -s $OUT/gt_$P/gt_modified.bed ] || python $REPO/scripts/ground_truth/motif_gt.py --ref $REFGZ --preset $P --outdir $OUT/gt_$P > /dev/null
done

# dorado 0.9.2 basecall of a small read subset, with move tables
BAM092=$OUT/$S.dorado092.bam
if [ ! -s $BAM092 ]; then
    $DOR092 basecaller $DMODEL $POD5 --emit-moves --reference $REFFA --max-reads $NREADS > $OUT/tmp092.bam
    $SAM sort -@ 8 -o $BAM092 $OUT/tmp092.bam && $SAM index $BAM092 && rm $OUT/tmp092.bam
fi
echo "dorado versions in BAM headers:"; for b in $BAM140 $BAM092; do $SAM view -H $b | grep "^@PG" | grep -o "VN:[^	]*" | head -1; done

run() {  # label bam model ctxflags gtpreset [extra unimeth flags]
    local L=$1 B=$2 M=$3 CTX=$4 G=$5; shift 5
    echo "==================== $L ===================="
    [ -s $OUT/$L.txt ] || unimeth-infer --pod5 $POD5 --bam $B --model $CKPT/$M --pore_type R10.4.1 \
        --frequency 5khz $CTX --output_format tsv --out $OUT/$L.txt --num_workers 8 --limit $LIMIT "$@" \
        2>&1 | grep -E "Dorado Ver|Dorado Source|Error|Traceback" || true
    python $REPO/scripts/benchmark/diag_summary.py --calls $OUT/$L.txt --gt $OUT/gt_$G/gt_modified.bed --label $L --out $SUMMARY | tail -1
}

run A_d140_auto_6mA    $BAM140 unimeth_r10.4.1_5kHz_6mA.pt "--m6A 1" ecoli_dam
run B_d140_legacy_6mA  $BAM140 unimeth_r10.4.1_5kHz_6mA.pt "--m6A 1" ecoli_dam --dorado_version 0.7.1
run C_d092_auto_6mA    $BAM092 unimeth_r10.4.1_5kHz_6mA.pt "--m6A 1" ecoli_dam
run D_d092_legacy_6mA  $BAM092 unimeth_r10.4.1_5kHz_6mA.pt "--m6A 1" ecoli_dam --dorado_version 0.7.1
run E_d140_auto_5mC    $BAM140 unimeth_r10.4.1_5kHz_5mC.pt "--cpg 1 --chg 1 --chh 1" ecoli_dcm
run F_d092_auto_5mC    $BAM092 unimeth_r10.4.1_5kHz_5mC.pt "--cpg 1 --chg 1 --chh 1" ecoli_dcm

echo; echo "=== SUMMARY (also written to $SUMMARY) ==="; column -t $SUMMARY
