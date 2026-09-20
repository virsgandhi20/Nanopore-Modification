#!/bin/bash
# Why is UniMeth blind on our E. coli data (mean P(6mA) = 0.0009 at Dam GATC
# sites, where methylation is ~100%)? Controlled 2x2 plus a positive control:
#
#   A  dorado 1.4.0   auto normalization (>0.7.1)   6mA    <- the failing config
#   B  dorado 1.4.0   forced legacy (0.7.1)         6mA
#   C  dorado 0.9.2   auto                          6mA    <- the paper's basecaller
#   D  dorado 0.9.2   forced legacy (0.7.1)         6mA
#   E  dorado 1.4.0   auto                          5mC vs Dcm CCWGG (6mA-specific?)
#   F  dorado 0.9.2   auto                          5mC vs Dcm CCWGG
#
# Built to survive unattended: NO `set -e`; every config is isolated, so one
# failure cannot take the others down; every failure is recorded with its
# reason in a status file that is tracked by git. If the 0.9.2 basecall fails,
# A/B/E still run. `MODE=preflight` checks everything that does not need a GPU.
set -uo pipefail

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
M6A=unimeth_r10.4.1_5kHz_6mA.pt; M5C=unimeth_r10.4.1_5kHz_5mC.pt
BAM140=$BASE/bam/$S.moves.bam
BAM092=$OUT/$S.dorado092.bam
NREADS=${NREADS:-2500}; LIMIT=${LIMIT:-600}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
RES=$REPO/benchmark_results/unimeth
SUMMARY=$RES/diag_summary.tsv
STATUS=$RES/diag_status.txt

mkdir -p $OUT $BASE/ref $RES
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null

# ---------------------------------------------------------------- preflight
if [ "${MODE:-}" = "preflight" ]; then
    FAILS=0
    chk() { if eval "$2" >/dev/null 2>&1; then echo "  PASS  $1"; else echo "  FAIL  $1"; FAILS=$((FAILS+1)); fi; }
    echo "=== preflight (no GPU needed) ==="
    chk "pod5 readable"                     "[ -r $POD5 ]"
    chk "reference readable"                "[ -r $REFGZ ]"
    chk "dorado-1.4.0 BAM exists"           "[ -s $BAM140 ]"
    chk "  ...and carries move tables"      "[ \$($SAM view $BAM140 | head -200 | grep -c 'mv:B') -gt 100 ]"
    chk "dorado 0.9.2 binary runs"          "$DOR092 --version"
    chk "  ...and supports --max-reads"     "$DOR092 basecaller --help 2>&1 | grep -q -- '--max-reads'"
    chk "  ...and supports --emit-moves"    "$DOR092 basecaller --help 2>&1 | grep -q -- '--emit-moves'"
    chk "dorado sup@v5.0.0 model present"   "[ -d $DMODEL ]"
    chk "UniMeth 6mA checkpoint"            "[ -s $CKPT/$M6A ]"
    chk "UniMeth 5mC checkpoint"            "[ -s $CKPT/$M5C ]"
    chk "unimeth-infer imports and runs"    "unimeth-infer --help"
    chk "python has numpy + scikit-learn"   "python -c 'import numpy, sklearn'"
    chk "can write to scratch output dir"   "touch $OUT/.w && rm $OUT/.w"
    chk "can write to repo results dir"     "touch $RES/.w && rm $RES/.w"
    for P in ecoli_dam ecoli_dcm; do
        [ -s $OUT/gt_$P/gt_modified.bed ] || python $REPO/scripts/ground_truth/motif_gt.py --ref $REFGZ --preset $P --outdir $OUT/gt_$P >/dev/null 2>&1
        chk "ground truth $P generated"     "[ \$(wc -l < $OUT/gt_$P/gt_modified.bed) -gt 1000 ]"
    done
    if [ -s $BASE/smoke/calls.txt ]; then
        chk "summary tool parses REAL UniMeth output" \
            "python $REPO/scripts/benchmark/diag_summary.py --calls $BASE/smoke/calls.txt --gt $OUT/gt_ecoli_dam/gt_modified.bed --label preflight --out $OUT/.pf.tsv"
        [ -s $OUT/.pf.tsv ] && { echo "        (expect the known-blind row:)"; tail -1 $OUT/.pf.tsv | sed 's/^/        /'; rm -f $OUT/.pf.tsv; }
    fi
    echo; [ $FAILS -eq 0 ] && echo "PREFLIGHT: ALL PASS - safe to submit" || echo "PREFLIGHT: $FAILS FAILURE(S) - do not submit, paste this output"
    exit $FAILS
fi

# ---------------------------------------------------------------- real run
rm -f $SUMMARY; : > $STATUS
note() { echo "$*" | tee -a $STATUS; }
note "diag started $(date)  host=$(hostname)  gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"

[ -s $REFFA ] || zcat $REFGZ > $REFFA
for P in ecoli_dam ecoli_dcm; do
    [ -s $OUT/gt_$P/gt_modified.bed ] || python $REPO/scripts/ground_truth/motif_gt.py --ref $REFGZ --preset $P --outdir $OUT/gt_$P >/dev/null 2>&1
    [ -s $OUT/gt_$P/gt_modified.bed ] && note "GT $P: $(wc -l < $OUT/gt_$P/gt_modified.bed) sites" || note "GT $P: FAILED to generate"
done

# dorado 0.9.2 basecall of a read subset (failure here must not stop A/B/E)
if [ ! -s $BAM092 ]; then
    $DOR092 basecaller $DMODEL $POD5 --emit-moves --reference $REFFA --max-reads $NREADS > $OUT/tmp092.bam 2> $OUT/dorado092.err
    if [ $? -eq 0 ] && [ -s $OUT/tmp092.bam ]; then
        $SAM sort -@ 8 -o $BAM092 $OUT/tmp092.bam 2>/dev/null && $SAM index $BAM092 && rm -f $OUT/tmp092.bam
    fi
fi
if [ -s $BAM092 ]; then note "dorado 0.9.2 basecall: OK ($($SAM view -c $BAM092) records)"
else note "dorado 0.9.2 basecall: FAILED -> C, D, F will be skipped. reason: $(tail -2 $OUT/dorado092.err 2>/dev/null | tr '\n' ' ')"; fi

run() {  # label bam model ctxflags gtpreset [extra unimeth flags]
    local L=$1 B=$2 M=$3 CTX=$4 G=$5; shift 5
    echo "==================== $L  ($(date +%H:%M:%S)) ===================="
    if [ ! -s "$B" ]; then note "$L: SKIPPED (no BAM)"; return; fi
    if [ ! -s $OUT/$L.txt ]; then
        unimeth-infer --pod5 $POD5 --bam $B --model $CKPT/$M --pore_type R10.4.1 --frequency 5khz $CTX \
            --output_format tsv --out $OUT/$L.txt --num_workers 8 --limit $LIMIT "$@" > $OUT/$L.runlog 2>&1
        local rc=$?
        if [ $rc -ne 0 ] || [ ! -s $OUT/$L.txt ]; then
            # --num_workers 8 is untested on a GPU node; fall back to UniMeth's
            # auto worker count, which is the configuration known to run.
            note "$L: first attempt failed (rc=$rc), retrying with default workers"
            rm -f $OUT/$L.txt
            unimeth-infer --pod5 $POD5 --bam $B --model $CKPT/$M --pore_type R10.4.1 --frequency 5khz $CTX \
                --output_format tsv --out $OUT/$L.txt --limit $LIMIT "$@" > $OUT/$L.runlog 2>&1
            rc=$?
        fi
        if [ $rc -ne 0 ] || [ ! -s $OUT/$L.txt ]; then
            note "$L: FAILED rc=$rc :: $(grep -E 'Error|error' $OUT/$L.runlog | tail -1)"
            rm -f $OUT/$L.txt; return
        fi
    fi
    local ver; ver=$(grep -E "Dorado Ver" $OUT/$L.runlog 2>/dev/null | head -1 | awk -F: '{print $2}' | tr -d ' ')
    if python $REPO/scripts/benchmark/diag_summary.py --calls $OUT/$L.txt --gt $OUT/gt_$G/gt_modified.bed --label $L --out $SUMMARY > $OUT/$L.sum 2>&1
    then note "$L: OK  (unimeth saw dorado=${ver:-?})  $(tail -1 $OUT/$L.sum | cut -f4,5,8 | awk '{print "meanP_gt="$1" meanP_bg="$2" auroc="$3}')"
    else note "$L: summary FAILED :: $(tail -1 $OUT/$L.sum)"; fi
}

run A_d140_auto_6mA    $BAM140 $M6A "--m6A 1" ecoli_dam
run B_d140_legacy_6mA  $BAM140 $M6A "--m6A 1" ecoli_dam --dorado_version 0.7.1
run C_d092_auto_6mA    $BAM092 $M6A "--m6A 1" ecoli_dam
run D_d092_legacy_6mA  $BAM092 $M6A "--m6A 1" ecoli_dam --dorado_version 0.7.1
run E_d140_auto_5mC    $BAM140 $M5C "--cpg 1 --chg 1 --chh 1" ecoli_dcm
run F_d092_auto_5mC    $BAM092 $M5C "--cpg 1 --chg 1 --chh 1" ecoli_dcm

note "diag finished $(date)"
echo; echo "=== SUMMARY ==="; column -t $SUMMARY 2>/dev/null; echo; echo "=== STATUS ==="; cat $STATUS
