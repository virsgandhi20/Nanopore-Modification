#!/bin/bash
# Spot-check of UniMeth's fixed 5 kHz path (v0.3.2+, their issue 21, fixed 24 Sep 2026) against the matrix numbers,
# which were produced with v0.3.1 and the --frequency 4khz workaround. Two rows, one per published model:
#   E. coli M.SssI 5mC (ecoli_mssi, unimeth_5mC)  and  E. coli Dam 6mA (ecoli_wt, unimeth_6mA)
# Same sub.bam, same pod5, same ground truth and candidates, same scorer; only the UniMeth version and the flag differ.
#
#   MODE=setup  bash unimeth_v033_check.sh      # login node: clone tag 0.3.3, build env $ME/envs/unimeth033 (~10 min)
#   MODE=run    bash unimeth_v033_check.sh      # two GPU jobs; each ends by scoring and appending to $OUT/result.txt
#   MODE=status bash unimeth_v033_check.sh      # result.txt next to the matrix values
set -uo pipefail
ME=/fs/nexus-scratch/vgandhi; W=/fs/cbcb-lab/storm/vgandhi/matrix; CAT=/fs/cbcb-lab/storm/shared/data
REPO=$(cd "$(dirname "$0")/../.." && pwd); HERE=$REPO/scripts/benchmark
SRC=$ME/Unimeth_0.3.3; ENV=$ME/envs/unimeth033; OUT=$ME/unimeth_v033_check; mkdir -p $OUT/logs
MODELS=$ME/unimeth_models/checkpoints
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
export PIP_CACHE_DIR=$ME/.cache/pip CONDA_PKGS_DIRS=$ME/.cache/conda_pkgs TMPDIR=$ME/tmp; mkdir -p $PIP_CACHE_DIR $CONDA_PKGS_DIRS $TMPDIR
SB="--account=scavenger --partition=scavenger --qos=scavenger --requeue --gres=gpu:rtxa5000:1"
MODE=${MODE:-status}
case $MODE in
setup)
    source $HOME/miniconda3/etc/profile.d/conda.sh
    [ -d $SRC ] || git clone --depth 1 --branch 0.3.3 https://github.com/sekeyWang/Unimeth.git $SRC || exit 1
    [ -d $ENV ] || conda create -y -p $ENV python=3.12 pip || exit 1
    conda activate $ENV
    python -m pip install --only-binary=:all: "pod5==0.3.44" "lib-pod5==0.3.44" "pyarrow>=22,<23" || exit 1   # README: pod5/pyarrow pins
    (cd $SRC && python -m pip install .) || exit 1
    unimeth infer --help > $OUT/infer_help.txt 2>&1 && echo "setup OK: $(unimeth --version 2>/dev/null || grep -m1 -i version $OUT/infer_help.txt)"; echo "help saved to $OUT/infer_help.txt" ;;
run)
    # ds  tool  row  model  flags  types
    while read ds tool row model flags types; do
        T=$OUT/$ds.$tool; mkdir -p $T; rm -f $T/FAILED
        J=$(sbatch --parsable $SB --cpus-per-task=8 --mem=48G --time=04:00:00 --job-name=um033_$ds --output=$OUT/logs/${ds}_%j.log \
            --wrap="set -uo pipefail; source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $ENV; set -x
                    POD5=\$(awk -F'\t' '\$1==\"$ds\"{print \$3}' $HERE/matrix/datasets.tsv | sed 's#\${CAT}#$CAT#')
                    rm -f $T/calls.txt; unimeth infer --pod5 \$POD5 --bam $W/$ds/sub.bam --model $MODELS/$model --out $T/calls.txt --output_format tsv ${flags//_/ } --batch_size 256 --pore_type R10.4.1 --frequency 5khz > $T/infer.log 2>&1 || { echo \"infer failed: \$(grep -iE error $T/infer.log | tail -1)\" > $T/FAILED; exit 1; }
                    python $HERE/matrix/sites_std.py --tool unimeth --types '$types' $T/calls.txt $T/sites.std.tsv > $T/std.log 2>&1 || { echo \"sites_std failed: \$(tail -1 $T/std.log)\" > $T/FAILED; exit 1; }
                    python $HERE/score_sites.py --calls $T/sites.std.tsv --gt $W/rows/$row/gt.bed --candidates $W/rows/$row/cand_cov10.bed --min-cov 10 --chrom-col 0 --pos-col 1 --cov-col 2 --freq-col 4 --fill-missing --label $row/$tool/v0.3.3_5khz > $T/score.txt 2> $T/score.err
                    old=\$(awk -F'\t' -v r=$row -v t=$tool '\$1==r && \$2==t && \$3==\"mean_P\"{print \$8}' $W/matrix_long.tsv)
                    echo -e \"\$(tail -1 $T/score.txt)\tmatrix_v0.3.1_4khz_meanP=\$old\t\$(grep fill-missing $T/score.err)\" >> $OUT/result.txt")
        echo "$ds / $tool: job $J"
    done <<< "ecoli_mssi unimeth_5mC ecoli_mssi_5mC unimeth_r10.4.1_5kHz_5mC.pt --5mCpG_1_--5mCHG_1_--5mCHH_1 [CpG],[CHG],[CHH]
ecoli_wt unimeth_6mA ecoli_dam_6mA unimeth_r10.4.1_5kHz_6mA.pt --6mA_1 [m6A]" ;;
status)
    echo "== jobs"; squeue -h -u $USER -o "%.9i %.16j %.3t %.10M" | grep um033 || echo none
    echo "== results (label  n_sites  n_pos  pos_rate  AUROC  AUPRC  F1)"; cat $OUT/result.txt 2>/dev/null || echo "none yet"
    for d in $OUT/*/FAILED; do [ -f $d ] && echo "FAILED $(dirname $d): $(cat $d)"; done; true ;;
*) echo "MODE must be setup | run | status"; exit 1 ;;
esac
