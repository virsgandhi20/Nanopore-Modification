#!/bin/bash
# Round 2 launcher (login node). Step 1 needs no GPU and answers in a minute:
# what input does UniMeth's network see on our data vs on its own demo data?
# Step 2 queues the GPU configs G H P Q without wiping the A-F rows.
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
BASE=/fs/nexus-scratch/vgandhi/unimeth_bench
RES=$REPO/benchmark_results/unimeth
S=Ecoli_WT_5kHz
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null

rm -f $RES/sigstats.tsv
{
python $REPO/scripts/benchmark/unimeth_sigstats.py --label ours_ecoli --out $RES/sigstats.tsv \
    --bam $BASE/bam/$S.moves.bam --pod5 /fs/cbcb-lab/storm/bds062/data/benchmark/bacteria/$S/pod5/$S.pod5
if [ -s $BASE/demo/demo.bam ] && [ -s $BASE/demo/subset_18.pod5 ]; then
    python $REPO/scripts/benchmark/unimeth_sigstats.py --label unimeth_demo --out $RES/sigstats.tsv \
        --bam $BASE/demo/demo.bam --pod5 $BASE/demo/subset_18.pod5
else
    echo "[unimeth_demo] NOT DOWNLOADED: $BASE/demo is missing demo.bam / subset_18.pod5 (P and Q will be skipped)"
fi
} 2>&1 | tee $RES/sigstats.txt

sbatch --account=cbcb --partition=cbcb --qos=high --gres=gpu:1 --exclude=cbcb25 --job-name=um_diag2 \
    --mem=64G --cpus-per-task=12 --time=02:00:00 --output=/fs/nexus-scratch/vgandhi/um_diag2_%j.log \
    --export=ALL,APPEND=1,CONFIGS="G H P Q" --wrap="bash $REPO/scripts/benchmark/unimeth_diag.sh"
