#!/bin/bash
# Login-node launcher for the two strand checks, plus an access probe for the
# pieces needed to time RawMod featurization (Table S3).
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
BENCH=/fs/cbcb-lab/storm/bds062/data/benchmark
RES=$REPO/benchmark_results/strand; mkdir -p $RES
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null
{
echo "################ access probe"
for P in /fs/nexus-scratch/bds062/rawhash2-env/rawhash2-storm/test/benchmark/scripts/7_refine_moves_remora.sh \
         /fs/nexus-scratch/bds062/rawhash2-env/rawhash2-storm/test/benchmark/scripts/3_run_dorado.sh \
         /fs/nexus-scratch/bds062/envs/mod/bin/python \
         /fs/cbcb-scratch/bds062/results/benchmark_results \
         /fs/cbcb-scratch/bds062/results/rawmod_full_pipeline4/features \
         /fs/cbcb-scratch/bds062/data/gt; do
    if [ -r "$P" ]; then echo "  readable   $P"; else echo "  NO ACCESS  $P"; fi
done
python -c "import h5py" 2>/dev/null && echo "  h5py: present in unimeth env" || echo "  h5py: missing in unimeth env"
python -c "import remora" 2>/dev/null && echo "  remora: present" || echo "  remora: missing"

echo; echo "################ k-mer level table"
LT=""
for C in /fs/nexus-scratch/bds062/rawhash2-env/rawhash2-storm/extern/local_kmer_models/uncalled_r1041_model_only_means.txt \
         $(find /fs/cbcb-lab/storm/shared/rawhash2 /fs/nexus-scratch/vgandhi $HOME -maxdepth 7 -name "uncalled_r1041_model_only_means.txt" 2>/dev/null | head -3); do
    [ -r "$C" ] && { LT=$C; break; }
done
echo "  using: ${LT:-NONE FOUND}"

if [ -n "$LT" ]; then
    echo; echo "################ 1. which sequence does a reverse read's signal follow?"
    python $REPO/scripts/test/strand_orientation_check.py --level-table $LT --ref /fs/nexus-scratch/vgandhi/unimeth_bench/ref/ecoli.fa \
        --bam /fs/nexus-scratch/vgandhi/unimeth_bench/bam/Ecoli_WT_5kHz.moves.bam \
        --pod5 $BENCH/bacteria/Ecoli_WT_5kHz/pod5/Ecoli_WT_5kHz.pod5
fi

echo; echo "################ 2. what forward-only featurization sees at GT sites"
for X in ecoli_dam:ecoli.fa.gz ecoli_dcm:ecoli.fa.gz anabaena:anabaena_sp_PCC7120_ATCC27893.fa.gz hpylori_26695:hpylori_26695.fa.gz \
         hpylori_j99:hpylori_J99_ATCC700824.fa.gz tdenticola:treponema_denticola_ATCC35405.fa.gz; do
    python $REPO/scripts/ground_truth/strand_label_audit.py --preset ${X%%:*} --ref $BENCH/references/${X##*:}
done
} 2>&1 | tee $RES/strand_check.txt
