#!/bin/bash
# Score UniMeth's 6mA site calls against a ground truth derived from an
# independent caller (Dorado 6mA modbam -> modkit pileup, already computed by
# unimeth_gt_check.sh): positive = Dorado >= 70% modified, negative = <= 10%,
# both at coverage >= 10; sites in between are left out. For organisms whose
# motif preset turned out not to be methylated (T. denticola, H. pylori J99).
# E. coli is the control: it should reproduce the motif-GT number (0.9998).
# Login-node safe (a couple of minutes).
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
OUT=/fs/nexus-scratch/vgandhi/unimeth_bench
RES=$REPO/benchmark_results/unimeth
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null
T=$OUT/table1_unimeth_vs_dorado.tsv; rm -f $T
for S in Ecoli_WT_5kHz Tdenticola_WT_5kHz HPJ99_WT_5kHz; do
    B=$OUT/gtcheck/${S}_6mA.bed; G=$OUT/gtcheck/$S
    [ -s $B ] || { echo "$S: no Dorado pileup at $B"; continue; }
    awk '$4=="a" && $10>=10 && $11>=70 {print $1"\t"$2"\t"$3}' $B > $G.dorado_pos.bed
    awk '$4=="a" && $10>=10 && ($11>=70 || $11<=10) {print $1"\t"$2"\t"$3}' $B > $G.dorado_cand.bed
    echo "$S: Dorado-derived GT: $(wc -l < $G.dorado_pos.bed) positives among $(wc -l < $G.dorado_cand.bed) confident sites"
    for MODE in freq meanprob; do
        EXTRA=""; [ $MODE = meanprob ] && EXTRA="--num-col 5"
        python $REPO/scripts/benchmark/score_sites.py --calls $OUT/$S.6mA/sites.tsv --gt $G.dorado_pos.bed --candidates $G.dorado_cand.bed \
            --min-cov 10 --chrom-col 0 --pos-col 1 --cov-col 8 $EXTRA --label "$S/6mA/$MODE" --out $T > /dev/null 2>&1 || echo "  scoring failed: $S $MODE"
    done
done
echo; column -t $T; cp $T $RES/ 2>/dev/null
echo; echo "=== HP26695 motifs (the run that crashed on the reference header) ==="
python $REPO/scripts/benchmark/motif_enrichment.py --sites $OUT/HP26695_WT_5kHz.6mA/sites.tsv \
    --ref /fs/cbcb-lab/storm/bds062/data/benchmark/references/hpylori_26695.fa.gz \
    --gt /fs/nexus-scratch/vgandhi/hp_labels/gt_6mA.bed --label "HP26695_WT_5kHz  UniMeth 6mA" --base A 2>&1 | tee $RES/gt_check_hp26695.txt
