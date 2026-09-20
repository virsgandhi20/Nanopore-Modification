#!/bin/bash
# UniMeth benchmark for the RawMod paper's Table 1 (site-level AUROC).
#
# Per dataset: (1) get a Dorado BAM WITH move tables (UniMeth needs --emit-moves;
# re-basecalls from pod5 if the available BAM lacks mv tags), (2) unimeth infer
# -> per-read TSV, (3) call_modification_frequency.py -> per-site frequencies,
# (4) ground truth from motif_gt.py presets (bacteria) or a supplied BED,
# (5) score_sites.py -> one row in $OUT/table1_unimeth.tsv.
#
# UniMeth covers 5mC (CpG/CHG/CHH) and 6mA only; 4mC and 5hmU rows are N/A.
#
# Usage: DATASETS="Ecoli_WT_5kHz:6mA:ecoli_dam Anabaena_WT_5kHz:6mA:anabaena" bash run_unimeth.sh
#   each entry = <benchmark sample>:<mod>:<gt>, gt = preset | bed path | motif:<IUPAC>:<offset>:<+|both>
#   mod = 6mA | 5mC (all contexts) | CpG (human/mouse CpG-only model)
set -euo pipefail

BENCH=${BENCH:-/fs/cbcb-lab/storm/bds062/data/benchmark}
REFS=${REFS:-$BENCH/references}
OUT=${OUT:-/fs/nexus-scratch/vgandhi/unimeth_bench}
MODELS=${MODELS:-/fs/nexus-scratch/vgandhi/unimeth_models/checkpoints}
DORADO=${DORADO:-/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-1.4.0-linux-x64/bin/dorado}
DORADO_MODEL=${DORADO_MODEL:-/fs/nexus-scratch/vgandhi/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0}
SAM=${SAM:-/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools}
UNIMETH_ENV=${UNIMETH_ENV:-/fs/nexus-scratch/vgandhi/envs/unimeth}
UNIMETH_SRC=${UNIMETH_SRC:-/fs/nexus-scratch/vgandhi/Unimeth}   # git clone: editable install + scripts/
# NOTE: install UniMeth from this clone (pip install --no-deps -e), NOT from PyPI:
# the PyPI wheel is v0.1.0 and ships without its configs/ directory.
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
MIN_COV=${MIN_COV:-10}
# per-mod model files (names as shipped on the UniMeth Google Drive; override if different)
MODEL_5mC=${MODEL_5mC:-$MODELS/unimeth_r10.4.1_5kHz_5mC.pt}
MODEL_6mA=${MODEL_6mA:-$MODELS/unimeth_r10.4.1_5kHz_6mA.pt}
MODEL_CPG=${MODEL_CPG:-$MODELS/unimeth_r10.4.1_5kHz_CpG.pt}   # human/mouse CpG-only rows

# sample -> reference fasta.gz basename
declare -A REF=( [Ecoli_WT_5kHz]=ecoli.fa.gz [Ecoli_DM_5kHz]=ecoli.fa.gz [Ecoli_DM_MSssI_5kHz]=ecoli.fa.gz
                 [Anabaena_WT_5kHz]=anabaena_sp_PCC7120_ATCC27893.fa.gz
                 [Tdenticola_WT_5kHz]=treponema_denticola_ATCC35405.fa.gz
                 [HPJ99_WT_5kHz]=hpylori_J99_ATCC700824.fa.gz [HP26695_WT_5kHz]=hpylori_26695.fa.gz )

mkdir -p $OUT $OUT/ref
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh
conda activate $UNIMETH_ENV

for entry in ${DATASETS:?set DATASETS}; do
    IFS=: read -r S MOD GT <<< "$entry"
    W=$OUT/$S.$MOD; mkdir -p $W
    echo "==================== $S  mod=$MOD  gt=$GT ===================="; date

    REFGZ=$REFS/${REF[$S]}; REFFA=$OUT/ref/${REF[$S]%.gz}
    [ -s $REFFA ] || zcat $REFGZ > $REFFA
    [ -s $REFFA.fai ] || $SAM faidx $REFFA
    POD5=$BENCH/bacteria/$S/pod5/$S.pod5

    # (1) BAM with move tables, ONE per sample (shared by the 6mA and 5mC rows).
    # Prefer an existing mv-tagged BAM; else basecall. Adopt a BAM left by an
    # earlier run that stored it under the per-mod workdir.
    mkdir -p $OUT/bam; BAM=$OUT/bam/$S.moves.bam
    for old in $OUT/$S.*/$S.moves.bam; do
        [ -s "$old" ] && [ ! -s $BAM ] && { mv $old $BAM; mv $old.bai $BAM.bai 2>/dev/null || true; }
    done
    if [ ! -s $BAM ]; then
        SRC=$(ls $BENCH/bacteria/$S/modbam/*.bam 2>/dev/null | head -1)
        if [ -n "$SRC" ] && $SAM view $SRC | head -500 | grep -q "mv:B"; then
            echo "$S: existing modbam carries move tables, reusing"; ln -sf $SRC $BAM
            [ -s $BAM.bai ] || $SAM index $BAM
        else
            echo "$S: basecalling with --emit-moves (GPU)"
            $DORADO basecaller $DORADO_MODEL $POD5 --emit-moves --reference $REFFA > $OUT/bam/$S.unsorted.bam
            $SAM sort -@ 8 -o $BAM $OUT/bam/$S.unsorted.bam && $SAM index $BAM && rm $OUT/bam/$S.unsorted.bam
        fi
    fi

    [ -s $BAM.bai ] || $SAM index $BAM

    # (2) unimeth-infer -> per-read TSV
    case "$MOD" in
        6mA)  MODEL=$MODEL_6mA; CTX="--m6A 1" ;;
        CpG)  MODEL=$MODEL_CPG; CTX="--cpg 1" ;;                   # human/mouse CpG model
        *)    MODEL=$MODEL_5mC; CTX="--cpg 1 --chg 1 --chh 1" ;;   # all-context 5mC
    esac
    # v0.3.1 CLI (the README's `unimeth infer` is stale): separate `unimeth-infer`
    # command, m6A is an explicit switch, TSV path must end in .txt
    [ -s $W/calls.txt ] || unimeth-infer --pod5 $POD5 --bam $BAM --model $MODEL \
        --pore_type R10.4.1 --frequency 5khz $CTX --output_format tsv --out $W/calls.txt \
        --batch_size ${BATCH:-256} ${LIMIT:+--limit $LIMIT}

    # (3) per-site frequency
    [ -s $W/sites.tsv ] || python $UNIMETH_SRC/scripts/call_modification_frequency.py \
        -i $W/calls.txt -o $W/sites.tsv --sort

    # (4) ground truth: preset name -> generate from the reference; else a BED path
    # GT spec: a BED path | a motif_gt preset name | motif:<IUPAC>:<offset>:<+|both>
    # (custom motif lets a mixed preset like hpylori_j99 be restricted to one mod)
    if [ -f "$GT" ]; then GTBED=$GT; else
        GTBED=$W/gt/gt_modified.bed
        if [ ! -s $GTBED ]; then
            if [[ "$GT" == motif:* ]]; then
                IFS=: read -r _ MOTIF OFFS STR <<< "$GT"
                python $REPO/scripts/ground_truth/motif_gt.py --ref $REFGZ --motif $MOTIF \
                    --mod-base ${MOTIF:$OFFS:1} --mod-offset $OFFS --strand $STR --outdir $W/gt
            else
                python $REPO/scripts/ground_truth/motif_gt.py --ref $REFGZ --preset $GT --outdir $W/gt
            fi
        fi
    fi

    # (5) score: UniMeth freq output columns: chrom=0 pos=1 ... coverage=8 freq=9
    python $REPO/scripts/benchmark/score_sites.py --calls $W/sites.tsv --gt $GTBED \
        --min-cov $MIN_COV --chrom-col 0 --pos-col 1 --cov-col 8 --freq-col 9 \
        --label "$S/$MOD" --out $OUT/table1_unimeth.tsv
done
echo "=== DONE $(date) ==="; cat $OUT/table1_unimeth.tsv
