#!/bin/bash
# DeepMod2 column, CpG 5mC (the only context it calls). First contact with this
# tool, so: MODE=install (login node) sets up an isolated env and prints the
# CLI; the default mode runs one row and records every step's outcome.
#
# Row: E. coli M.SssI (every CpG methylated = positives) vs the dam-/dcm- strain
# without M.SssI (same CpGs unmethylated = negatives). DeepMod2's
# bilstm_r10.4.1_5khz_v5.0 model matches our basecalling exactly
# (dna_r10.4.1_e8.2_400bps_sup@v5.0.0, move tables, aligned).
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
ME=/fs/nexus-scratch/vgandhi
OUT=${OUT:-$ME/deepmod2_bench}; SRC=$OUT/DeepMod2; ENVD=$ME/envs/deepmod2
BENCH=/fs/cbcb-lab/storm/bds062/data/benchmark
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
MODEL=${MODEL:-bilstm_r10.4.1_5khz_v5.0}; NREADS=${NREADS:-15000}
POS_BAM=${POS_BAM:-$ME/unimeth_bench/bam/Ecoli_DM_MSssI_5kHz.moves.bam}
NEG_BAM=${NEG_BAM:-$ME/typing_overnight/bam/Ecoli_DM_5kHz.moves.bam}
REFFA=$ME/unimeth_bench/ref/ecoli.fa
mkdir -p $OUT; STATUS=$OUT/status.txt
note() { echo "$*" | tee -a $STATUS; }

if [ "${MODE:-run}" = install ]; then
    source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth      # only to borrow a python >= 3.10
    export PIP_CACHE_DIR=$ME/pip_cache
    [ -d $SRC/.git ] || git clone -q --depth 1 https://github.com/WGLab/DeepMod2.git $SRC
    [ -x $ENVD/bin/python ] || python -m venv $ENVD
    $ENVD/bin/pip -q install --upgrade pip
    $ENVD/bin/pip -q install torch==2.5.1 --index-url https://download.pytorch.org/whl/cu124
    $ENVD/bin/pip -q install numpy numba pysam h5py tqdm ont-fast5-api pod5
    echo "--- versions"; $ENVD/bin/python -c "import torch, numpy, numba, pysam, pod5, h5py; print('torch', torch.__version__, '| numpy', numpy.__version__, '| numba', numba.__version__, '| pod5', pod5.__version__)"
    echo "--- models shipped in the repository"; ls $SRC/models 2>/dev/null | head -20; find $SRC -maxdepth 3 -iname "*5khz_v5*" | head
    echo "--- CLI"; $ENVD/bin/python $SRC/deepmod2 detect --help 2>&1 | head -60
    du -sh $ENVD $SRC; df -h $ME | tail -1
    exit 0
fi

: > $STATUS; note "deepmod2 run started $(date) host=$(hostname) gpu=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1)"
PY=$ENVD/bin/python
for tag in pos neg; do
    [ $tag = pos ] && { BAM=$POS_BAM; S=Ecoli_DM_MSssI_5kHz; } || { BAM=$NEG_BAM; S=Ecoli_DM_5kHz; }
    W=$OUT/$tag; mkdir -p $W
    if [ ! -s "$BAM" ]; then note "$tag ($S): SKIPPED, no move-table BAM at $BAM"; continue; fi
    # first N reads of the coordinate-sorted BAM = a contiguous region at full depth (as for UniMeth)
    [ -s $W/sub.bam ] || { $SAM view -h $BAM | awk -v n=$NREADS '/^@/ {print; next} c<n {print; c++}' | $SAM view -b -o $W/sub.bam - && $SAM index $W/sub.bam; }
    if ! ls $W/calls/*per_site* > /dev/null 2>&1; then
        $PY $SRC/deepmod2 detect --bam $W/sub.bam --input $BENCH/bacteria/$S/pod5 --file_type pod5 --model $MODEL \
            --ref $REFFA --threads 12 --output $W/calls > $W/detect.log 2>&1
        rc=$?; [ $rc -eq 0 ] || { note "$tag ($S): deepmod2 detect FAILED rc=$rc :: $(grep -iE 'error|exception' $W/detect.log | tail -1)"; continue; }
    fi
    PS=$(ls $W/calls/*per_site 2>/dev/null | head -1)
    [ -s "$PS" ] && note "$tag ($S): OK, $(wc -l < $PS) per-site rows; header: $(head -1 $PS | tr '\t' ' ')" || note "$tag ($S): no per_site output in $W/calls: $(ls $W/calls 2>/dev/null | tr '\n' ' ')"
done

# standardize by header name (chrom, 0-based pos, coverage, modified fraction); negatives get NEG_ contigs
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth
python - $OUT <<'PY' 2>&1 | tee -a $STATUS
import glob, sys
out = sys.argv[1]; rows = 0
with open(f"{out}/combined_sites.tsv", "w") as fo:
    for tag, prefix in (("pos", ""), ("neg", "NEG_")):
        fs = [f for f in glob.glob(f"{out}/{tag}/calls/*per_site") ]
        if not fs: print(f"standardize: no per_site file for {tag}"); continue
        with open(fs[0]) as f:
            hdr = f.readline().lstrip("#").rstrip("\n").split("\t"); col = {h.strip().lower(): i for i, h in enumerate(hdr)}
            need = {"chrom": ["chromosome", "chrom", "contig"], "pos": ["position_before", "start", "pos0"], "cov": ["coverage", "cov"],
                    "frac": ["mod_fraction", "mod_percentage", "fraction", "modified_fraction"]}
            idx = {k: next((col[c] for c in v if c in col), None) for k, v in need.items()}
            if None in idx.values(): print(f"standardize: unexpected {tag} header {hdr} -> {idx}; stopping so the columns can be set by hand"); sys.exit(0)
            scale = 100.0 if "percentage" in hdr[idx["frac"]].lower() else 1.0
            for line in f:
                c = line.rstrip("\n").split("\t")
                try: fo.write(f"{prefix}{c[idx['chrom']]}\t{int(c[idx['pos']])}\t{c[idx['cov']]}\t{float(c[idx['frac']])/scale}\n"); rows += 1
                except (ValueError, IndexError): pass
print(f"standardize: {rows} site rows written")
PY
GTBED=$OUT/gt/gt_modified.bed
[ -s $GTBED ] || python $REPO/scripts/ground_truth/motif_gt.py --ref $BENCH/references/ecoli.fa.gz --preset ecoli_msssi --outdir $OUT/gt > /dev/null 2>&1
if [ -s $OUT/combined_sites.tsv ] && [ -s $GTBED ]; then
    rm -f $OUT/table1_deepmod2.tsv
    python $REPO/scripts/benchmark/score_sites.py --calls $OUT/combined_sites.tsv --gt $GTBED --min-cov 10 --chrom-col 0 --pos-col 1 --cov-col 2 --freq-col 3 \
        --label "Ecoli_MSssI_vs_DM/5mC_CpG/DeepMod2_$MODEL" --out $OUT/table1_deepmod2.tsv > $OUT/score.log 2>&1 \
        && note "score: OK $(tail -1 $OUT/table1_deepmod2.tsv)" || note "score: FAILED $(tail -1 $OUT/score.log)"
fi
note "deepmod2 run finished $(date)"
