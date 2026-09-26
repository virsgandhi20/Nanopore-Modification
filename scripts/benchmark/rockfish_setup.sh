#!/bin/bash
# Rockfish (lbcb-sci/rockfish, branch r10.4.1) for the RawMod paper's Table 1: CpG 5mC only.
# Login node only (needs internet for pip and the model download). Installs an isolated conda env,
# clones the R10.4.1 branch, downloads the 5 kHz model, prints the CLI, and runs a 200-read CPU
# smoke test on the E. coli M.SssI basecalls from the collection so the output format is known
# before the GPU driver is written. Re-runnable; every step is skipped when its output exists.
set -uo pipefail
ME=/fs/nexus-scratch/vgandhi
ENVD=$ME/envs/rockfish; SRC=$ME/rockfish_bench/rockfish; MODELS=$ME/rockfish_bench/models; T=$ME/rockfish_bench/smoke
CAT=/fs/cbcb-lab/storm/shared/data; DS=ecoli_k12_damdcm_mssi_r10.4.1_ontbasemod_2024
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
STATUS=$ME/rockfish_bench/setup_status.txt; mkdir -p $ME/rockfish_bench $T; : > $STATUS
note() { echo "$*" | tee -a $STATUS; }
source $HOME/miniconda3/etc/profile.d/conda.sh
[ -x $ENVD/bin/python ] || { conda create -y -q -p $ENVD python=3.10 pip > $ME/rockfish_bench/conda_create.log 2>&1 && note "env created" || { note "conda create FAILED (see conda_create.log)"; exit 1; }; }
conda activate $ENVD
[ -d $SRC/.git ] || { git clone -q --branch r10.4.1 https://github.com/lbcb-sci/rockfish.git $SRC && note "cloned r10.4.1 branch: $(cd $SRC && git log --oneline -1)" || { note "clone FAILED"; exit 1; }; }
# mappy (minimap2 binding) has no wheel and the login node has no compiler: take the prebuilt one from bioconda first
python -c "import mappy" 2>/dev/null || { conda install -y -q -c conda-forge -c bioconda mappy > $ME/rockfish_bench/conda_mappy.log 2>&1 && note "mappy installed from bioconda" || { note "mappy install FAILED (see conda_mappy.log)"; exit 1; }; }
python -c "import rockfish" 2>/dev/null || { (cd $SRC && pip install -q . > $ME/rockfish_bench/pip.log 2>&1) && note "pip install OK: torch $(python -c 'import torch; print(torch.__version__, torch.cuda.is_available())' 2>&1 | tail -1)" || { note "pip install FAILED (see pip.log): $(tail -3 $ME/rockfish_bench/pip.log | tr '\n' ' ')"; exit 1; }; }
ls $MODELS/* > /dev/null 2>&1 || { mkdir -p $MODELS && rockfish download -m 5kHz -s $MODELS > $ME/rockfish_bench/download.log 2>&1 && note "model downloaded: $(ls $MODELS | tr '\n' ' ')" || { note "model download FAILED: $(tail -2 $ME/rockfish_bench/download.log | tr '\n' ' ')"; exit 1; }; }
note "--- rockfish inference --help"; rockfish inference --help 2>&1 | head -60 | tee -a $STATUS
# 200-read CPU smoke test: subset the collection's M.SssI basecalls, run, show the output's head
[ -s $T/sub.bam ] || { $SAM view -h $CAT/$DS/basecalled/reads.bam | awk '/^@/ {print; next} c<200 {print; c++}' | $SAM view -b -o $T/sub.bam - && $SAM index $T/sub.bam; }
MODEL=$(ls $MODELS/* | head -1)
rockfish inference -i $CAT/$DS/pod5_files --bam_path $T/sub.bam --model_path $MODEL -r -t 4 -b 512 -o $T/out 2>&1 | tail -5 | tee -a $STATUS
note "--- smoke test outputs:"; ls -la $T | tee -a $STATUS; for f in $(ls $T | grep -v bam); do note "== $f"; head -5 $T/$f | tee -a $STATUS; done
note "setup finished $(date)"
