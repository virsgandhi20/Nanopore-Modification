#!/bin/bash
# Reproducibility test of Bhargav's public RawMod repository (bds062/Nanopore-Modification) on the UMIACS cluster.
# Runs the smallest real benchmark pair (ONT oligos syn_5mC vs syn_control, ~3 GB) through the repo's own scripts.
#
#   bash cluster_test.sh check     # every path this script relies on exists?        (login node, seconds)
#   bash cluster_test.sh env       # clone from the bundle + two conda envs + selftest     (login node, ~10 min)
#   bash cluster_test.sh gpu       # rawmod-selftest on a GPU node                    (srun, ~1 min)
#   bash cluster_test.sh torch     # swap torch for the CUDA 12.6 build (all GPU generations) (login node, ~5 min)
#   bash cluster_test.sh config    # write test/benchmark/config/site.env, show it    (login node, seconds)
#   bash cluster_test.sh data      # download the oligo dataset with the repo script  (login node, ~3 GB)
#   bash cluster_test.sh run       # submit the benchmark for syn_5mC + syn_control   (sbatch chain)
#   bash cluster_test.sh status    # what has run, and the queue                      (any time)
#   bash cluster_test.sh score     # tables + summary once the jobs are done
#
# Every stage ends with "STAGE <name> OK"; paste the last ~20 lines of each stage back.
set -euo pipefail

ME=/fs/nexus-scratch/vgandhi
ROOT=/fs/cbcb-lab/storm/vgandhi/rawmod-test   # clone + env + caches (lab share: scratch is 97% full)
REPO_OLD=$ROOT/repo                            # public bds062 commit 5f6e00d; the Oct 28 job chain ran from here
REPO=$ROOT/RawMod                              # STORMgroup/RawMod e309935 ("use this instead", Oct 28 12:53): private, so it
BUNDLE=$HOME/RawMod.bundle                     #   arrives as a git bundle copied over with scp and is cloned from that
ENV=$ROOT/env                                  # RawMod env: exactly the README's `pip install -e .`
ENV_REFINE=$ROOT/env_refine                    # Remora 3.3.0 for the move-table refinement (README: separate env; it pins numpy<2)
WORK=/fs/cbcb-lab/storm/vgandhi/rawmod-test   # outputs + downloaded data (lab share, not scratch)
CONDA_SH=$HOME/miniconda3/etc/profile.d/conda.sh
DORADO=/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-1.4.0-linux-x64/bin/dorado
DORADO_MODELS=$ME/dorado_models                # holds dna_r10.4.1_e8.2_400bps_sup@v5.0.0
SAMTOOLS=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
EXCL=--exclude=legacy43                        # legacy43 has neither /fs/cbcb-lab nor /fs/nexus-scratch mounted (rows job died there, exit 53, Oct 28)
SLURM_CPU="--account=cbcb --partition=cbcb --qos=high --cpus-per-task=8 --mem=32G --time=06:00:00 $EXCL"
SLURM_GPU="--account=cbcb --partition=cbcb --qos=high --gres=gpu:1 --cpus-per-task=8 --mem=32G --time=06:00:00 $EXCL"

# keep conda/pip caches off the home quota
export CONDA_PKGS_DIRS=$ROOT/.conda_pkgs PIP_CACHE_DIR=$ROOT/.pip_cache
mkdir -p "$CONDA_PKGS_DIRS" "$PIP_CACHE_DIR" "$ROOT/tmp"
# TMPDIR on the lab share is set only inside the env stage (pip unpacks 3 GB of wheels there, home is too small).
# It must NOT reach the sbatch jobs (--export=ALL): Python's multiprocessing then puts its temp dir on NFS and its
# exit-time cleanup fails on .nfs lock files (seen in both tool jobs, Oct 28).
unset TMPDIR

activate() { source "$CONDA_SH"; conda activate "$ENV"; [ -d "$REPO/.git" ] || REPO=$REPO_OLD; cd "$REPO"; }   # the new clone exists only after env

case "${1:-}" in
  check)
    for p in "$BUNDLE" "$CONDA_SH" "$DORADO" "$DORADO_MODELS/dna_r10.4.1_e8.2_400bps_sup@v5.0.0" "$SAMTOOLS" /fs/cbcb-lab/storm/vgandhi "$ME"; do
      [ -e "$p" ] && echo "ok       $p" || echo "MISSING  $p"
    done
    echo "--- quotas"; df -h "$ME" /fs/cbcb-lab/storm 2>/dev/null | tail -2
    echo "--- sbatch: $(command -v sbatch || echo not on PATH)"
    echo "STAGE check OK";;

  env)
    export TMPDIR=$ROOT/tmp
    mkdir -p "$ROOT" "$WORK"
    if [ ! -d "$REPO/.git" ]; then
      [ -f "$BUNDLE" ] || { echo "MISSING $BUNDLE: scp it from the Mac first"; exit 1; }
      git clone --quiet -b main "$BUNDLE" "$REPO"
    fi
    cd "$REPO" && git log --oneline -1
    source "$CONDA_SH"
    # RawMod env = exactly the README line. A first attempt put Remora in here too; Remora 3.3.0 pins numpy<2
    # and pulled numpy 1.26 in, so an env whose numpy is not 2.x is rebuilt from scratch (wheels are cached).
    if ! "$ENV/bin/python" -c "import numpy; assert numpy.__version__.startswith('2')" 2>/dev/null; then
      [ -d "$ENV" ] && conda remove -y -p "$ENV" --all
      conda create -y -p "$ENV" python=3.11
    fi
    "$ENV/bin/pip" install -e .
    # refinement env, as the README says (separate); awscli lives here too, for the download script
    [ -x "$ENV_REFINE/bin/python" ] || conda create -y -p "$ENV_REFINE" python=3.11
    "$ENV_REFINE/bin/pip" install "ont-remora==3.3.0" awscli
    "$ENV/bin/python" -c "import numpy, pandas, pod5, pysam, torch; print('rawmod env: numpy', numpy.__version__, 'pandas', pandas.__version__, 'pod5', pod5.__version__, 'torch', torch.__version__)"
    "$ENV_REFINE/bin/python" -c "import remora, numpy, torch; print('refine env: remora', remora.__version__, 'numpy', numpy.__version__, 'torch', torch.__version__)"
    "$ENV/bin/rawmod-selftest"
    echo "STAGE env OK";;

  gpu)
    activate
    srun $SLURM_GPU --time=00:15:00 --mem=8G --cpus-per-task=2 bash -c "nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv; for p in $ENV $ENV_REFINE; do \$p/bin/python -c 'import torch; print(\"torch\", torch.__version__, \"sees cuda:\", torch.cuda.is_available())'; done; $ENV/bin/rawmod-selftest"
    echo "STAGE gpu OK (both 'sees cuda' lines must be True and the self-test line must say 'device cuda')";;

  torch)
    # The cbcb partition mixes GPU generations. `pip install -e .` pulls the CUDA 13 build of torch, which has no kernels
    # for the GTX 1080 Ti nodes (compute capability 6.1): tool:syn_control died there with "no kernel image is available"
    # (Oct 28) while tool:syn_5mC passed on an RTX 2080 Ti. A CUDA 12.6 build of torch covers both.
    activate
    export TMPDIR=$ROOT/tmp
    pip install "torch==2.14.1+cu126" --index-url https://download.pytorch.org/whl/cu126   # the cu126 index had no 2.15.0 on Oct 28
    python -c "import torch; print(torch.__version__, 'archs:', ' '.join(torch.cuda.get_arch_list()))"
    rawmod-selftest | tail -2
    echo "STAGE torch OK (archs must include sm_61 and sm_75)";;

  config)
    activate
    mkdir -p "$WORK/work" "$WORK/data"
    cat > test/benchmark/config/site.env <<EOF
BENCH_WORK=$WORK/work
RAWMOD_PYTHON=$ENV/bin/python
REFINE_PYTHON=$ENV_REFINE/bin/python
POD5=$ENV/bin/pod5
SAMTOOLS=$SAMTOOLS
DORADO=$DORADO
DORADO_MODELS=$DORADO_MODELS
BENCH_SLURM_CPU="$SLURM_CPU"
BENCH_SLURM_GPU="$SLURM_GPU"
EOF
    cat test/benchmark/config/site.env; echo
    python test/benchmark/rawmod_bench.py config | grep -E '^\s+(ok|--)\s+(BENCH_WORK|BENCH_DATA|DORADO|DORADO_MODELS|SAMTOOLS|POD5|RAWMOD_PYTHON|REFINE_PYTHON|BENCH_SLURM)'
    echo "STAGE config OK";;

  data)
    activate
    # the download script writes next to itself; keep the 3 GB on the lab share via a symlink
    [ -e test/data/d8_ont_oligos ] || ln -s "$WORK/data/d8_ont_oligos" test/data/d8_ont_oligos
    mkdir -p "$WORK/data/d8_ont_oligos"
    if [ "$(ls "$WORK/data/d8_ont_oligos/pod5_files"/*.pod5 2>/dev/null | wc -l)" -ge 8 ]; then
      echo "oligo pod5 files already present, download skipped"
    else
      (cd test/data && PATH="$(dirname "$SAMTOOLS"):$PATH:$ENV_REFINE/bin" bash download_d8_ont_oligos.sh)   # samtools + aws
    fi
    ls -la "$WORK/data/d8_ont_oligos/pod5_files" "$WORK/data/d8_ont_oligos/ref.fa"
    python test/benchmark/rawmod_bench.py datasets | grep -E 'sample|syn_'
    echo "STAGE data OK";;

  run)
    activate
    bash test/evaluation/1_run_benchmark.sh --samples syn_5mC syn_control --tools rawmod --slurm
    squeue -u "$USER" -o "%.10i %.20j %.8T %.10M %R" | grep -E 'JOBID|rb_' || true
    echo "STAGE run OK (jobs submitted; run 'status' until every column shows done)";;

  status)
    activate; echo "repo: $(git log --oneline -1)"
    python test/benchmark/rawmod_bench.py status | grep -E 'sample|syn_5mC|syn_control'
    echo "--- queue"; squeue -u "$USER" -o "%.10i %.20j %.8T %.10M %R" | grep -E 'JOBID|rb_' || echo "(no rb_ jobs queued)"
    echo "--- rb_ jobs of the last 12 h that did not finish cleanly (sacct)"
    sacct -X -u "$USER" -S "$(date -d '12 hours ago' +%Y-%m-%dT%H:%M)" --format=JobID,JobName%28,State%22,ExitCode,Elapsed 2>/dev/null \
      | grep -E '^[0-9]+ +rb_' | grep -vE 'COMPLETED|RUNNING|PENDING' || echo "(none)"
    echo "--- logs mentioning an error, newest first (logs in $WORK/work/logs/bench)"
    for f in $(grep -lE 'Error|Traceback|error:|FAILED' "$WORK"/work/logs/bench/* 2>/dev/null | xargs -r ls -t | head -3); do echo "== $f"; tail -n 15 "$f"; done
    echo "STAGE status OK";;

  score)
    activate
    # The score step runs the sequence-only nulls over EVERY row and crashes on rows whose dataset is not
    # downloaded (finding for Bhargav). Pre-build the nulls for our row under the three tags 2_score.sh uses
    # (region = Table 1, positions = supplementary, none = all sites) so score skips that step.
    (cd test/benchmark && for h in region positions ""; do BENCH_HOLDOUT=$h python -m bench.nulls --rows syn_5mC; done)
    ls "$WORK"/work/nulls_*.tsv
    bash test/evaluation/2_score.sh
    echo "--- summaries"; ls "$WORK"/work/results_*/rawmod/summary.md
    grep -hE 'syn_5mC|^\|.*row' "$WORK"/work/results_*/rawmod/summary.md | head -12
    echo "STAGE score OK";;

  *) sed -n '2,15p' "$0"; exit 1;;
esac
