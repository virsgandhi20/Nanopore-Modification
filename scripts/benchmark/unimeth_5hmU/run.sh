#!/bin/bash
# UniMeth fine-tuned for 5hmU on SPO1 (Bhargav, 2026-09-22: "run unimeth fine tuned on 5hmU").
#
# Labels: SPO1 native DNA carries 5hmU at every T (barcodes 06, 07); the PCR
# amplicon libraries of the same genome carry none (barcodes 01-05), so every T
# in a read is labelled by its sample. Train on bc06 (modified) + bc02
# (unmodified); test on bc07 + bc01, which the model never sees (the RawMod
# draft also holds out barcode01). Start from the published 6mA checkpoint,
# whose [m6A] token row seeds the new [5hmU] token.
#
# Modes (login node unless stated):
#   setup     clone UniMeth into its own directory and apply the 5hmU patch
#   preflight checks + a 3-step finetune on 200 reads (GPU job, ~10 min)
#   prep      submits: bc01 basecalling (GPU) -> tagged BAMs, val pod5s, BAM
#             indexes (CPU)
#   train     submits the finetune (GPU, depends on prep)
#   eval      submits inference on bc07 + bc01 and scoring (GPU, depends on train)
#   status    what finished, with the numbers
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)
HERE=$REPO/scripts/benchmark/unimeth_5hmU
ME=/fs/nexus-scratch/vgandhi
W=${W:-$ME/unimeth_5hmU}                      # everything this makes
UM=$ME/Unimeth_5hmU                           # patched clone (shadows the installed package via PYTHONPATH)
DATA=/fs/cbcb-lab/storm/shared/umbc-ont-data
POD5=$DATA/pod5_by_barcode/run1_jan31/single_end
REF=$DATA/ref/SPO1_FJ230960.1.fasta
SPO1=$ME/spo1_typing                          # bc02/06/07 move-table BAMs from the typing work (Sep 19)
BASE_CKPT=$ME/unimeth_models/checkpoints/unimeth_r10.4.1_5kHz_6mA.pt
DORADO=/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-1.4.0-linux-x64/bin/dorado
DMODEL=$ME/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
ENVACT="source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth; export PYTHONPATH=$UM"
SB="--account=scavenger --partition=scavenger --qos=scavenger --requeue"
GPU="--gres=gpu:rtxa5000:1"
MODE=${MODE:-status}
MAX_STEPS=${MAX_STEPS:-3000}; BATCH=${BATCH:-32}; VAL_READS=${VAL_READS:-800}; TEST_LIMIT=${TEST_LIMIT:-15000}
mkdir -p $W/{bam,pod5,logs,status,runs,eval} 2>/dev/null
declare -A BAM=( [bc06]=$SPO1/barcode06.mod.sorted.bam [bc02]=$SPO1/barcode02.mod.sorted.bam [bc07]=$SPO1/barcode07.mod.sorted.bam [bc01]=$W/bam/barcode01.moves.bam )
declare -A STATE=( [bc06]=modified [bc02]=unmodified [bc07]=modified [bc01]=unmodified )
COMMON="--pore_type R10.4.1 --frequency 4khz --hmU 1"      # 4khz = the calibrated normalization path (see the UniMeth bug write-up)
sub() { local id; id=$(sbatch --parsable $SB "$@") || { echo "sbatch failed: $*" >&2; echo ""; return; }; echo ${id%%;*}; }
note() { echo "$*" | tee -a $W/status/$MODE.txt; }

case $MODE in
setup)
    eval "$ENVACT"
    [ -d $UM/.git ] || git clone -q --depth 30 https://github.com/sekeyWang/Unimeth.git $UM
    (cd $UM && git checkout -q 11215d4 2>/dev/null; git status --short | grep -q . && { echo "clone already modified; resetting"; git checkout -q -- .; })
    python $HERE/patch_unimeth_5hmU.py $UM
    python -c "import unimeth, unimeth.config as c; print('patched clone imports from', unimeth.__file__); print('vocab', len(c.VOCAB), c.VOCAB[-1])"
    python -m unimeth.training --help 2>&1 | grep -q -- "--hmU" && echo "training CLI has --hmU" || echo "FAIL: --hmU missing from training CLI"
    python -m unimeth.inference --help 2>&1 | grep -q -- "--hmU" && echo "inference CLI has --hmU" || echo "FAIL: --hmU missing from inference CLI" ;;

preflight)
    eval "$ENVACT"; FAILS=0; ok() { echo "  PASS  $1"; }; bad() { echo "  FAIL  $1"; FAILS=$((FAILS+1)); }
    python -c "import unimeth.config as c; assert c.VOCAB[-1] == '[5hmU]'" 2>/dev/null && ok "patched clone on PYTHONPATH" || bad "patched clone missing: run MODE=setup"
    for b in bc06 bc02 bc07; do [ -s ${BAM[$b]} ] && ok "BAM $b" || bad "BAM ${BAM[$b]}"; done
    for b in 01 02 06 07; do [ -r $POD5/barcode$b.pod5 ] && ok "pod5 barcode$b ($(du -h $POD5/barcode$b.pod5 | cut -f1))" || bad "pod5 barcode$b"; done
    [ -s $BASE_CKPT ] && ok "base checkpoint (6mA)" || bad "base checkpoint $BASE_CKPT"
    [ -r $REF ] && ok "reference" || bad "reference $REF"
    NMV=$($SAM view ${BAM[bc06]} 2>/dev/null | head -200 | grep -c "mv:B" || true)   # grep -q under pipefail false-fails on a hit
    [ "${NMV:-0}" -gt 0 ] && ok "bc06 BAM has move tables ($NMV of first 200 records)" || bad "bc06 BAM lacks mv tags"
    FREE=$(df --output=avail -BG $ME | tail -1 | tr -dc 0-9); [ "${FREE:-0}" -ge 8 ] && ok "free space ${FREE}G" || bad "only ${FREE}G free (need ~8G)"
    echo "--- tiny real run as a GPU job (200 reads each, 3 steps): $W/pf"
    mkdir -p $W/pf
    J=$(sub $GPU --cpus-per-task=6 --mem=48G --time=00:40:00 --job-name=hmu_pf --output=$W/logs/pf_%j.log --wrap="$ENVACT; set -uo pipefail; cd $W/pf
      python $HERE/add_5hmU_tags.py --in ${BAM[bc06]} --out $W/pf/bc06.bam --state modified --limit 400 && python $HERE/add_5hmU_tags.py --in ${BAM[bc02]} --out $W/pf/bc02.bam --state unmodified --limit 400 || exit 1
      python - <<'PY'
import pysam, pod5
from unimeth.ioutils.reader import BamReader
for b, src in (('bc06', '$POD5/barcode06.pod5'), ('bc02', '$POD5/barcode02.pod5')):
    ids = [r.query_name for r in pysam.AlignmentFile(f'$W/pf/{b}.bam')]
    with pod5.Reader(src) as rd, pod5.Writer(f'$W/pf/{b}_train.pod5') as wt, pod5.Writer(f'$W/pf/{b}_val.pod5') as wv:
        for i, rec in enumerate(rd.reads(selection=ids, missing_ok=True)):
            (wv if i % 4 == 0 else wt).add_read(rec.to_read())
    BamReader(f'$W/pf/{b}.bam'); print(b, len(ids), 'reads prepared')
PY
      UNIMETH_OUT_DIR=$W/pf/run UNIMETH_LOG_STEPS=1 UNIMETH_EVAL_STEPS=100 UNIMETH_SAVE_STEPS=100 UNIMETH_DL_WORKERS=2 python -m unimeth.training --mode finetune --bam_dir $W/pf/bc06.bam,$W/pf/bc02.bam --train_pod5_dir $W/pf/bc06_train.pod5,$W/pf/bc02_train.pod5 --val_pod5_dir $W/pf/bc06_val.pod5,$W/pf/bc02_val.pod5 --model_dir $BASE_CKPT $COMMON --dorado_version 1.4 --max_steps 3 --batch_size 8 --run_name pf && ls -la $W/pf/run/final.pt
      python -m unimeth.inference --pod5 $W/pf/bc06_val.pod5 --bam $W/pf/bc06.bam --model $W/pf/run/final.pt $COMMON --output_format tsv --out $W/pf/calls.txt --num_workers 2 && echo \"preflight inference rows: \$(wc -l < $W/pf/calls.txt) types: \$(cut -f7 $W/pf/calls.txt | sort | uniq -c | tr '\\n' ' ')\" && echo PREFLIGHT_JOB_OK")
    [ -n "$J" ] && ok "preflight GPU job submitted: $J (check: grep -E 'PREFLIGHT_JOB_OK|Error|error' $W/logs/pf_$J.log)" || bad "could not submit preflight job"
    echo; [ $FAILS -eq 0 ] && echo "PREFLIGHT CHECKS: ALL PASS (now wait for job $J)" || echo "PREFLIGHT: $FAILS FAILURE(S)"; exit $FAILS ;;

prep)
    : > $W/status/prep.txt; PRE="$ENVACT; set -uo pipefail; mkdir -p $W/status"
    # bc01 has never been basecalled with move tables; GPU
    JBC=""
    if [ ! -s ${BAM[bc01]} ]; then
        JBC=$(sub $GPU --cpus-per-task=8 --mem=48G --time=04:00:00 --job-name=hmu_bc01 --output=$W/logs/bc01_%j.log --wrap="$PRE; $DORADO basecaller $DMODEL $POD5/barcode01.pod5 --emit-moves --reference $REF > $W/bam/barcode01.unsorted.bam && $SAM sort -@ 8 -o ${BAM[bc01]} $W/bam/barcode01.unsorted.bam && $SAM index ${BAM[bc01]} && rm -f $W/bam/barcode01.unsorted.bam && echo \"bc01 basecall OK: \$($SAM flagstat ${BAM[bc01]} | grep -m1 'mapped (')\" >> $W/status/prep.txt || echo 'bc01 basecall FAILED' >> $W/status/prep.txt")
        echo "basecall bc01: job $JBC"
    fi
    # tags + val pod5 + indexes; CPU; bc01 part waits for the basecall
    JP=$(sub ${JBC:+--dependency=afterany:$JBC} --cpus-per-task=4 --mem=32G --time=06:00:00 --job-name=hmu_prep --output=$W/logs/prep_%j.log --wrap="$PRE
      for b in bc06 bc02 bc07 bc01; do
        src=\$(case \$b in bc06) echo ${BAM[bc06]};; bc02) echo ${BAM[bc02]};; bc07) echo ${BAM[bc07]};; bc01) echo ${BAM[bc01]};; esac)
        st=\$(case \$b in bc06|bc07) echo modified;; *) echo unmodified;; esac)
        [ -s \$src ] || { echo \"\$b: source BAM missing, skipped\" >> $W/status/prep.txt; continue; }
        [ -s $W/bam/\$b.tagged.bam ] || python $HERE/add_5hmU_tags.py --in \$src --out $W/bam/\$b.tagged.bam --state \$st >> $W/status/prep.txt 2>&1 || echo \"\$b: tagging FAILED\" >> $W/status/prep.txt
        python -c \"from unimeth.ioutils.reader import BamReader; BamReader('$W/bam/\$b.tagged.bam')\" > /dev/null 2>&1 && echo \"\$b: BAM read-id index OK\" >> $W/status/prep.txt || echo \"\$b: index FAILED\" >> $W/status/prep.txt
      done
      for b in bc06 bc02; do
        num=\$(case \$b in bc06) echo 06;; bc02) echo 02;; esac)
        [ -s $W/pod5/\${b}_val.pod5 ] || python - <<PY >> $W/status/prep.txt 2>&1
import numpy as np, pod5
src = '$POD5/barcode\$num.pod5'
with pod5.Reader(src) as rd: ids = [str(x) for x in rd.read_ids]
sel = set(np.random.default_rng(0).choice(ids, min($VAL_READS, len(ids)), replace=False).tolist())
n = 0
with pod5.Reader(src) as rd, pod5.Writer('$W/pod5/\${b}_val.pod5') as wv:
    for rec in rd.reads(selection=list(sel)): wv.add_read(rec.to_read()); n += 1
print('\$b: val pod5 with', n, 'of', len(ids), 'reads (training streams the full pod5; these', n, 'reads overlap it)')
PY
      done; echo 'prep finished' >> $W/status/prep.txt")
    echo "prep: job $JP"; echo -e "bc01\t$JBC\nprep\t$JP" > $W/jobs_prep.tsv ;;

train)
    JP=$(awk -F'\t' '$1=="prep"{print $2}' $W/jobs_prep.tsv 2>/dev/null)
    [ -n "$JP" ] && [ -z "$(squeue -h -j $JP 2>/dev/null)" ] && JP=""      # prep already finished and left the queue
    RUN=${RUN:-$W/runs/hmu_$(date +%m%d_%H%M)}; mkdir -p $RUN           # RUN=<existing run dir> continues from its last checkpoint-* (weights only, see below)
    python3 $HERE/patch2_init_weights.py $UM || exit 1                    # idempotent; adds UNIMETH_INIT_WEIGHTS to the patched clone
    # The Trainer's own resume (optimizer/rng state via torch.load) is refused by transformers on the cluster's torch 2.5.1,
    # so a continuation loads the last checkpoint's weights, restarts optimizer + LR schedule, runs the REMAINING steps in a new dir.
    CK=$( { ls -d $RUN/checkpoint-* 2>/dev/null || true; } | sed 's/.*checkpoint-//' | sort -n | tail -1 )
    STEPS=$MAX_STEPS; OUT=$RUN; INITENV="UNIMETH_RESUME=0"; NOTE="fresh run from $(basename $BASE_CKPT)"
    if [ -n "$CK" ]; then
        INIT=$(ls $RUN/checkpoint-$CK/pytorch_model.bin $RUN/checkpoint-$CK/model.safetensors 2>/dev/null | head -1)
        [ -n "$INIT" ] || { echo "checkpoint-$CK has no weights file"; exit 1; }
        STEPS=$((MAX_STEPS - CK)); [ $STEPS -gt 0 ] || { echo "checkpoint-$CK already reached MAX_STEPS=$MAX_STEPS"; exit 1; }
        OUT=$W/runs/$(basename $RUN)_from$CK; mkdir -p $OUT
        INITENV="UNIMETH_INIT_WEIGHTS=$INIT UNIMETH_RESUME=0"; NOTE="continuing from $INIT for $STEPS more steps (optimizer and LR schedule restart)"
    fi
    : > $W/status/train.txt; echo "$NOTE -> $OUT" | tee -a $W/status/train.txt
    J=$(sub ${JP:+--dependency=afterany:$JP} $GPU --cpus-per-task=8 --mem=64G --time=08:00:00 --job-name=hmu_train --output=$W/logs/train_%j.log --wrap="$ENVACT; set -uo pipefail; cd $OUT
      for f in $W/bam/bc06.tagged.bam $W/bam/bc02.tagged.bam $W/pod5/bc06_val.pod5 $W/pod5/bc02_val.pod5; do [ -s \$f ] || { echo \"train: missing input \$f\" >> $W/status/train.txt; exit 1; }; done
      echo \"train started \$(date) run=$OUT steps=$STEPS batch=$BATCH\" >> $W/status/train.txt
      PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True UNIMETH_OUT_DIR=$OUT $INITENV UNIMETH_LOG_STEPS=50 UNIMETH_EVAL_STEPS=500 UNIMETH_SAVE_STEPS=500 UNIMETH_DL_WORKERS=6 python -m unimeth.training --mode finetune --bam_dir $W/bam/bc06.tagged.bam,$W/bam/bc02.tagged.bam --train_pod5_dir $POD5/barcode06.pod5,$POD5/barcode02.pod5 --val_pod5_dir $W/pod5/bc06_val.pod5,$W/pod5/bc02_val.pod5 --model_dir $BASE_CKPT $COMMON --dorado_version 1.4 --max_steps $STEPS --batch_size $BATCH --run_name hmu > $OUT/train.log 2>&1
      rc=\$?; [ -s $OUT/final.pt ] && echo \"train OK \$(date): $OUT/final.pt; last eval: \$(grep -o \"'eval_\\[5hmU\\]': {[^}]*}\" $OUT/train.log | tail -1 | cut -c1-200)\" >> $W/status/train.txt || echo \"train FAILED rc=\$rc: \$(grep -iE 'error|Traceback' $OUT/train.log | tail -2 | tr '\\n' ' ')\" >> $W/status/train.txt")
    echo "train: job $J -> $OUT"; echo -e "train\t$J\t$OUT" >> $W/jobs_train.tsv ;;

eval)
    JT=$(tail -1 $W/jobs_train.tsv | cut -f2); RUN=$(tail -1 $W/jobs_train.tsv | cut -f3); E=$W/eval/$(basename $RUN); mkdir -p $E; : > $W/status/eval.txt
    [ -n "$JT" ] && [ -z "$(squeue -h -j $JT 2>/dev/null)" ] && JT=""
    J=$(sub ${JT:+--dependency=afterany:$JT} $GPU --cpus-per-task=8 --mem=48G --time=04:00:00 --job-name=hmu_eval --output=$W/logs/eval_%j.log --wrap="$ENVACT; set -uo pipefail
      [ -s $RUN/final.pt ] || { echo 'eval: no final.pt' >> $W/status/eval.txt; exit 1; }
      for b in bc07 bc01; do num=\$(case \$b in bc07) echo 07;; bc01) echo 01;; esac
        [ -s $E/\$b.calls.txt ] || python -m unimeth.inference --pod5 $POD5/barcode\$num.pod5 --bam $W/bam/\$b.tagged.bam --model $RUN/final.pt $COMMON --output_format tsv --out $E/\$b.calls.txt --num_workers 8 --limit $TEST_LIMIT > $E/\$b.infer.log 2>&1 || { echo \"eval: inference \$b FAILED: \$(grep -iE 'error' $E/\$b.infer.log | tail -1)\" >> $W/status/eval.txt; exit 1; }
        python $ME/Unimeth/scripts/call_modification_frequency.py -i $E/\$b.calls.txt -o $E/\$b.sites.tsv --sort
      done
      python $HERE/score_5hmU.py --pos-calls $E/bc07.calls.txt --neg-calls $E/bc01.calls.txt --pos-sites $E/bc07.sites.tsv --neg-sites $E/bc01.sites.tsv --scorer $REPO/scripts/benchmark/score_sites.py --out $E/result.txt >> $W/status/eval.txt 2>&1
      echo \"eval finished \$(date)\" >> $W/status/eval.txt")
    echo "eval: job $J -> $E" ;;

status)
    squeue -u $USER -o "%.9i %.14j %.3t %.9M %R" | grep -E "hmu_|JOBID"; echo
    for f in prep train eval; do [ -s $W/status/$f.txt ] && { echo "== $f"; cat $W/status/$f.txt; echo; }; done
    RUN=$(tail -1 $W/jobs_train.tsv 2>/dev/null | cut -f3); [ -n "$RUN" ] && [ -s $RUN/train.log ] && { echo "== training progress"; grep -E "^\{'loss'|eval_\[5hmU\]" $RUN/train.log | tail -4 | cut -c1-230; }
    ls $W/eval/*/result.txt 2>/dev/null | while read f; do echo "== $f"; cat $f; done ;;
*) echo "MODE must be setup | preflight | prep | train | eval | status";;
esac
