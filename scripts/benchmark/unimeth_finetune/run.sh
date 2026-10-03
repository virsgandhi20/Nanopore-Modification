#!/bin/bash
# Fine-tune UniMeth for a chemistry it has no model for (Bhargav, Sep 26: every tool fills every row, so the 4mC and
# 5hmC columns need models). Same recipe as the 5hmU run: labels per read from matched samples, start from the
# published all-context 5mC checkpoint, held-out evaluation through the benchmark matrix.
#
#   CHEM=5hmC   positives = ONT all-5mers 5hmC oligo sample (the 256 designed C sites), control = the unmodified oligo
#               sample; hold-out = 20% of the reads of both samples (the matrix oligo rows then score every tool on
#               those reads)
#   CHEM=4mC    positives = H. pylori 26695 WT (4mC at the WT-vs-WGA differential sites), control = the WGA sample;
#               hold-out = the reads of the matrix's sub.bam region (the rest of the genome trains)
#
#   CHEM=5hmC MODE=setup|prep|train|status bash run.sh        (setup on the login node; prep CPU job; train GPU job)
#   Continue a preempted run:  CHEM=5hmC MODE=train RUN=<run dir> bash run.sh      FORCE=1 MODE=prep redoes the split + tags
# After training, run the matrix with TOOLS="unimeth_5hmC" (or unimeth_4mC): M_5hmC / M_4mC point at <run>/final.pt.
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd); HERE=$REPO/scripts/benchmark/unimeth_finetune; MX=$REPO/scripts/benchmark/matrix
ME=/fs/nexus-scratch/vgandhi; CAT=/fs/cbcb-lab/storm/shared/data; MATRIX=/fs/cbcb-lab/storm/vgandhi/matrix
CHEM=${CHEM:?set CHEM=5hmC, 4mC or 4mC_smrt}; W=${W:-$ME/unimeth_ft/$CHEM}; mkdir -p $W/{logs,status,runs,pod5,bam}
UM=$ME/Unimeth_5hmU                                     # the patched clone (5hmU patch + label patch); PYTHONPATH shadows the package
BASE_CKPT=$ME/unimeth_models/checkpoints/unimeth_r10.4.1_5kHz_5mC.pt
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
ENVACT="source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth; export PYTHONPATH=$UM"
SB="--account=scavenger --partition=scavenger --qos=scavenger --requeue"; GPU="--gres=gpu:rtxa5000:1"; CPU="--account=cbcb --partition=cbcb --qos=high"
[ "${PART:-}" = cbcb ] && { SB="--account=cbcb --partition=cbcb --qos=high"; GPU="--gres=gpu:1"; }
MAX_STEPS=${MAX_STEPS:-3000}; BATCH=${BATCH:-32}; VAL_READS=${VAL_READS:-800}; MODE=${MODE:-status}
case $CHEM in
5hmC) SYN=$CAT/synthetic_all5mers_r10.4.1_ont_open_data
      # the collection's reads_refined.bam carry RawMod-refined move tables that no longer match the basecalls (UniMeth:
      # SignalSequenceMismatchError), so the oligo samples use the matrix's own Dorado --emit-moves basecalls (MODE=basecall)
      POS_BAM=$MATRIX/syn_5hmC/moves.bam; POS_POD5=$SYN/pod5_files/5hmC_rep1.pod5
      NEG_BAM=$MATRIX/syn_control/moves.bam; NEG_POD5=$SYN/pod5_files/control_rep1.pod5
      SITES=$SYN/ground_truth/all_5mers_5hmC_sites.bed; CODE=h; FLAG="--hmC 1"; HOLDOUT=reads ;;
4mC)  HP=$CAT/hpylori_26695_wt_r10.4.1_ontbasemod_2024; HPW=$CAT/hpylori_26695_wga_r10.4.1_ontbasemod_2024
      POS_BAM=$HP/basecalled/reads.bam; POS_POD5=$HP/pod5_files; NEG_BAM=$HPW/basecalled/reads.bam; NEG_POD5=$HPW/pod5_files
      SITES=$ME/hp_labels/gt_4mC.bed; CODE=21839; FLAG="--m4C 1"; HOLDOUT=region ;;
4mC_smrt)  # same samples, labels from the SMRT methylome motif (GAAGA / TCTTC) instead of Dorado's differential calls
      HP=$CAT/hpylori_26695_wt_r10.4.1_ontbasemod_2024; HPW=$CAT/hpylori_26695_wga_r10.4.1_ontbasemod_2024
      POS_BAM=$HP/basecalled/reads.bam; POS_POD5=$HP/pod5_files; NEG_BAM=$HPW/basecalled/reads.bam; NEG_POD5=$HPW/pod5_files
      SITES=$ME/hp_labels/smrt_4mC_26695.bed; CODE=21839; FLAG="--m4C 1"; HOLDOUT=region ;;
*) echo "CHEM must be 5hmC, 4mC or 4mC_smrt"; exit 1 ;;
esac
COMMON="--pore_type R10.4.1 --frequency 4khz --cpg 1 --chg 1 --chh 1 $FLAG"
sub() { local id; id=$(sbatch --parsable $SB "$@") || { echo "sbatch failed: $*" >&2; echo ""; return; }; echo ${id%%;*}; }

case $MODE in
setup)
    eval "$ENVACT"
    [ -d $UM/.git ] || { echo "patched clone $UM missing: run unimeth_5hmU/run.sh MODE=setup first"; exit 1; }
    python $HERE/patch_labels.py $UM && python $REPO/scripts/benchmark/unimeth_5hmU/patch2_init_weights.py $UM
    python -m unimeth.training --help 2>&1 | grep -q -- "--hmC" && echo "training CLI has --hmC/--m4C" || { echo "FAIL: --hmC missing (clone does not import?)"; python -c "import unimeth.training.__main__" 2>&1 | tail -3; exit 1; }
    for f in $POS_BAM $POS_POD5 $NEG_BAM $NEG_POD5 $SITES $BASE_CKPT; do [ -e $f ] && echo "  ok $f" || echo "  MISSING $f"; done ;;

prep)
    : > $W/status/prep.txt
    [ -n "${FORCE:-}" ] && rm -rf $W/pod5/* $W/bam/*                     # FORCE=1: redo the read split and the tagged BAMs
    DEP=""; for d in syn_5hmC syn_control hp26695 hp26695_wga; do        # wait for a matrix basecall/prep of the source samples if one is queued
        j=$(awk -v k=basecall:$d -F'\t' '$1==k{print $2}' $MATRIX/jobs.tsv 2>/dev/null | tail -1); [ -n "$j" ] && [ -n "$(squeue -h -j $j 2>/dev/null)" ] && DEP="$DEP:$j"; done
    J=$(sub ${DEP:+--dependency=afterany$DEP} $CPU --cpus-per-task=8 --mem=48G --time=08:00:00 --job-name=ft_${CHEM}_prep --output=$W/logs/prep_%j.log --wrap="MODE=_prep CHEM=$CHEM W=$W bash $HERE/run.sh")
    echo "prep: job $J -> $W"; echo -e "prep\t$J" >> $W/jobs.tsv ;;

_prep)
    eval "$ENVACT"; set -x
    for S in pos neg; do
        if [ $S = pos ]; then BAM=$POS_BAM; POD5=$POS_POD5; STATE=modified; else BAM=$NEG_BAM; POD5=$NEG_POD5; STATE=unmodified; fi
        python $MX/pod5_ids.py $POD5 2>>$W/status/prep.txt | sort -u > $W/pod5/$S.pod5_ids.txt || exit 1
        $SAM view -F 0x904 $BAM | cut -f1 | sort -u | comm -12 - $W/pod5/$S.pod5_ids.txt > $W/pod5/$S.all_ids.txt
        if [ $HOLDOUT = reads ]; then
            shuf --random-source=<(yes 0) $W/pod5/$S.all_ids.txt > $W/pod5/$S.shuf.txt; n=$(wc -l < $W/pod5/$S.shuf.txt); t=$((n / 5))
            head -n $t $W/pod5/$S.shuf.txt | sort > $W/pod5/$S.test_ids.txt; tail -n +$((t + 1)) $W/pod5/$S.shuf.txt | sort > $W/pod5/$S.train_ids.txt
        else   # region hold-out: the reads of the matrix subset are the test set, everything else trains
            MSUB=$( [ $S = pos ] && echo $MATRIX/hp26695/sub.bam || echo $MATRIX/hp26695_wga/sub.bam )
            $SAM view $MSUB | cut -f1 | sort -u > $W/pod5/$S.test_ids.txt
            comm -23 $W/pod5/$S.all_ids.txt $W/pod5/$S.test_ids.txt > $W/pod5/$S.train_ids.txt
        fi
        head -n $VAL_READS $W/pod5/$S.train_ids.txt > $W/pod5/$S.val_ids.txt
        PFOPT=""; pod5 filter --help 2>/dev/null | grep -q -- "--missing-ok" && PFOPT="--missing-ok"
        if [ -d $POD5 ]; then PIN="-r $POD5"; else PIN="$POD5"; fi
        for part in train test val; do
            [ -s $W/pod5/$S.$part.pod5 ] || pod5 filter --ids $W/pod5/$S.${part}_ids.txt --output $W/pod5/$S.$part.pod5 $PFOPT -t 8 $PIN > $W/logs/pod5_${S}_$part.log 2>&1 || { echo "$S $part pod5 FAILED" >> $W/status/prep.txt; exit 1; }
        done
        [ -s $W/bam/$S.tagged.bam ] || python $HERE/add_site_tags.py --in $BAM --out $W/bam/$S.tagged.bam --base C --code $CODE --sites $SITES --state $STATE >> $W/status/prep.txt 2>&1 || { echo "$S tagging FAILED" >> $W/status/prep.txt; exit 1; }
        python -c "from unimeth.ioutils.reader import BamReader; BamReader('$W/bam/$S.tagged.bam')" > /dev/null 2>&1 && echo "$S: BAM read-id index OK" >> $W/status/prep.txt || { echo "$S: index FAILED" >> $W/status/prep.txt; exit 1; }
        echo "$S: $(wc -l < $W/pod5/$S.train_ids.txt) train / $(wc -l < $W/pod5/$S.test_ids.txt) test reads (of $(wc -l < $W/pod5/$S.all_ids.txt) with signal)" >> $W/status/prep.txt
    done
    echo "prep finished $(date)" >> $W/status/prep.txt ;;

train)
    JP=$(awk -F'\t' '$1=="prep"{print $2}' $W/jobs.tsv 2>/dev/null | tail -1); [ -n "$JP" ] && [ -z "$(squeue -h -j $JP 2>/dev/null)" ] && JP=""
    RUN=${RUN:-$W/runs/${CHEM}_$(date +%m%d_%H%M)}; mkdir -p $RUN
    CK=$( { ls -d $RUN/checkpoint-* 2>/dev/null || true; } | sed 's/.*checkpoint-//' | sort -n | tail -1 )
    STEPS=$MAX_STEPS; OUT=$RUN; INITENV="UNIMETH_RESUME=0"; NOTE="fresh run from $(basename $BASE_CKPT)"
    if [ -n "$CK" ]; then
        INIT=$(ls $RUN/checkpoint-$CK/pytorch_model.bin $RUN/checkpoint-$CK/model.safetensors 2>/dev/null | head -1); [ -n "$INIT" ] || { echo "checkpoint-$CK has no weights"; exit 1; }
        STEPS=$((MAX_STEPS - CK)); [ $STEPS -gt 0 ] || { echo "already at MAX_STEPS"; exit 1; }
        OUT=$W/runs/$(basename $RUN)_from$CK; mkdir -p $OUT; INITENV="UNIMETH_INIT_WEIGHTS=$INIT UNIMETH_RESUME=0"; NOTE="continuing from $INIT for $STEPS more steps"
    fi
    : > $W/status/train.txt; echo "$NOTE -> $OUT" | tee -a $W/status/train.txt
    J=$(sub ${JP:+--dependency=afterany:$JP} $GPU --cpus-per-task=8 --mem=64G --time=08:00:00 --job-name=ft_${CHEM}_train --output=$W/logs/train_%j.log \
          --wrap="MODE=_train CHEM=$CHEM W=$W OUT=$OUT STEPS=$STEPS INITENV='$INITENV' bash $HERE/run.sh")
    echo "train: job $J -> $OUT"; echo -e "train\t$J\t$OUT" >> $W/jobs.tsv ;;

_train)
    eval "$ENVACT"; cd $OUT
    for f in $W/bam/pos.tagged.bam $W/bam/neg.tagged.bam $W/pod5/pos.train.pod5 $W/pod5/neg.train.pod5 $W/pod5/pos.val.pod5 $W/pod5/neg.val.pod5; do [ -s $f ] || { echo "train: missing $f" >> $W/status/train.txt; exit 1; }; done
    echo "train started $(date) run=$OUT steps=$STEPS batch=$BATCH" >> $W/status/train.txt
    env $INITENV PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True UNIMETH_OUT_DIR=$OUT UNIMETH_LOG_STEPS=50 UNIMETH_EVAL_STEPS=500 UNIMETH_SAVE_STEPS=500 UNIMETH_DL_WORKERS=6 \
        python -m unimeth.training --mode finetune --bam_dir $W/bam/pos.tagged.bam,$W/bam/neg.tagged.bam --train_pod5_dir $W/pod5/pos.train.pod5,$W/pod5/neg.train.pod5 \
        --val_pod5_dir $W/pod5/pos.val.pod5,$W/pod5/neg.val.pod5 --model_dir $BASE_CKPT $COMMON --dorado_version 1.4 --max_steps $STEPS --batch_size $BATCH --run_name $CHEM > $OUT/train.log 2>&1
    rc=$?
    if [ -s $OUT/final.pt ]; then ln -sfn $OUT $W/runs/latest; echo "train OK $(date): $OUT/final.pt (runs/latest)" >> $W/status/train.txt
    else echo "train FAILED rc=$rc: $(grep -iE 'error|Traceback' $OUT/train.log | tail -2 | tr '\n' ' ')" >> $W/status/train.txt; fi ;;

status)
    squeue -u $USER -o "%.9i %.18j %.3t %.9M %R" | grep -E "ft_${CHEM}|JOBID"; echo
    for f in prep train; do [ -s $W/status/$f.txt ] && { echo "== $f"; cat $W/status/$f.txt; echo; }; done
    RUN=$(awk -F'\t' '$1=="train"{print $3}' $W/jobs.tsv 2>/dev/null | tail -1); [ -n "$RUN" ] && [ -s $RUN/train.log ] && { echo "== training progress"; grep -E "^\{'loss'|eval_" $RUN/train.log | tail -4 | cut -c1-230; } ;;
*) echo "MODE must be setup | prep | train | status"; exit 1 ;;
esac
