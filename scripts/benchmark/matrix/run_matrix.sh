#!/bin/bash
# Every tool on every row of Table 1 (Bhargav, Sep 25-26: no N/A cells; a tool with no model for a row's chemistry
# scores 0.5 through the fill rule). Samples come from datasets.tsv, rows from rows.tsv, both next to this script.
#
#   MODE=check  bash run_matrix.sh                 # login node: every path exists, ground truth sits on the right base
#   MODE=prep   [DS="ecoli_wt anabaena"] ...       # CPU jobs: reference copy, sorted BAM, first-NREADS subset (sub.bam);
#                                                  # a BAM without move tables is first re-basecalled (GPU job, --emit-moves)
#   MODE=infer  [TOOLS="unimeth_6mA"] [DS=...]     # GPU jobs, one per tool x sample, wait for that sample's prep; a scoring job is queued behind them (NOSCORE=1 to skip)
#   MODE=score                                     # CPU job: score_matrix.py over everything that has finished (waits for queued cells; NOWAIT=1 to score now)
#   MODE=all                                       # prep + infer for every sample and tool, then score
#   MODE=status                                    # queue, what finished, the grid
#   MODE=audit                                     # reads seen by each finished tool vs reads in the subset (catches truncated outputs)
#   FORCE=1 with prep or infer redoes a finished sample / tool (prep also re-applies the pod5 read filter); with basecall, redoes the BAM
#   NOTRIM=1 with basecall passes --no-trim to Dorado (needed for the short oligo reads, see _basecall)
#   PART=cbcb submits to the lab partition (qos high, any GPU) instead of scavenger; UniMeth then runs with --batch_size 64
# Tools: unimeth_5mC (all-context model), unimeth_6mA, unimeth_5hmU / unimeth_5hmC / unimeth_4mC (fine-tuned, patched clone), deepmod2 (CpG),
# rockfish (CpG; wired once its output format is known). Plant/human samples reuse the sub.bam, UniMeth 5mC and
# DeepMod2 outputs already produced under $EUK by run_catalog_rows.sh (datasets.tsv `reuse`).
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd); HERE=$REPO/scripts/benchmark/matrix
export CAT=/fs/cbcb-lab/storm/shared/data ME=/fs/nexus-scratch/vgandhi EUK=/fs/cbcb-lab/storm/vgandhi/euk UMBC=/fs/cbcb-lab/storm/shared/umbc-ont-data
W=${WBASE:-/fs/cbcb-lab/storm/vgandhi/matrix}; mkdir -p $W/logs $W/status $W/rows
ALLTOOLS="unimeth_5mC unimeth_6mA unimeth_5hmU unimeth_5hmC unimeth_4mC deepmod2 rockfish"     # every tool a scoring pass covers (missing outputs stay pending)
TOOLS=${TOOLS:-"unimeth_5mC unimeth_6mA unimeth_5hmU deepmod2"}
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
ENVACT="source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth033"   # UniMeth 0.3.3 (normalization fix, their issue 21); built by unimeth_v033_check.sh MODE=setup
UM_PATCHED=$ME/Unimeth_0.3.3_patched                      # v0.3.3 clone + patch_unimeth_033_infer.py (18-token vocabulary of the fine-tuned checkpoints, --hmU); shadows the package via PYTHONPATH
MODELS=$ME/unimeth_models/checkpoints
M_5mC=$MODELS/unimeth_r10.4.1_5kHz_5mC.pt; M_6mA=$MODELS/unimeth_r10.4.1_5kHz_6mA.pt
M_5hmU=${M_5hmU:-$ME/unimeth_5hmU/runs/hmu_0922_1446_from1000_from1500/final.pt}     # the 3,000-step model (0.7725 held out)
M_5hmC=${M_5hmC:-$ME/unimeth_ft/5hmC/runs/latest/final.pt}; M_4mC=${M_4mC:-$ME/unimeth_ft/4mC/runs/latest/final.pt}   # unimeth_finetune/run.sh
DORADO=/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-1.4.0-linux-x64/bin/dorado; DMODEL=$ME/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
RF_ENV=$ME/envs/rockfish; RF_MODEL=${RF_MODEL:-$ME/rockfish_bench/models/rf_5kHz.ckpt}; RF_ORIENT=${RF_ORIENT:-read}; RF_SHIFT=${RF_SHIFT:-0}   # Rockfish positions index the read as sequenced (M.SssI check, Sep 26): plus calls land on the C, minus calls on the G = the minus-strand C
DM2_ENV=$ME/envs/deepmod2; DM2_SRC=$ME/deepmod2_bench/DeepMod2; DM2_MODEL=${DM2_MODEL:-bilstm_r10.4.1_5khz_v5.0}
SB="--account=scavenger --partition=scavenger --qos=scavenger --requeue"; GPU="--gres=gpu:rtxa5000:1"
CPU="--account=cbcb --partition=cbcb --qos=high"
UM_BATCH=${UM_BATCH:-256}                                                   # UniMeth inference batch (default 256 fits the 24 GB A5000s)
[ "${PART:-}" = cbcb ] && { SB="--account=cbcb --partition=cbcb --qos=high"; GPU="--gres=gpu:1"; UM_BATCH=${UM_BATCH_CBCB:-64}; }   # PART=cbcb: lab partition (11 GB cards: UniMeth OOMs at 256, 64 fits)
MODE=${MODE:-status}

# ---- tables -> bash arrays (paths expanded by matrix_common.py so both sides read them the same way)
declare -A GTDIR POD5 BAM REF NREADS REUSE RFROM; DSIDS=()
while IFS=$'\t' read -r id gtdir pod5 bam ref nreads reuse rfrom; do
    [ -n "$id" ] && [ "${id:0:1}" != "#" ] || continue
    DSIDS+=($id); GTDIR[$id]=$gtdir; POD5[$id]=$pod5; BAM[$id]=$bam; REF[$id]=$ref; NREADS[$id]=$nreads; REUSE[$id]=$reuse; RFROM[$id]=${rfrom:--}
done < <(python -c "
import sys; sys.path.insert(0, '$HERE'); from matrix_common import load_datasets
for d in load_datasets('$HERE/datasets.tsv').values(): print('\t'.join(str(d[k]) for k in ('id','gtdir','pod5','bam','ref','nreads','reuse','region_from')))")
[ ${#DSIDS[@]} -gt 0 ] || { echo "datasets.tsv could not be read (python + matrix_common.py?)"; exit 1; }
DS=${DS:-${DSIDS[*]}}

sub() { local id; id=$(sbatch --parsable $SB "$@") || { echo "sbatch failed: $*" >&2; echo ""; return; }; echo ${id%%;*}; }
job_of() { awk -v k="$1" -F'\t' '$1==k{print $2}' $W/jobs.tsv 2>/dev/null | tail -1; }
dep() { local j; j=$(job_of "$1"); [ -n "$j" ] && [ -n "$(squeue -h -j $j 2>/dev/null)" ] && echo "--dependency=afterany:$j"; }
record() { echo -e "$1\t$2" >> $W/jobs.tsv; }
big_ref() { [ $(stat -Lc %s ${REF[$1]}) -gt 1000000000 ]; }                 # > 1 Gb reference: DeepMod2 workers each hold it
has_moves() { [ $($SAM view $1 2>/dev/null | head -200 | grep -c 'mv:B') -gt 0 ]; }   # first 200 records carry move tables?

case $MODE in
check)
    eval "$ENVACT" 2>/dev/null; python $HERE/check_rows.py --datasets $HERE/datasets.tsv --rows $HERE/rows.tsv --work $W --euk $EUK ${DETAIL:+--detail $DETAIL} ;;

basecall)   # GPU: Dorado sup v5 with --emit-moves, for samples whose collection BAM has no move tables (the oligo set)
    for d in $DS; do
        D=$W/$d; mkdir -p $D; [ -n "${FORCE:-}" ] && rm -f $D/moves.bam $D/moves.bam.bai $D/status/prep.done
        [ -s $D/moves.bam ] && { echo "basecall $d: $D/moves.bam exists"; continue; }
        J=$(sub $GPU --cpus-per-task=8 --mem=48G --time=08:00:00 --job-name=mtx_bc_$d --output=$W/logs/basecall_${d}_%j.log \
              --wrap="MODE=_basecall DSID=$d WBASE=$W NOTRIM=${NOTRIM:-} bash $HERE/run_matrix.sh")
        [ -n "$J" ] && { record basecall:$d $J; echo "basecall $d: job $J"; }
    done ;;

prep)
    for d in $DS; do
        D=$W/$d; mkdir -p $D
        [ -n "${FORCE:-}" ] && rm -f $D/status/prep.done $D/status/pod5_filtered $D/reads.sorted.bam $D/reads.sorted.bam.bai $D/sub.bam $D/sub.bam.bai $D/pod5_ids.txt $D/reads.inpod5.bam   # rebuild from the current source BAM
        [ -s $D/status/prep.done ] 2>/dev/null && { echo "prep $d: done already"; continue; }
        if [ ! -s $D/moves.bam ] && [ -z "$(job_of basecall:$d)" ] && ! has_moves ${BAM[$d]}; then
            echo "prep $d: source BAM has no move tables, submitting a basecall first"; MODE=basecall DS=$d WBASE=$W bash $HERE/run_matrix.sh; fi
        J=$(sub $(dep basecall:$d) $CPU --cpus-per-task=8 --mem=32G --time=06:00:00 --job-name=mtx_prep_$d --output=$W/logs/prep_${d}_%j.log \
              --wrap="MODE=_prep DSID=$d WBASE=$W bash $HERE/run_matrix.sh")
        [ -n "$J" ] && { record prep:$d $J; echo "prep $d: job $J"; }
    done ;;

infer)
    for t in $TOOLS; do for d in $DS; do
        T=$W/$d/$t; [ -n "${FORCE:-}" ] && rm -rf $T; mkdir -p $T
        [ -s $T/sites.std.tsv ] && { echo "$t on $d: done already"; continue; }
        case $t in
            deepmod2) if big_ref $d; then RES="--cpus-per-task=4 --mem=120G"; else RES="--cpus-per-task=12 --mem=48G"; fi ;;
            *)        RES="--cpus-per-task=8 --mem=48G" ;;
        esac
        J=$(sub $(dep prep:$d) $GPU $RES --time=08:00:00 --job-name=mtx_${t}_$d --output=$W/logs/${t}_${d}_%j.log \
              --wrap="MODE=_infer DSID=$d TOOL=$t WBASE=$W RF_ORIENT=$RF_ORIENT RF_SHIFT=$RF_SHIFT UM_BATCH=$UM_BATCH bash $HERE/run_matrix.sh")
        [ -n "$J" ] && { record infer:$t:$d $J; echo "$t on $d: job $J $(dep prep:$d)"; NEWJOBS=1; }
    done; done
    # every inference submission queues a scoring pass behind it (cells that finish with no scorer waiting stay unscored)
    [ -n "${NEWJOBS:-}" ] && [ -z "${NOSCORE:-}" ] && MODE=score TOOLS="$ALLTOOLS" DS="${DSIDS[*]}" WBASE=$W bash $HERE/run_matrix.sh ;;

score)
    DEPS=$(for t in $TOOLS; do for d in $DS; do j=$(job_of infer:$t:$d); [ -n "$j" ] && [ -n "$(squeue -h -j $j 2>/dev/null)" ] && echo -n ":$j"; done; done)
    [ -n "${NOWAIT:-}" ] && DEPS=""                                          # NOWAIT=1: score what is finished now; pending cells stay pending
    J=$(sub ${DEPS:+--dependency=afterany${DEPS}} $CPU --cpus-per-task=4 --mem=48G --time=06:00:00 --job-name=mtx_score --output=$W/logs/score_%j.log \
          --wrap="MODE=_score WBASE=$W TOOLS='$TOOLS' CPGMERGE='${CPGMERGE:-}' bash $HERE/run_matrix.sh")
    [ -n "$J" ] && { record score $J; echo "score: job $J${DEPS:+ (after$DEPS)}"; } ;;

all)
    MODE=prep bash $HERE/run_matrix.sh; MODE=infer TOOLS="$TOOLS" DS="$DS" bash $HERE/run_matrix.sh; MODE=score TOOLS="$TOOLS" DS="$DS" bash $HERE/run_matrix.sh ;;

audit)   # reads seen by each finished tool vs reads in the subset: a cell far below 100% is a truncated output (preempted job)
    printf "%-12s %8s" sample reads; for t in $TOOLS; do printf " %-13s" $t; done; echo
    for d in $DS; do
        [ -s $W/$d/sub.bam ] || continue; nb=$($SAM view $W/$d/sub.bam | cut -f1 | sort -u | wc -l); printf "%-12s %8s" $d $nb
        for t in $TOOLS; do
            if [ ! -s $W/$d/$t/sites.std.tsv ]; then s=-
            else case $t in
                unimeth_*) n=$(cut -f5 $W/$d/$t/calls.txt | sort -u | wc -l);;
                rockfish)  n=$(cut -f1 $W/$d/$t/calls.tsv | sort -u | wc -l);;
                deepmod2)  n=$(cat $W/$d/$t/calls/*per_read* 2>/dev/null | cut -f1 | sort -u | wc -l);;
                *) n=0;; esac; s="$n ($((100 * n / (nb > 0 ? nb : 1)))%)"; fi
            printf " %-13s" "$s"
        done; echo
    done ;;

status)
    squeue -u $USER -o "%.9i %.28j %.3t %.9M %R" | grep -E "mtx_|JOBID"; echo
    printf "%-12s %-5s" sample prep; for t in $TOOLS; do printf " %-13s" $t; done; echo
    for d in $DS; do
        printf "%-12s %-5s" $d "$([ -s $W/$d/status/prep.done ] && echo ok || echo -)"
        for t in $TOOLS; do
            if [ -s $W/$d/$t/sites.std.tsv ]; then s=ok; elif [ -n "$(dep infer:$t:$d)" ]; then s=queued; elif [ -s $W/$d/$t/FAILED ]; then s=FAILED; else s=-; fi
            printf " %-13s" $s
        done; echo
    done
    [ -s $W/matrix_grid.tsv ] && { echo; echo "== grid (AUROC; mean P for UniMeth, call frequency otherwise)"; column -t -s $'\t' $W/matrix_grid.tsv; }
    for d in $DS; do for t in $TOOLS; do [ -s $W/$d/$t/FAILED ] && [ -z "$(dep infer:$t:$d)" ] && echo "FAILED $t on $d: $(cat $W/$d/$t/FAILED)"; done; done; true ;;

# ---------------------------------------------------------------- job bodies (run inside sbatch)
_prep)
    d=$DSID; D=$W/$d; mkdir -p $D/status; : > $D/status/prep.txt; eval "$ENVACT"; set -x
    RU=${REUSE[$d]}; src=${BAM[$d]}; n=${NREADS[$d]}; [ -s $D/moves.bam ] && src=$D/moves.bam
    if [ "$RU" != "-" ] && [ -s $EUK/$RU/ref.fa.fai ]; then ln -sf $EUK/$RU/ref.fa $D/ref.fa; ln -sf $EUK/$RU/ref.fa.fai $D/ref.fa.fai
    elif [ ! -s $D/ref.fa.fai ]; then cp -L ${REF[$d]} $D/ref.fa && $SAM faidx $D/ref.fa || exit 1; fi
    [ -s $D/pod5_ids.txt ] || python $HERE/pod5_ids.py ${POD5[$d]} 2>> $D/status/prep.txt | sort -u > $D/pod5_ids.txt || exit 1
    if [ "$RU" != "-" ] && [ -s $EUK/$RU/sub.bam ]; then ln -sf $EUK/$RU/sub.bam $D/sub.bam; ln -sf $EUK/$RU/sub.bam.bai $D/sub.bam.bai; echo "sub.bam reused from $EUK/$RU" >> $D/status/prep.txt
    elif [ ! -s $D/sub.bam ]; then
        if [ $(stat -Lc %s $src) -gt 20000000000 ]; then   # a whole-run BAM: keep only the pod5's reads first, or the first N alignments hold almost no read with signal (rice: 1,093 of 15,000)
            [ -s $D/reads.inpod5.bam ] || $SAM view -@ 8 -b -N $D/pod5_ids.txt -o $D/reads.inpod5.bam $src || exit 1
            echo "source BAM cut to the pod5's reads: $($SAM view -c $D/reads.inpod5.bam) records" >> $D/status/prep.txt; src=$D/reads.inpod5.bam
        fi
        if $SAM view -H $src | grep -q 'SO:coordinate'; then ln -sf $src $D/reads.sorted.bam; { [ -s $src.bai ] && ln -sf $src.bai $D/reads.sorted.bam.bai; } || $SAM index $D/reads.sorted.bam
        else $SAM sort -@ 8 -m 2G -T $D/tmp_sort -o $D/reads.sorted.bam $src && $SAM index $D/reads.sorted.bam || exit 1; fi
        RF=${RFROM[$d]}
        if [ "$RF" != "-" ]; then      # same region as the other sample's subset (two-sample rows need matching coverage)
            [ -s $W/$RF/sub.bam ] || { echo "region_from $RF has no sub.bam yet: prep $RF first" >> $D/status/prep.txt; exit 1; }
            REG=$($SAM view $W/$RF/sub.bam | awk 'NR==1{c=$3; s=$4} {if($3==c){e=$4+length($10)}} END{print c":"s"-"e}')
            $SAM view -b -o $D/sub.bam $D/reads.sorted.bam $REG && $SAM index $D/sub.bam && echo "subset = region $REG of $RF" >> $D/status/prep.txt || exit 1
        elif [ $n -eq 0 ]; then ln -sf $D/reads.sorted.bam $D/sub.bam; ln -sf $D/reads.sorted.bam.bai $D/sub.bam.bai
        else $SAM view -h $D/reads.sorted.bam | awk -v n=$n '/^@/ {print; next} c<n {print; c++}' | $SAM view -b -o $D/sub.bam - && $SAM index $D/sub.bam || exit 1; fi
    fi
    # keep only reads whose signal is in the pod5 (the collection's oligo pod5s are per-replicate subsets of a BAM that
    # covers the whole run; every tool skips or stalls on reads it cannot find, so all five must see the same reads)
    if [ ! -s $D/status/pod5_filtered ]; then
        n0=$($SAM view -c $D/sub.bam)
        $SAM view -b -N $D/pod5_ids.txt -o $D/sub.inpod5.bam $D/sub.bam && $SAM index $D/sub.inpod5.bam || exit 1
        rm -f $D/sub.bam $D/sub.bam.bai; mv $D/sub.inpod5.bam $D/sub.bam; mv $D/sub.inpod5.bam.bai $D/sub.bam.bai
        echo "reads with signal in the pod5: $($SAM view -c $D/sub.bam) of $n0 alignments kept" | tee -a $D/status/prep.txt > $D/status/pod5_filtered
    fi
    echo "sub.bam: $($SAM view -c $D/sub.bam) alignments, mv tags in first 200: $($SAM view $D/sub.bam | head -200 | grep -c 'mv:B'), span: $($SAM view $D/sub.bam | awk 'NR==1{c=$3; s=$4} {e=$4} END{print c":"s"-"e}')" >> $D/status/prep.txt
    echo "prep finished $(date)" >> $D/status/prep.txt; cp $D/status/prep.txt $D/status/prep.done ;;

_basecall)
    d=$DSID; D=$W/$d; mkdir -p $D/status; set -x
    [ -s $D/ref.fa.fai ] || { cp -L ${REF[$d]} $D/ref.fa && $SAM faidx $D/ref.fa || exit 1; }
    # NOTRIM=1: Dorado's adapter/primer trimming leaves the move table shorter than the trimmed sequence on short reads
    # (oligos: 121 move bases for 142 bases), and UniMeth silently drops every such read; untrimmed reads keep them aligned
    $DORADO basecaller $DMODEL ${POD5[$d]} --emit-moves ${NOTRIM:+--no-trim} --reference $D/ref.fa > $D/moves.unsorted.bam \
        && $SAM sort -@ 8 -m 2G -T $D/tmp_bc -o $D/moves.bam $D/moves.unsorted.bam && $SAM index $D/moves.bam && rm -f $D/moves.unsorted.bam \
        && rm -f $D/reads.sorted.bam $D/reads.sorted.bam.bai $D/sub.bam $D/sub.bam.bai $D/status/prep.done \
        && echo "basecalled with move tables: $($SAM flagstat $D/moves.bam | grep -m1 'mapped (')" >> $D/status/prep.txt || { rm -f $D/moves.bam; echo 'basecall FAILED' >> $D/status/prep.txt; exit 1; } ;;

_infer)
    d=$DSID; t=$TOOL; D=$W/$d; T=$D/$t; mkdir -p $T; rm -f $T/FAILED; set -x
    [ -s $D/sub.bam ] || { echo "prep outputs missing" > $T/FAILED; exit 1; }
    RU=${REUSE[$d]}
    case $t in
    unimeth_5mC|unimeth_6mA|unimeth_5hmU|unimeth_5hmC|unimeth_4mC)
        eval "$ENVACT"
        case $t in
            unimeth_5mC)  M=$M_5mC;  FLAGS="--cpg 1 --chg 1 --chh 1"; TYPES='[CpG],[CHG],[CHH]'; CMD="unimeth infer" ;;
            unimeth_6mA)  M=$M_6mA;  FLAGS="--cpg 0 --chg 0 --chh 0 --m6A 1"; TYPES='[m6A]'; CMD="unimeth infer" ;;
            unimeth_5hmU) M=$M_5hmU; FLAGS="--cpg 0 --chg 0 --chh 0 --hmU 1"; TYPES='[5hmU]'; CMD="python -m unimeth.inference"; export PYTHONPATH=$UM_PATCHED ;;
            unimeth_5hmC|unimeth_4mC) M=$( [ $t = unimeth_5hmC ] && echo $M_5hmC || echo $M_4mC ); FLAGS="--cpg 1 --chg 1 --chh 1"; TYPES='[CpG],[CHG],[CHH]'; CMD="python -m unimeth.inference"; export PYTHONPATH=$UM_PATCHED ;;   # fine-tuned all-context C models
        esac
        [ -s $M ] || { echo "model missing: $M" > $T/FAILED; exit 1; }
        [ -s $T/calls.txt ] || { rm -f $T/part.txt; PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True $CMD --pod5 ${POD5[$d]} --bam $D/sub.bam --model $M --pore_type R10.4.1 --frequency 5khz $FLAGS --batch_size $UM_BATCH --output_format tsv --out $T/part.txt --num_workers 8 --signal_index $T/signal-index.sqlite > $T/infer.log 2>&1 \
            && mv $T/part.txt $T/calls.txt || { echo "inference failed: $(grep -iE 'error' $T/infer.log | tail -1)" > $T/FAILED; exit 1; }; }
        python $HERE/sites_std.py --tool unimeth --types "$TYPES" $T/calls.txt $T/sites.std.tsv > $T/std.log 2>&1 || { echo "sites_std failed: $(tail -1 $T/std.log)" > $T/FAILED; exit 1; } ;;
    deepmod2)
        if [ "$RU" != "-" ] && ls $EUK/$RU/deepmod2/calls/*per_site* > /dev/null 2>&1 && [ ! -d $T/calls ]; then ln -sfn $EUK/$RU/deepmod2/calls $T/calls; fi
        if big_ref $d; then TH=4; else TH=12; fi
        ls $T/calls/*per_site* > /dev/null 2>&1 || { rm -rf $T/calls.part; $DM2_ENV/bin/python $DM2_SRC/deepmod2 detect --bam $D/sub.bam --input ${POD5[$d]} --file_type pod5 --model $DM2_MODEL --seq_type dna --ref $D/ref.fa --threads $TH --output $T/calls.part > $T/detect.log 2>&1 \
            && mv $T/calls.part $T/calls || { echo "deepmod2 detect failed: $(grep -iE 'error|exception' $T/detect.log | tail -1)" > $T/FAILED; exit 1; }; }
        eval "$ENVACT"; python $HERE/sites_std.py --tool deepmod2 $T/calls $T/sites.std.tsv > $T/std.log 2>&1 || { echo "sites_std failed: $(tail -1 $T/std.log)" > $T/FAILED; exit 1; } ;;
    rockfish)
        source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $RF_ENV
        [ -s $RF_MODEL ] || { echo "model missing: $RF_MODEL" > $T/FAILED; exit 1; }
        # Rockfish walks every read of the pod5 input (~20 reads/s on CPU), so the pod5 is first cut down to the reads of sub.bam
        if [ ! -s $T/sub.pod5 ]; then
            $SAM view $D/sub.bam | cut -f1 | sort -u > $T/ids.txt
            if [ -d ${POD5[$d]} ]; then PIN="-r ${POD5[$d]}"; else PIN="${POD5[$d]}"; fi
            PFOPT=""; pod5 filter --help 2>/dev/null | grep -q -- "--missing-ok" && PFOPT="--missing-ok"      # flag names differ between pod5 versions
            pod5 filter --ids $T/ids.txt --output $T/sub.pod5 $PFOPT -t 8 $PIN > $T/pod5_filter.log 2>&1 || { rm -f $T/sub.pod5; echo "pod5 filter failed: $(tail -1 $T/pod5_filter.log)" > $T/FAILED; exit 1; }
        fi
        [ -s $T/calls.tsv ] || { rm -f $T/part.tsv; rockfish inference -i $T/sub.pod5 --bam_path $D/sub.bam --model_path $RF_MODEL -d 0 -t 8 -b 512 -o $T/part.tsv > $T/infer.log 2>&1 \
            && mv $T/part.tsv $T/calls.tsv || { echo "rockfish inference failed: $(grep -iE 'error' $T/infer.log | tail -1)" > $T/FAILED; exit 1; }; }
        eval "$ENVACT"; python $HERE/sites_std.py --tool rockfish --bam $D/sub.bam --orient $RF_ORIENT --minus-shift $RF_SHIFT $T/calls.tsv $T/sites.std.tsv > $T/std.log 2>&1 || { echo "sites_std failed: $(tail -1 $T/std.log)" > $T/FAILED; exit 1; } ;;
    *)  echo "unknown tool $t" > $T/FAILED; exit 1 ;;
    esac
    echo "$t on $d finished $(date): $(wc -l < $T/sites.std.tsv) site rows" ;;

_score)
    # CPGMERGE=row1,row2: variant scoring of CpG rows with both strands of a CpG summed before the floor; written to
    # matrix_{long,grid}_cpgmerge.tsv and status/score_cpgmerge.txt, the main results are untouched
    EXTRA=""; SFX=""; [ -n "${CPGMERGE:-}" ] && { EXTRA="--cpg-merge $CPGMERGE --only $CPGMERGE --suffix _cpgmerge"; SFX=_cpgmerge; }
    eval "$ENVACT"; python $HERE/score_matrix.py --work $W --datasets $HERE/datasets.tsv --rows $HERE/rows.tsv --tools "$TOOLS" --repo $REPO --samtools $SAM $EXTRA 2>&1 | tee $W/status/score$SFX.txt ;;

*)  echo "MODE must be check | basecall | prep | infer | score | audit | all | status"; exit 1 ;;
esac
