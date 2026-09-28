#!/bin/bash
# Basecall a pod5 collection with move tables, aligned to a reference, as N parallel GPU jobs plus one merge job.
# Made for the two datasets the table still needs at floor 10: the full public HG002 flowcell and Ernest's rice files.
#
#   NAME=hg002_10x POD5=<dir with .pod5> REF=<fasta> OUT=<dir> [NJOBS=24] [PART=cbcb] bash basecall_moves.sh
#   (Bhargav, Sep 27: up to 100 concurrent jobs are fine, so use many shards)
#   MODE=status NAME=... OUT=... bash basecall_moves.sh
#   SHARDS="16 17 18" MERGE=0 ... bash basecall_moves.sh     submit only those shards (e.g. the pending half on the other partition)
#   MODE=merge NAME=... OUT=... bash basecall_moves.sh        merge job that waits for every queued shard of this NAME
# Result: $OUT/$NAME.moves.bam (+ .bai), coordinate-sorted, Dorado 1.4 sup v5 with --emit-moves, ready for datasets.tsv.
# Each shard is written as .part and renamed on success, so a preempted shard restarts cleanly; finished shards are skipped.
set -uo pipefail
NAME=${NAME:?}; OUT=${OUT:?}; NJOBS=${NJOBS:-24}; MODE=${MODE:-run}; POD5=${POD5:-}; REF=${REF:-}
[ $MODE = run ] && { : ${POD5:?set POD5=<dir with .pod5>}; : ${REF:?set REF=<indexed fasta>}; }
ME=/fs/nexus-scratch/vgandhi
DORADO=/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-1.4.0-linux-x64/bin/dorado; DMODEL=$ME/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
SB="--account=scavenger --partition=scavenger --qos=scavenger --requeue"; GPU="--gres=gpu:rtxa5000:1"; CPU="--account=cbcb --partition=cbcb --qos=high"
[ "${PART:-}" = cbcb ] && { SB="--account=cbcb --partition=cbcb --qos=high"; GPU="--gres=gpu:1"; }   # PART=cbcb: lab partition instead of scavenger
mkdir -p $OUT/shards $OUT/logs $OUT/lists
sub() { local id; id=$(sbatch --parsable $SB "$@") || { echo "sbatch failed" >&2; echo ""; return; }; echo ${id%%;*}; }

case $MODE in
run)
    [ -s $REF.fai ] || { echo "reference index missing: $REF.fai (copy the reference somewhere writable and run samtools faidx)"; exit 1; }
    if [ -s $OUT/lists/njobs.txt ] && [ "$(cat $OUT/lists/njobs.txt)" != "$NJOBS" ]; then echo "this OUT was sharded with NJOBS=$(cat $OUT/lists/njobs.txt); pass the same NJOBS (a different split would overlap and skip files)"; exit 1; fi
    echo $NJOBS > $OUT/lists/njobs.txt
    find -L $POD5 -name '*.pod5' | sort > $OUT/lists/all.txt; N=$(wc -l < $OUT/lists/all.txt)
    [ $N -gt 0 ] || { echo "no pod5 files under $POD5"; exit 1; }
    echo "$N pod5 files -> $NJOBS shards"; JOBS=""
    for i in ${SHARDS:-$(seq 0 $((NJOBS - 1)))}; do
        awk -v i=$i -v n=$NJOBS 'NR % n == i' $OUT/lists/all.txt > $OUT/lists/shard_$i.txt
        [ -s $OUT/lists/shard_$i.txt ] || continue
        [ -s $OUT/shards/shard_$i.bam ] && { echo "shard $i: done already"; continue; }
        J=$(sub $GPU --cpus-per-task=8 --mem=48G --time=12:00:00 --job-name=bc_${NAME}_$i --output=$OUT/logs/shard_${i}_%j.log \
              --wrap="set -uo pipefail; mkdir -p $OUT/shards/in_$i; rm -f $OUT/shards/in_$i/*; while read f; do ln -sf \$f $OUT/shards/in_$i/; done < $OUT/lists/shard_$i.txt
                      rm -f $OUT/shards/shard_$i.part.bam
                      $DORADO basecaller $DMODEL $OUT/shards/in_$i --emit-moves --reference $REF > $OUT/shards/shard_$i.part.bam && mv $OUT/shards/shard_$i.part.bam $OUT/shards/shard_$i.bam && echo \"shard $i OK: \$($SAM view -c $OUT/shards/shard_$i.bam) records\" || { rm -f $OUT/shards/shard_$i.part.bam; echo \"shard $i FAILED\"; exit 1; }")
        [ -n "$J" ] && { echo "shard $i: job $J ($(wc -l < $OUT/lists/shard_$i.txt) files)"; JOBS="$JOBS:$J"; }
    done
    [ "${MERGE:-1}" = 0 ] && { echo "no merge submitted (MERGE=0): run MODE=merge once every shard is queued"; exit 0; }
    J=$(sub ${JOBS:+--dependency=afterok$JOBS} $CPU --cpus-per-task=8 --mem=48G --time=12:00:00 --job-name=bc_${NAME}_merge --output=$OUT/logs/merge_%j.log \
          --wrap="set -uo pipefail; ls $OUT/shards/shard_*.bam > $OUT/lists/done.txt; [ \$(wc -l < $OUT/lists/done.txt) -gt 0 ] || exit 1
                  $SAM cat -o $OUT/$NAME.unsorted.bam \$(cat $OUT/lists/done.txt) && $SAM sort -@ 8 -m 3G -T $OUT/tmp_sort -o $OUT/$NAME.moves.bam $OUT/$NAME.unsorted.bam && $SAM index $OUT/$NAME.moves.bam && rm -f $OUT/$NAME.unsorted.bam
                  echo \"merged: \$($SAM flagstat $OUT/$NAME.moves.bam | grep -m2 -E 'primary mapped|in total')\" > $OUT/status.txt; echo \"$NAME.moves.bam ready \$(date)\" >> $OUT/status.txt")
    echo "merge: job $J (after the shards); result: $OUT/$NAME.moves.bam" ;;
merge)   # waits for every shard job of this NAME still in the queue (any partition), then merges
    JOBS=$(squeue -h -u $USER -o "%i %j" | awk -v n="bc_${NAME}_" '$2 ~ "^"n"[0-9]+$"{printf ":%s", $1}')
    J=$(sub ${JOBS:+--dependency=afterok$JOBS} $CPU --cpus-per-task=8 --mem=48G --time=12:00:00 --job-name=bc_${NAME}_merge --output=$OUT/logs/merge_%j.log \
          --wrap="set -uo pipefail; ls $OUT/shards/shard_*.bam > $OUT/lists/done.txt; [ \$(wc -l < $OUT/lists/done.txt) -gt 0 ] || exit 1
                  $SAM cat -o $OUT/$NAME.unsorted.bam \$(cat $OUT/lists/done.txt) && $SAM sort -@ 8 -m 3G -T $OUT/tmp_sort -o $OUT/$NAME.moves.bam $OUT/$NAME.unsorted.bam && $SAM index $OUT/$NAME.moves.bam && rm -f $OUT/$NAME.unsorted.bam
                  echo \"merged: \$($SAM flagstat $OUT/$NAME.moves.bam | grep -m2 -E 'primary mapped|in total')\" > $OUT/status.txt; echo \"$NAME.moves.bam ready \$(date)\" >> $OUT/status.txt")
    echo "merge: job $J (after shards${JOBS:-: none queued})" ;;

status)
    squeue -u $USER -o "%.9i %.22j %.3t %.9M %R" | grep -E "bc_${NAME}|JOBID"
    echo "shards done: $(ls $OUT/shards/shard_*.bam 2>/dev/null | grep -vc part) of $(ls $OUT/lists/shard_*.txt 2>/dev/null | wc -l)"; cat $OUT/status.txt 2>/dev/null
    grep -h "FAILED\|rror" $OUT/logs/*.log 2>/dev/null | tail -5 ;;
*) echo "MODE=run|merge|status"; exit 1 ;;
esac
