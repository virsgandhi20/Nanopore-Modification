#!/bin/bash
# Table 1 rows for datasets in the lab's data collection (plants, human): UniMeth 5mC and DeepMod2,
# scored against the collection's OWN ground truth (ground_truth/gt_{plus,minus}.bed = positives,
# cand_{plus,minus}.bed = candidates), the files RawMod itself was scored on, so the columns are
# comparable by construction. Positives = gt, negatives = candidates minus gt, coverage floor 10,
# strands collapsed per position (score_sites.py). One dataset per invocation.
#
#   DS=arabidopsis_col0_r10.4.1_ontbasemod_2024 MODE=prep     bash run_catalog_rows.sh   # CPU: sorted BAM, ref copy, first-N-reads subset, GT by context
#   DS=... MODE=unimeth  bash run_catalog_rows.sh   # GPU: all-context 5mC model on the subset; CpG and non-CpG rows scored
#   DS=... MODE=deepmod2 bash run_catalog_rows.sh   # GPU: DeepMod2 CpG model on the same subset
#   DS=... MODE=status   bash run_catalog_rows.sh
# The subset = the first NREADS (15,000) alignments of the coordinate-sorted BAM = one contiguous region
# at full depth, the same convention as the bacterial UniMeth rows.
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
CAT=/fs/cbcb-lab/storm/shared/data
ME=/fs/nexus-scratch/vgandhi
DS=${DS:?set DS=<dataset folder name under $CAT>}
SHORT=${SHORT:-${DS%%_*}}
W=${WBASE:-/fs/cbcb-lab/storm/vgandhi/euk}/$DS; mkdir -p $W/logs $W/status
# Defaults follow the collection's layout; each can be overridden for a dataset outside it (rice):
#   SRC_BAM=<move-table BAM> POD5=<pod5 dir> REF_SRC=<fasta> GT_POS="<positives bed(s)>" GT_CAND="<candidate bed(s)>"
SRC_BAM=${SRC_BAM:-$CAT/$DS/basecalled/reads.bam}; POD5=${POD5:-$CAT/$DS/pod5_files}; GT=$CAT/$DS/ground_truth
REF_SRC=${REF_SRC:-$CAT/$DS/ref.fa}; GT_POS=${GT_POS:-"$GT/gt_plus.bed $GT/gt_minus.bed"}; GT_CAND=${GT_CAND:-"$GT/cand_plus.bed $GT/cand_minus.bed"}
MINCOV=${MINCOV:-10}          # coverage floor for scoring; the human set is ~1x, so its row needs MINCOV=1 and more reads
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
ENVACT="source $HOME/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth"
MODEL_5mC=$ME/unimeth_models/checkpoints/unimeth_r10.4.1_5kHz_5mC.pt; UNIMETH_SRC=$ME/Unimeth
DM2_ENV=$ME/envs/deepmod2; DM2_SRC=$ME/deepmod2_bench/DeepMod2; DM2_MODEL=${DM2_MODEL:-bilstm_r10.4.1_5khz_v5.0}
DM2_THREADS=${DM2_THREADS:-12}; DM2_MEM=${DM2_MEM:-48G}     # human: DM2_THREADS=4 DM2_MEM=120G (each worker holds the 3 Gb reference; 12 workers OOM-killed at 48G)
NREADS=${NREADS:-15000}
SB="--account=scavenger --partition=scavenger --qos=scavenger --requeue"; GPU="--gres=gpu:rtxa5000:1"
CPU="--account=cbcb --partition=cbcb --qos=high"
MODE=${MODE:-status}

sub() {  # syntax-check the job body, then submit; prints the job id (empty on failure)
    local id a; for a in "$@"; do case "$a" in --wrap=*) bash -n <(printf '%s\n' "${a#--wrap=}") || { echo "NOT submitted: shell syntax error in the job body (above)" >&2; echo ""; return; };; esac; done
    id=$(sbatch --parsable $SB "$@") || { echo "sbatch failed: $*" >&2; echo ""; return; }; echo ${id%%;*}; }
dep() { local j; j=$(awk -v k=$1 -F'\t' '$1==k{print $2}' $W/jobs.tsv 2>/dev/null | tail -1); [ -n "$j" ] && [ -n "$(squeue -h -j $j 2>/dev/null)" ] && echo "--dependency=afterany:$j"; }

for p in $SRC_BAM $POD5 $REF_SRC $GT_POS $GT_CAND; do [ -e $p ] || { echo "missing: $p"; exit 1; }; done

case $MODE in
prep)
    : > $W/status/prep.txt
    J=$(sub $CPU --cpus-per-task=8 --mem=32G --time=08:00:00 --job-name=cat_prep_$SHORT --output=$W/logs/prep_%j.log --wrap="$ENVACT; set -uo pipefail; cd $W
      [ -s $W/ref.fa.fai ] || { cp -L $REF_SRC $W/ref.fa && $SAM faidx $W/ref.fa && echo 'reference copied and indexed' >> $W/status/prep.txt; }
      if [ ! -s $W/reads.sorted.bam ]; then
        if $SAM view -H $SRC_BAM | grep -q 'SO:coordinate'; then ln -sf $SRC_BAM $W/reads.sorted.bam; { [ -s $SRC_BAM.bai ] && ln -sf $SRC_BAM.bai $W/reads.sorted.bam.bai; } || $SAM index $W/reads.sorted.bam; echo 'source BAM already coordinate-sorted' >> $W/status/prep.txt
        else $SAM sort -@ 8 -m 2G -T $W/tmp_sort -o $W/reads.sorted.bam $SRC_BAM && $SAM index $W/reads.sorted.bam && echo 'BAM sorted and indexed' >> $W/status/prep.txt; fi
      fi
      if [ $NREADS -eq 0 ]; then ln -sf $W/reads.sorted.bam $W/sub.bam; ln -sf $W/reads.sorted.bam.bai $W/sub.bam.bai; else [ -s $W/sub.bam ] || { $SAM view -h $W/reads.sorted.bam | awk -v n=$NREADS '/^@/ {print; next} c<n {print; c++}' | $SAM view -b -o $W/sub.bam - && $SAM index $W/sub.bam; }; fi
      echo \"sub.bam: \$($SAM view -c $W/sub.bam) alignments, mv tags in first 200: \$($SAM view $W/sub.bam | head -200 | grep -c 'mv:B'), span: \$($SAM view $W/sub.bam | awk 'NR==1{c=\$3; s=\$4} {e=\$4} END{print c\":\"s\"-\"e}')\" >> $W/status/prep.txt
      python $REPO/scripts/ground_truth/split_by_context.py --ref $W/ref.fa --bed $GT_POS --out-cpg $W/gt_cpg.bed --out-noncpg $W/gt_noncpg.bed >> $W/status/prep.txt 2>&1
      python $REPO/scripts/ground_truth/split_by_context.py --ref $W/ref.fa --bed $GT_CAND --out-cpg $W/cand_cpg.bed --out-noncpg $W/cand_noncpg.bed >> $W/status/prep.txt 2>&1
      wc -l $W/gt_cpg.bed $W/gt_noncpg.bed $W/cand_cpg.bed $W/cand_noncpg.bed | sed 's/^/  /' >> $W/status/prep.txt
      echo \"prep finished \$(date)\" >> $W/status/prep.txt")
    echo "prep: job $J -> $W"; echo -e "prep\t$J" >> $W/jobs.tsv ;;

unimeth)
    : > $W/status/unimeth.txt
    J=$(sub $(dep prep) $GPU --cpus-per-task=8 --mem=48G --time=06:00:00 --job-name=cat_um_$SHORT --output=$W/logs/unimeth_%j.log --wrap="$ENVACT; set -uo pipefail; U=$W/unimeth; mkdir -p \$U
      [ -s $W/sub.bam ] && [ -s $W/cand_cpg.bed ] || { echo 'unimeth: prep outputs missing, run MODE=prep first' >> $W/status/unimeth.txt; exit 1; }
      [ -s \$U/calls.txt ] || unimeth-infer --pod5 $POD5 --bam $W/sub.bam --model $MODEL_5mC --pore_type R10.4.1 --frequency 4khz --cpg 1 --chg 1 --chh 1 --output_format tsv --out \$U/calls.txt --num_workers 8 --signal_index \$U/signal-index.sqlite > \$U/infer.log 2>&1 || { echo \"unimeth: inference FAILED: \$(grep -iE 'error' \$U/infer.log | tail -1)\" >> $W/status/unimeth.txt; exit 1; }
      [ -s \$U/sites.tsv ] || python $UNIMETH_SRC/scripts/call_modification_frequency.py -i \$U/calls.txt -o \$U/sites.tsv --sort
      for ctx in cpg noncpg; do
        [ -s $W/cand_\$ctx.bed ] || { echo \"no candidates for \$ctx, row skipped\" >> $W/status/unimeth.txt; continue; }
        python $REPO/scripts/benchmark/score_sites.py --calls \$U/sites.tsv --gt $W/gt_\$ctx.bed --candidates $W/cand_\$ctx.bed --min-cov $MINCOV --label \"UniMeth 5mC \$ctx $SHORT cov$MINCOV (call freq)\" --out $W/status/table.tsv >> $W/status/unimeth.txt 2>&1
        python $REPO/scripts/benchmark/score_sites.py --calls \$U/sites.tsv --gt $W/gt_\$ctx.bed --candidates $W/cand_\$ctx.bed --min-cov $MINCOV --num-col 5 --label \"UniMeth 5mC \$ctx $SHORT cov$MINCOV (mean P)\" --out $W/status/table.tsv >> $W/status/unimeth.txt 2>&1
      done
      echo \"unimeth finished \$(date)\" >> $W/status/unimeth.txt")
    echo "unimeth: job $J -> $W/unimeth"; echo -e "unimeth\t$J" >> $W/jobs.tsv ;;

deepmod2)
    : > $W/status/deepmod2.txt
    J=$(sub $(dep prep) $GPU --cpus-per-task=$DM2_THREADS --mem=$DM2_MEM --time=08:00:00 --job-name=cat_dm2_$SHORT --output=$W/logs/deepmod2_%j.log --wrap="set -uo pipefail; D=$W/deepmod2; mkdir -p \$D
      [ -s $W/sub.bam ] && [ -s $W/cand_cpg.bed ] || { echo 'deepmod2: prep outputs missing, run MODE=prep first' >> $W/status/deepmod2.txt; exit 1; }
      ls \$D/calls/*per_site* > /dev/null 2>&1 || $DM2_ENV/bin/python $DM2_SRC/deepmod2 detect --bam $W/sub.bam --input $POD5 --file_type pod5 --model $DM2_MODEL --seq_type dna --ref $W/ref.fa --threads $DM2_THREADS --output \$D/calls > \$D/detect.log 2>&1 || { echo \"deepmod2 detect FAILED: \$(grep -iE 'error|exception' \$D/detect.log | tail -1)\" >> $W/status/deepmod2.txt; exit 1; }
      $ENVACT
      python $REPO/scripts/benchmark/deepmod2_sites.py \$D/calls \$D/sites.tsv >> $W/status/deepmod2.txt 2>&1 || exit 1
      python $REPO/scripts/benchmark/score_sites.py --calls \$D/sites.tsv --gt $W/gt_cpg.bed --candidates $W/cand_cpg.bed --min-cov $MINCOV --chrom-col 0 --pos-col 1 --cov-col 2 --freq-col 3 --label \"DeepMod2 CpG $SHORT cov$MINCOV\" --out $W/status/table.tsv >> $W/status/deepmod2.txt 2>&1
      echo \"deepmod2 finished \$(date)\" >> $W/status/deepmod2.txt")
    echo "deepmod2: job $J -> $W/deepmod2"; echo -e "deepmod2\t$J" >> $W/jobs.tsv ;;

status)
    squeue -u $USER -o "%.9i %.16j %.3t %.9M %R" | grep -E "cat_|JOBID"; echo
    for f in prep unimeth deepmod2; do [ -s $W/status/$f.txt ] && { echo "== $f"; cat $W/status/$f.txt; echo; }; done
    [ -s $W/status/table.tsv ] && { echo "== table rows"; cat $W/status/table.tsv; } ;;
*)  echo "MODE must be prep | unimeth | deepmod2 | status"; exit 1 ;;
esac
