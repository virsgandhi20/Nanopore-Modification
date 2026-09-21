#!/bin/bash
# Overnight per-read typing runs. One entry point, four modes (login node):
#
#   MODE=disk       what is filling scratch; gzip UniMeth's per-read call files
#                   (lossless, frees ~25 GB); never deletes anything
#   MODE=preflight  build site lists, run the extractor on 150 real reads of each
#                   existing BAM, run the trainer for one epoch on that for every
#                   kind of experiment, validate the sbatch templates. No queue.
#   MODE=submit     snapshot the code, submit basecall -> extract -> train jobs
#   MODE=status     what finished, with the headline numbers (use in the morning)
#
# Design for "must not fail overnight": scripts are copied to $RUN/code at
# submit time (a later git pull cannot change a running job); train jobs depend
# on their inputs with afterany and tolerate a missing input; every job writes a
# status line; nothing uses `set -e`.
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
ME=/fs/nexus-scratch/vgandhi
RUN=${RUN:-$ME/typing_overnight}
UB=$ME/unimeth_bench
BENCH=/fs/cbcb-lab/storm/bds062/data/benchmark
ONT=/fs/nexus-scratch/bds062/data/ont-os
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
DORADO=/fs/cbcb-lab/storm/shared/rawhash2/basecallers/dorado-1.4.0-linux-x64/bin/dorado
DMODEL=$ME/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.0.0
ENVACT="source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth"
MODE=${MODE:-status}
mkdir -p $RUN/{bam,sites,npz,runs,logs,code} 2>/dev/null
N=$RUN/npz; S=$RUN/sites
# SBASE can be overridden to send a second, independent copy of the whole run to
# another queue (use a different RUN so the two copies never share a file), e.g.
#   RUN=$ME/typing_overnight_scav SBASE="--account=scavenger --partition=scavenger --qos=scavenger --requeue"
SBASE=${SBASE:-"--account=cbcb --partition=cbcb --qos=high --exclude=cbcb25"}

# ------------------------------------------------------------------ samples
# name | bam | pod5 | ref | sites (group=bed,...) | label map | mod-base for unstranded beds | reads per site
samples() { cat <<EOF
ecoli_wt|$UB/bam/Ecoli_WT_5kHz.moves.bam|$BENCH/bacteria/Ecoli_WT_5kHz/pod5|$UB/ref/ecoli.fa|dam=$S/ecoli/dam.bed,dcm=$S/ecoli/dcm.bed,cpg=$S/ecoli/cpg.bed,bgA=$S/ecoli/bgA.bed,bgC=$S/ecoli/bgC.bed|dam=6mA,dcm=5mC,cpg=none,bgA=none,bgC=none||12
ecoli_dm|$RUN/bam/Ecoli_DM_5kHz.moves.bam|$BENCH/bacteria/Ecoli_DM_5kHz/pod5|$UB/ref/ecoli.fa|dam=$S/ecoli/dam.bed,dcm=$S/ecoli/dcm.bed,cpg=$S/ecoli/cpg.bed,bgA=$S/ecoli/bgA.bed,bgC=$S/ecoli/bgC.bed|dam=none,dcm=none,cpg=none,bgA=none,bgC=none||12
ecoli_msssi|$UB/bam/Ecoli_DM_MSssI_5kHz.moves.bam|$BENCH/bacteria/Ecoli_DM_MSssI_5kHz/pod5|$UB/ref/ecoli.fa|cpg=$S/ecoli/cpg.bed|cpg=5mC||12
hp_wt|$UB/bam/HP26695_WT_5kHz.moves.bam|$BENCH/bacteria/HP26695_WT_5kHz/pod5|$UB/ref/hpylori_26695.fa|hp6mA=$ME/hp_labels/gt_6mA.bed,hp5mC=$ME/hp_labels/gt_5mC.bed,hp4mC=$ME/hp_labels/gt_4mC.bed,bgA=$S/hp/bgA.bed,bgC=$S/hp/bgC.bed|hp6mA=6mA,hp5mC=5mC,hp4mC=4mC,bgA=none,bgC=none||12
hp_wga|$RUN/bam/HP26695_WGA_5kHz.moves.bam|$BENCH/bacteria/HP26695_WGA_5kHz/pod5|$UB/ref/hpylori_26695.fa|hp6mA=$ME/hp_labels/gt_6mA.bed,hp5mC=$ME/hp_labels/gt_5mC.bed,hp4mC=$ME/hp_labels/gt_4mC.bed,bgA=$S/hp/bgA.bed,bgC=$S/hp/bgC.bed|hp6mA=none,hp5mC=none,hp4mC=none,bgA=none,bgC=none||12
anabaena|$UB/bam/Anabaena_WT_5kHz.moves.bam|$BENCH/bacteria/Anabaena_WT_5kHz/pod5|$UB/ref/anabaena_sp_PCC7120_ATCC27893.fa|dam=$S/anabaena/dam.bed,bgA=$S/anabaena/bgA.bed|dam=6mA,bgA=none||12
tdent|$UB/bam/Tdenticola_WT_5kHz.moves.bam|$BENCH/bacteria/Tdenticola_WT_5kHz/pod5|$UB/ref/treponema_denticola_ATCC35405.fa|d6mA=$S/tdent/pos.bed,dneg=$S/tdent/neg.bed|d6mA=6mA,dneg=none||12
j99|$UB/bam/HPJ99_WT_5kHz.moves.bam|$BENCH/bacteria/HPJ99_WT_5kHz/pod5|$UB/ref/hpylori_J99_ATCC700824.fa|d6mA=$S/j99/pos.bed,dneg=$S/j99/neg.bed|d6mA=6mA,dneg=none||12
EOF
for R in rep1 rep2; do
echo "syn_control_$R|$RUN/bam/syn_control_$R.moves.bam|$ONT/subset/control_$R.pod5|$ONT/references/all_5mers.fa|synC=$ONT/references/all_5mers_C_sites.bed,synA=$ONT/references/all_5mers_A_sites.bed|synC=none,synA=none||400 --min-mapq ${SYN_MAPQ:-5}"
echo "syn_5mC_$R|$RUN/bam/syn_5mC_$R.moves.bam|$ONT/subset/5mC_$R.pod5|$ONT/references/all_5mers.fa|synC=$ONT/references/all_5mers_5mC_sites.bed|synC=5mC||400 --min-mapq ${SYN_MAPQ:-5}"
echo "syn_5hmC_$R|$RUN/bam/syn_5hmC_$R.moves.bam|$ONT/subset/5hmC_$R.pod5|$ONT/references/all_5mers.fa|synC=$ONT/references/all_5mers_5hmC_sites.bed|synC=5hmC||400 --min-mapq ${SYN_MAPQ:-5}"
echo "syn_6mA_$R|$RUN/bam/syn_6mA_$R.moves.bam|$ONT/subset/6mA_$R.pod5|$ONT/references/all_5mers.fa|synA=$ONT/references/all_5mers_6mA_sites.bed|synA=6mA||400 --min-mapq ${SYN_MAPQ:-5}"
done; }
extract_cmd() {  # sample-row  out-npz  extra-args
    IFS='|' read -r NAME BAM POD REF SITES LMAP MBASE RPS <<< "$1"
    echo "python $RUN/code/extract_read_windows.py --sample $NAME --bam $BAM --pod5 $POD --ref $REF --sites $SITES --label-map $LMAP ${MBASE:+--mod-base $MBASE} --max-reads-per-site $RPS --out $2 ${3:-}"
}

# ------------------------------------------------------------------ experiments
ECO="$N/ecoli_wt.npz,$N/ecoli_dm.npz,$N/ecoli_msssi.npz"; HP="$N/hp_wt.npz,$N/hp_wga.npz"
SYN=$(for R in rep1 rep2; do for C in control 5mC 5hmC 6mA; do printf "$N/syn_${C}_$R.npz,"; done; done | sed 's/,$//')
NOVEL="anabaena=$N/anabaena.npz,tdent=$N/tdent.npz,j99=$N/j99.npz"
# name | train | classes | inputs | extra trainer args | note
experiments() { cat <<EOF
E01_syn_type4_all|$SYN|none,5mC,5hmC,6mA|sig,dwell,seq||ONT oligos, exact per-read labels, test = unseen 5-mer contexts
E02_syn_type4_signal|$SYN|none,5mC,5hmC,6mA|sig,dwell||same without sequence input
E03_syn_type4_seqonly|$SYN|none,5mC,5hmC,6mA|seq||control: sequence alone cannot tell none/5mC/5hmC apart (same sites)
E04_syn_5mC_vs_5hmC|$SYN|none,5mC,5hmC|sig,dwell,seq|--groups synC|the hard pair, cytosine sites only
E05_syn_cross_replicate|$SYN|none,5mC,5hmC,6mA|sig,dwell,seq|--train-samples rep1|train on flow cell 1, test on flow cell 2
E06_syn_detect_binary|$SYN|none,mod|sig,dwell,seq|--merge mod=5mC+5hmC+6mA|per-read detection, chemistry-agnostic
E07_ecoli_type3_all|$ECO|none,6mA,5mC|sig,dwell,seq||E. coli with matched dam-/dcm- controls at the same sites
E08_ecoli_type3_signal|$ECO|none,6mA,5mC|sig,dwell||same without sequence input
E09_ecoli_type3_seqonly|$ECO|none,6mA,5mC|seq||control: how much of bacterial 'typing' is motif recognition
E10_bact_type4_all|$ECO,$HP|none,6mA,5mC,4mC|sig,dwell,seq|--eval $NOVEL|adds H. pylori (4mC, many motifs); scored on unseen organisms too
E11_xfer_ecoli_to_hp_all|$ECO,$HP|none,6mA,5mC|sig,dwell,seq|--train-samples ecoli|cross-organism, cross-motif; 4mC arrives as an unknown class
E12_xfer_ecoli_to_hp_signal|$ECO,$HP|none,6mA,5mC|sig,dwell|--train-samples ecoli|same without sequence input
E13_xfer_syn_to_bacteria|$SYN|none,5mC,5hmC,6mA|sig,dwell,seq|--eval ecoli=${ECO//,/+},hp=${HP//,/+},$NOVEL|does chemistry learned on oligos transfer to genomes
E14_xfer_bacteria_to_syn|$ECO,$HP|none,6mA,5mC,4mC|sig,dwell,seq|--eval syn=${SYN//,/+}|reverse direction; 5hmC arrives as an unknown class
E15_openset_4mC|$ECO,$HP|none,6mA,5mC|sig,dwell,seq||4mC never trained: is it flagged, and what is it called
E16_openset_5hmC|$SYN|none,5mC,6mA|sig,dwell,seq||5hmC never trained: is it flagged or silently called 5mC
E17_pooled_all5|$SYN,$ECO,$HP|none,5mC,5hmC,6mA,4mC|sig,dwell,seq|--eval $NOVEL|everything, five classes
E18_pooled_detect_binary|$SYN,$ECO,$HP|none,mod|sig,dwell,seq|--merge mod=5mC+5hmC+6mA+4mC --eval $NOVEL|per-read detection across all data (RawMod's task, per read)
EOF
}
train_cmd() {  # experiment-row  out-dir  extra
    IFS='|' read -r NAME TRAIN CLASSES INPUTS EXTRA NOTE <<< "$1"
    echo "python $RUN/code/train_read_typing.py --name $NAME --out $2 --train $TRAIN --classes $CLASSES --inputs $INPUTS $EXTRA --note \"$NOTE\" ${3:-}"
}

build_sites() {
    eval "$ENVACT"
    [ -s $S/ecoli/dam.bed ] || python $REPO/scripts/typing/make_typing_sites.py --ref $UB/ref/ecoli.fa --out $S/ecoli
    [ -s $S/anabaena/dam.bed ] || python $REPO/scripts/typing/make_typing_sites.py --ref $UB/ref/anabaena_sp_PCC7120_ATCC27893.fa --out $S/anabaena --motifs dam
    mkdir -p $S/hp $S/tdent $S/j99
    # H. pylori background: positions Dorado calls essentially unmodified in the wild type
    [ -s $S/hp/bgA.bed ] || awk '$4=="a" && $10>=15 && $11<=3 {print $1"\t"$2"\t"$3"\tbgA\t0\t"$6}' $ME/hp_labels/HP26695_WT_5kHz_6mA.bed | shuf -n 30000 --random-source=<(yes) > $S/hp/bgA.bed
    [ -s $S/hp/bgC.bed ] || awk '$4=="m" && $10>=15 && $11<=3 {print $1"\t"$2"\t"$3"\tbgC\t0\t"$6}' $ME/hp_labels/HP26695_WT_5kHz_5mC.bed | shuf -n 30000 --random-source=<(yes) > $S/hp/bgC.bed
    # T. denticola / J99: caller-derived 6mA (the motif presets are not methylated in these samples)
    for X in tdent:Tdenticola_WT_5kHz j99:HPJ99_WT_5kHz; do
        D=$S/${X%%:*}; B=$UB/gtcheck/${X##*:}_6mA.bed
        [ -s $D/pos.bed ] || awk '$4=="a" && $10>=15 && $11>=80 {print $1"\t"$2"\t"$3"\td6mA\t0\t"$6}' $B | shuf -n 12000 --random-source=<(yes) > $D/pos.bed
        [ -s $D/neg.bed ] || awk '$4=="a" && $10>=15 && $11<=3  {print $1"\t"$2"\t"$3"\tdneg\t0\t"$6}' $B | shuf -n 12000 --random-source=<(yes) > $D/neg.bed
    done
    for f in $S/*/*.bed; do printf "  %-28s %8d sites\n" ${f#$S/} $(wc -l < $f); done
}

case $MODE in
# =========================================================================== disk
disk)
    echo "=== scratch usage"; df -h $ME | tail -1
    du -sh $ME/* 2>/dev/null | sort -rh | head -15
    echo; echo "=== gzip UniMeth per-read calls (lossless; the site tables stay as they are)"
    for f in $UB/*/calls.txt; do [ -s "$f" ] || continue
        echo "  $(du -h $f | cut -f1)  $f"; gzip -1 "$f" && echo "     -> $(du -h $f.gz | cut -f1)"; done
    echo; echo "=== largest regenerable intermediates (NOT deleted; your call)"
    find $ME/orca_feat $ME/orca_feat_spo1 $ME/orca_feat_bench $ME/spo1_typing -type f -size +1G \
        \( -name "*.blow5" -o -name "*eventalign*" -o -name "*.fastq" -o -name "*.pileup" -o -name "*.bam" -o -name "*.tsv" -o -name "*.txt" \) \
        -printf "%s\t%p\n" 2>/dev/null | sort -rn | head -25 | awk '{printf "  %6.1f GB  %s\n", $1/1e9, $2}'
    echo; df -h $ME | tail -1 ;;

# =========================================================================== preflight
preflight)
    FAILS=0; ok() { echo "  PASS  $1"; }; bad() { echo "  FAIL  $1"; FAILS=$((FAILS+1)); }
    FREE=$(df --output=avail -BG $ME | tail -1 | tr -dc 0-9)
    [ "${FREE:-0}" -ge 10 ] && ok "free space ${FREE}G (need 10G)" || bad "only ${FREE}G free on scratch; run MODE=disk first"
    cp $REPO/scripts/typing/{extract_read_windows.py,train_read_typing.py,make_typing_sites.py} $RUN/code/ && ok "code snapshot" || bad "code snapshot"
    echo "--- site lists"; build_sites
    for f in ecoli/dam ecoli/dcm ecoli/cpg ecoli/bgA ecoli/bgC anabaena/dam hp/bgA hp/bgC tdent/pos tdent/neg j99/pos j99/neg; do
        [ "$(wc -l < $S/$f.bed 2>/dev/null || echo 0)" -ge 500 ] && ok "sites $f" || bad "sites $f has under 500 rows"; done
    eval "$ENVACT"
    echo "--- extractor on 150 real reads per existing BAM"
    mkdir -p $RUN/pf/npz
    while IFS= read -r ROW; do NAME=${ROW%%|*}; BAM=$(echo "$ROW" | cut -d'|' -f2)
        [ -s "$BAM" ] || { echo "  ....  $NAME: BAM not basecalled yet (will be made by a GPU job)"; continue; }
        eval "$(extract_cmd "$ROW" $RUN/pf/npz/$NAME.npz "--max-bam-reads 150 --max-sites 4000")" > $RUN/pf/$NAME.extract.log 2>&1 \
            && ok "extract $NAME: $(grep -o 'wrote [0-9,]* windows' $RUN/pf/$NAME.extract.log) $(grep -o 'label counts.*' $RUN/pf/$NAME.extract.log)" \
            || { bad "extract $NAME"; tail -5 $RUN/pf/$NAME.extract.log | sed 's/^/        /'; }
    done < <(samples)
    echo "--- trainer, one epoch on the tiny real data, one run per kind of experiment"
    P=$RUN/pf/npz; PECO="$P/ecoli_wt.npz,$P/ecoli_msssi.npz"; PNOV="anabaena=$P/anabaena.npz,tdent=$P/tdent.npz,j99=$P/j99.npz"
    T="python $RUN/code/train_read_typing.py --epochs 1 --batch 128"
    pf() { local nm=$1; shift; eval "$T --name pf_$nm --out $RUN/pf/runs/$nm $*" > $RUN/pf/$nm.train.log 2>&1
           local st=$(cut -f2 $RUN/pf/runs/$nm/status.txt 2>/dev/null)
           case "$st" in OK*) ok "train $nm: $st  $(grep -m1 '\[test\]' $RUN/pf/$nm.train.log | cut -c1-110)";; *) bad "train $nm: ${st:-no status}"; tail -6 $RUN/pf/$nm.train.log | sed 's/^/        /';; esac; }
    pf type3_all      --train $PECO --classes none,6mA,5mC --inputs sig,dwell,seq
    pf type3_seqonly  --train $PECO --classes none,6mA,5mC --inputs seq
    pf type4_eval     --train $PECO,$P/hp_wt.npz --classes none,6mA,5mC,4mC --inputs sig,dwell --eval $PNOV
    pf xfer_openset   --train $PECO,$P/hp_wt.npz --classes none,6mA,5mC --train-samples ecoli
    pf binary_merge   --train $PECO,$P/hp_wt.npz --classes none,mod --merge mod=5mC+5hmC+6mA+4mC --eval $PNOV
    pf missing_inputs --train $PECO,$P/not_there.npz --classes none,6mA,5mC --eval gone=$P/nope.npz
    echo "--- basecalling inputs"
    [ -x $DORADO ] && ok "dorado 1.4.0" || bad "dorado binary"
    [ -d $DMODEL ] && ok "sup@v5.0.0 model" || bad "dorado model dir"
    H=$($DORADO basecaller --help 2>&1); for F in emit-moves reference; do [[ "$H" == *"--$F"* ]] && ok "dorado --$F" || bad "dorado --$F"; done
    for P5 in $BENCH/bacteria/Ecoli_DM_5kHz/pod5 $BENCH/bacteria/HP26695_WGA_5kHz/pod5 $ONT/subset/control_rep1.pod5 $ONT/subset/5hmC_rep2.pod5; do [ -r $P5 ] && ok "readable $P5" || bad "cannot read $P5"; done
    echo "  info  ONT oligo BAM from the lab (dorado-aligned): $($SAM flagstat $ONT/basecalls/control_rep1.bam 2>/dev/null | grep -m1 'mapped (' )"
    echo "  info  oligo read lengths (first 2000): $($SAM view $ONT/basecalls/control_rep1.bam | head -2000 | awk '{n++; s+=length($10)} END{printf "mean %.0f bp over %d reads", s/n, n}')"
    echo "--- SLURM templates (--test-only, nothing is queued)"
    sbatch --test-only $SBASE --gres=gpu:1 --cpus-per-task=8 --mem=48G --time=04:00:00 --wrap="true" > $RUN/pf/sb_gpu.txt 2>&1 && ok "GPU job template: $(tr '\n' ' ' < $RUN/pf/sb_gpu.txt | cut -c1-90)" || bad "GPU job template: $(cat $RUN/pf/sb_gpu.txt)"
    sbatch --test-only $SBASE --cpus-per-task=4 --mem=32G --time=04:00:00 --wrap="true" > $RUN/pf/sb_cpu.txt 2>&1 && ok "CPU job template: $(tr '\n' ' ' < $RUN/pf/sb_cpu.txt | cut -c1-90)" || bad "CPU job template: $(cat $RUN/pf/sb_cpu.txt)"
    echo; [ $FAILS -eq 0 ] && echo "PREFLIGHT: ALL PASS - safe to run MODE=submit" || echo "PREFLIGHT: $FAILS FAILURE(S) - do not submit, paste this output"
    exit $FAILS ;;

# =========================================================================== submit
submit)
    cp $REPO/scripts/typing/{extract_read_windows.py,train_read_typing.py,make_typing_sites.py} $RUN/code/
    build_sites > $RUN/logs/sites.txt 2>&1
    : > $RUN/jobs.tsv
    sub() { local id; id=$(sbatch --parsable $SBASE "$@") || { echo "sbatch failed: $*" >&2; echo ""; return; }; echo ${id%%;*}; }
    PRE="$ENVACT; set -uo pipefail; mkdir -p $RUN/status"   # activate before -u: conda hooks are not -u clean
    basecall() {  # tag  "name:pod5:ref ..."
        local tag=$1 body="$PRE; ok=1"
        for T3 in $2; do IFS=: read -r NM POD REF <<< "$T3"; local B=$RUN/bam/$NM.moves.bam
            body+="; if [ ! -s $B ]; then $DORADO basecaller $DMODEL $POD --emit-moves --reference $REF > $RUN/bam/$NM.unsorted.bam 2> $RUN/logs/dorado_$NM.err && $SAM sort -@ 8 -o $B $RUN/bam/$NM.unsorted.bam && $SAM index $B; rm -f $RUN/bam/$NM.unsorted.bam; fi"
            body+="; if [ -s $B ]; then echo \"$NM basecall OK: \$($SAM flagstat $B | grep -m1 'mapped (')\" >> $RUN/status/basecall_$tag.txt; else ok=0; echo \"$NM basecall FAILED: \$(tail -2 $RUN/logs/dorado_$NM.err | tr '\\n' ' ')\" >> $RUN/status/basecall_$tag.txt; fi"
        done
        sub --gres=gpu:1 --cpus-per-task=8 --mem=48G --time=03:00:00 --job-name=ty_bc_$tag --output=$RUN/logs/bc_${tag}_%j.log --wrap="$body; [ \$ok = 1 ]"
    }
    J_DM=$(basecall ecoli_dm "Ecoli_DM_5kHz:$BENCH/bacteria/Ecoli_DM_5kHz/pod5:$UB/ref/ecoli.fa")
    J_WGA=$(basecall hp_wga "HP26695_WGA_5kHz:$BENCH/bacteria/HP26695_WGA_5kHz/pod5:$UB/ref/hpylori_26695.fa")
    SYNL=""; for R in rep1 rep2; do for C in control 5mC 5hmC 6mA; do SYNL+="syn_${C}_$R:$ONT/subset/${C}_$R.pod5:$ONT/references/all_5mers.fa "; done; done
    J_SYN=$(basecall syn "$SYNL")
    echo -e "basecall\tecoli_dm\t$J_DM\nbasecall\thp_wga\t$J_WGA\nbasecall\tsyn_x8\t$J_SYN" >> $RUN/jobs.tsv
    EX_SYN=""; EX_BACT=""; EX_NOVEL=""
    while IFS= read -r ROW; do NAME=${ROW%%|*}; DEP=""
        case $NAME in ecoli_dm) DEP=$J_DM;; hp_wga) DEP=$J_WGA;; syn_*) DEP=$J_SYN;; esac
        CMD="$PRE; if [ -s $N/$NAME.npz ]; then echo '$NAME extract REUSED' > $RUN/status/extract_$NAME.txt; else $(extract_cmd "$ROW" $N/$NAME.npz) && echo \"$NAME extract OK: \$(du -h $N/$NAME.npz | cut -f1)\" > $RUN/status/extract_$NAME.txt || { echo '$NAME extract FAILED' > $RUN/status/extract_$NAME.txt; exit 1; }; fi"
        J=$(sub ${DEP:+--dependency=afterany:$DEP} --cpus-per-task=4 --mem=32G --time=02:00:00 --job-name=ty_ex_$NAME --output=$RUN/logs/ex_${NAME}_%j.log --wrap="$CMD")
        echo -e "extract\t$NAME\t$J" >> $RUN/jobs.tsv
        if [ -n "$J" ]; then case $NAME in syn_*) EX_SYN+=":$J";; anabaena|tdent|j99) EX_NOVEL+=":$J";; *) EX_BACT+=":$J";; esac; fi
    done < <(samples)
    while IFS= read -r ROW; do NAME=${ROW%%|*}
        DEPS=""                       # wait only for the inputs this experiment reads
        [[ "$ROW" == *syn_* ]] && DEPS+=$EX_SYN
        [[ "$ROW" == *ecoli_* || "$ROW" == *hp_w* ]] && DEPS+=$EX_BACT
        [[ "$ROW" == *anabaena* ]] && DEPS+=$EX_NOVEL
        J=$(sub ${DEPS:+--dependency=afterany$DEPS} --gres=gpu:1 --cpus-per-task=6 --mem=64G --time=03:00:00 --job-name=ty_$NAME --output=$RUN/logs/tr_${NAME}_%j.log --wrap="$PRE; $(train_cmd "$ROW" $RUN/runs/$NAME)")
        echo -e "train\t$NAME\t$J" >> $RUN/jobs.tsv
    done < <(experiments)
    echo "submitted $(grep -c . $RUN/jobs.tsv) jobs ($(awk -F'\t' '$3==""' $RUN/jobs.tsv | wc -l) failed to submit):"; column -t -s$'\t' $RUN/jobs.tsv
    echo; echo "=== scheduler start estimates (worst case; backfill usually beats them)"; squeue -u $USER --start -o "%.9i %.26j %.3t %.20S %R" | head -45 ;;

# =========================================================================== status
status)
    eval "$ENVACT"
    echo "=== queue"; squeue -u $USER -o "%.9i %.28j %.3t %.10M %R" | head -40
    echo; echo "=== data jobs"; cat $RUN/status/basecall_*.txt $RUN/status/extract_*.txt 2>/dev/null
    echo; echo "=== experiments"
    mkdir -p $REPO/benchmark_results/typing_overnight
    python - $RUN/runs $REPO/benchmark_results/typing_overnight <<'PY'
import glob, json, os, shutil, sys
runs, dest = sys.argv[1], sys.argv[2]
rows = []
for d in sorted(glob.glob(runs + "/E*")):
    nm = os.path.basename(d); f = d + "/metrics.json"
    if not os.path.exists(f): rows.append((nm, "not finished", "", "", "", "", "")); continue
    m = json.load(open(f)); shutil.copy(f, f"{dest}/{nm}.json")
    t = m["sets"].get("test", {})
    rows.append((nm, m["status"][:34], t.get("read_macro_f1", ""), t.get("read_auroc_mod_vs_none", ""), t.get("site_acc", ""),
                 (t.get("stoichiometry") or {}).get("mae", ""), json.dumps(t.get("read_recall_by_class", ""))[:70]))
hdr = ("experiment", "status", "read_macroF1", "read_modAUROC", "site_acc", "stoich_MAE", "per-read recall by class (test split)")
w = [max(len(str(r[i])) for r in rows + [hdr]) for i in range(7)]
for r in [hdr] + rows: print("  ".join(str(v).ljust(w[i]) for i, v in enumerate(r)))
open(dest + "/summary.tsv", "w").write("\n".join("\t".join(str(v) for v in r) for r in [hdr] + rows) + "\n")
print(f"\nmetrics copied to {dest} (git add benchmark_results/typing_overnight to send them back)")
PY
    ;;
print)   # every command that submit would run, for review and for the record
    while IFS= read -r ROW; do extract_cmd "$ROW" $N/${ROW%%|*}.npz; done < <(samples)
    while IFS= read -r ROW; do train_cmd "$ROW" $RUN/runs/${ROW%%|*}; done < <(experiments) ;;
*) echo "MODE must be disk | preflight | submit | status | print";;
esac
