#!/bin/bash
# Pull the first N pod5_pass files of the public GIAB HG002 R10.4.1 flowcell (ONT open data, giab_2023.05,
# 20230424_1302_3H_PAO89685_2264ba8c) over HTTPS, then hand them to basecall_moves.sh. Bhargav (Sep 27): 10x
# genome-wide is enough, data lives under storm/shared/data. Each file is ~0.03x (16 files = 0.46x) and
# the run has 1,912 files (~55x); 10x is ~345 files, so the default NFILES=400 (~11.6x, ~260 GB) is the target.
# Login node (needs internet). Resumable: wget -c skips finished files; run again after any interruption.
#   DEST=/fs/cbcb-lab/storm/shared/data/human_hg002_r10.4.1_giab_2023.05_PAO89685_10x bash download_hg002.sh
#   MODE=status DEST=... bash download_hg002.sh
set -uo pipefail
DEST=${DEST:?}; NFILES=${NFILES:-400}; MODE=${MODE:-run}; PAR=${PAR:-4}
PREFIX="giab_2023.05/flowcells/hg002/20230424_1302_3H_PAO89685_2264ba8c/pod5_pass/"
BUCKET="https://ont-open-data.s3.amazonaws.com"
REF=/fs/cbcb-lab/storm/shared/data/human_hg002_r10.4.1_ont_open_data/ref.fa       # GRCh38, the reference the ground truth is on
mkdir -p $DEST/pod5_files $DEST/logs
case $MODE in
run)
    if [ ! -s $DEST/lists/all_keys.txt ]; then
        mkdir -p $DEST/lists; : > $DEST/lists/all_keys.txt; TOKEN=""
        while :; do                                                    # S3 list-objects-v2 over plain HTTPS, 1000 keys a page
            URL="$BUCKET/?list-type=2&prefix=$PREFIX&max-keys=1000${TOKEN:+&continuation-token=$TOKEN}"
            curl -sf "$URL" > $DEST/lists/page.xml || { echo "listing failed: $URL"; exit 1; }
            grep -o '<Key>[^<]*\.pod5</Key>' $DEST/lists/page.xml | sed 's/<Key>//; s/<\/Key>//' >> $DEST/lists/all_keys.txt
            TOKEN=$(grep -o '<NextContinuationToken>[^<]*' $DEST/lists/page.xml | sed 's/<NextContinuationToken>//' | sed 's/+/%2B/g; s/=/%3D/g; s/\//%2F/g')
            [ -n "$TOKEN" ] || break
        done
        sort -o $DEST/lists/all_keys.txt $DEST/lists/all_keys.txt
    fi
    TOTAL=$(wc -l < $DEST/lists/all_keys.txt); head -n $NFILES $DEST/lists/all_keys.txt > $DEST/lists/wanted.txt
    echo "$TOTAL pod5 files in the run, taking the first $NFILES -> $DEST/pod5_files"
    cat $DEST/lists/wanted.txt | xargs -P $PAR -I{} sh -c 'wget -q -c -O "'$DEST'/pod5_files/$(basename {})" "'$BUCKET'/{}" || echo "FAILED {}"' 2>&1 | tee -a $DEST/logs/download.log
    n=$(ls $DEST/pod5_files/*.pod5 2>/dev/null | wc -l); echo "downloaded: $n of $NFILES files, $(du -sh $DEST/pod5_files | cut -f1)" | tee -a $DEST/logs/download.log
    [ $n -eq $NFILES ] && echo "next: NAME=hg002_10x POD5=$DEST/pod5_files REF=<indexed copy of $REF> OUT=$DEST/basecalled NJOBS=8 bash $(dirname $0)/basecall_moves.sh" ;;
status)
    echo "files: $(ls $DEST/pod5_files/*.pod5 2>/dev/null | wc -l) of $(wc -l < $DEST/lists/wanted.txt 2>/dev/null || echo ?), $(du -sh $DEST/pod5_files 2>/dev/null | cut -f1)"; tail -3 $DEST/logs/download.log 2>/dev/null ;;
*) echo "MODE=run|status"; exit 1 ;;
esac
