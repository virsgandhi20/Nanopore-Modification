#!/bin/bash
# Pull the first N pod5_pass files of the public GIAB HG002 R10.4.1 flowcell (ONT open data, giab_2023.05,
# 20230424_1302_3H_PAO89685_2264ba8c) over HTTPS, then hand them to basecall_moves.sh. Bhargav (Sep 27): 10x
# genome-wide is enough, data lives under storm/shared/data. Each file is ~0.02-0.05x (24 GB of pod5 ~ 1x); the run has
# 1,912 files (~55x). NFILES=500 (~240 GB) lands around 10x. Login node (needs internet).
# Safe by construction: one instance at a time (lock), every file is fetched to a .part name and renamed only when
# its size equals the size in the bucket listing, oversized or corrupt files are deleted and refetched, rerunning
# resumes. Status/verify count only files whose size is right.
#   DEST=/fs/cbcb-lab/storm/shared/data/human_hg002_r10.4.1_giab_2023.05_PAO89685_10x [NFILES=500] bash download_hg002.sh
#   MODE=status|verify DEST=... bash download_hg002.sh
set -uo pipefail
DEST=${DEST:?}; NFILES=${NFILES:-500}; MODE=${MODE:-run}; PAR=${PAR:-4}
PREFIX="giab_2023.05/flowcells/hg002/20230424_1302_3H_PAO89685_2264ba8c/pod5_pass/"
BUCKET="https://ont-open-data.s3.amazonaws.com"
REF=/fs/cbcb-lab/storm/shared/data/human_hg002_r10.4.1_ont_open_data/ref.fa       # GRCh38, the reference the ground truth is on
mkdir -p $DEST/pod5_files $DEST/logs $DEST/lists
list() {   # key <TAB> size for every pod5 in the run, from the S3 list-objects-v2 XML (1000 keys a page)
    [ -s $DEST/lists/keys_sizes.tsv ] && return
    : > $DEST/lists/keys_sizes.tsv; TOKEN=""
    while :; do
        URL="$BUCKET/?list-type=2&prefix=$PREFIX&max-keys=1000${TOKEN:+&continuation-token=$TOKEN}"
        curl -sf "$URL" > $DEST/lists/page.xml || { echo "listing failed: $URL"; exit 1; }
        grep -o '<Contents>.*</Contents>' $DEST/lists/page.xml | sed 's/<\/Contents>/\n/g' | grep -o '<Key>[^<]*\.pod5</Key><LastModified>[^<]*</LastModified><ETag>[^<]*</ETag><Size>[0-9]*' \
            | sed 's/<Key>//; s/<\/Key>.*<Size>/\t/' >> $DEST/lists/keys_sizes.tsv
        TOKEN=$(grep -o '<NextContinuationToken>[^<]*' $DEST/lists/page.xml | sed 's/<NextContinuationToken>//' | sed 's/+/%2B/g; s/=/%3D/g; s/\//%2F/g')
        [ -n "$TOKEN" ] || break
    done
    sort -o $DEST/lists/keys_sizes.tsv $DEST/lists/keys_sizes.tsv
}
wanted() { list; head -n $NFILES $DEST/lists/keys_sizes.tsv > $DEST/lists/wanted.tsv; }
count_good() {   # files present with the right size; prints "good total bytes"
    local good=0 bytes=0 k s f
    while IFS=$'\t' read -r k s; do f=$DEST/pod5_files/$(basename $k); [ "$(stat -c %s $f 2>/dev/null)" = "$s" ] && { good=$((good + 1)); bytes=$((bytes + s)); }; done < $DEST/lists/wanted.tsv
    echo "$good $bytes"
}
fetch_one() {   # $1 key, $2 expected size
    local k=$1 s=$2 f=$DEST/pod5_files/$(basename $1)
    [ "$(stat -c %s $f 2>/dev/null)" = "$s" ] && return 0
    [ -e $f ] && rm -f $f                                                   # wrong size: corrupt or oversized, refetch
    [ -e $f.part ] && [ "$(stat -c %s $f.part)" -gt "$s" ] && rm -f $f.part  # a .part larger than the file cannot be resumed
    wget -q -c -O $f.part "$BUCKET/$k" && [ "$(stat -c %s $f.part 2>/dev/null)" = "$s" ] && mv $f.part $f && return 0
    echo "FAILED $(basename $k) (size $(stat -c %s $f.part 2>/dev/null || echo none) of $s)"; return 1
}
export -f fetch_one; export DEST BUCKET
case $MODE in
run)
    exec 9> $DEST/.download.lock; flock -n 9 || { echo "another download_hg002.sh is running on $DEST (lock held); not starting a second one"; exit 1; }
    wanted; read good bytes < <(count_good)
    echo "$(wc -l < $DEST/lists/keys_sizes.tsv) pod5 files in the run, taking the first $NFILES; $good already complete" | tee -a $DEST/logs/download.log
    cut -f1,2 $DEST/lists/wanted.tsv | xargs -P $PAR -n 2 bash -c 'fetch_one "$0" "$1"' 2>&1 | tee -a $DEST/logs/download.log
    read good bytes < <(count_good); echo "complete: $good of $NFILES files, $((bytes / 1000000000)) GB (~$((bytes / 24000000000))x)" | tee -a $DEST/logs/download.log
    [ $good -eq $NFILES ] && echo "next: NAME=hg002_10x POD5=$DEST/pod5_files REF=/fs/cbcb-lab/storm/vgandhi/matrix/hg002/ref.fa OUT=$DEST/basecalled NJOBS=32 PART=cbcb bash $(dirname $0)/basecall_moves.sh" ;;
status|verify)
    wanted; read good bytes < <(count_good)
    echo "complete: $good of $NFILES files, $((bytes / 1000000000)) GB (~$((bytes / 24000000000))x); other files present: $(ls $DEST/pod5_files 2>/dev/null | wc -l) total, $(ls $DEST/pod5_files/*.part 2>/dev/null | wc -l) partial"
    pgrep -fa "download_hg002|wget -q -c" | grep -v pgrep | head -3; tail -2 $DEST/logs/download.log 2>/dev/null ;;
*) echo "MODE=run|status|verify"; exit 1 ;;
esac
