#!/bin/bash
# Step 0 for the overnight typing runs: what am I allowed to run, and which
# inputs exist? Read-only, login-node safe, under a minute.
set -uo pipefail
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
BENCH=/fs/cbcb-lab/storm/bds062/data/benchmark
ONT=/fs/nexus-scratch/bds062/data/ont-os
ME=/fs/nexus-scratch/vgandhi
has_mv() { $SAM view "$1" 2>/dev/null | head -200 | grep -c "mv:B" ; }

echo "################ 1. SLURM limits"
echo "-- my associations (account|partition|qos|maxjobs|maxsubmit|grptres)"
sacctmgr -nP show assoc user=$USER format=account,partition,qos,maxjobs,maxsubmit,grptres 2>&1 | head -12
echo "-- qos limits (name|MaxTRESPerUser|MaxJobsPerUser|MaxSubmitPerUser|MaxWall|MaxTRESPerJob|GrpTRES)"
sacctmgr -nP show qos format=name,maxtrespu,maxjobspu,maxsubmitpu,maxwall,maxtres,grptres 2>&1 | grep -E "^(high|medium|default|huge-long|highmem|scavenger|cpu|gpu)" | head -12
echo "-- cbcb partition nodes (state gres cpus mem)"
sinfo -p cbcb -N -o "%N %T %G %c %m" 2>&1 | sort -u | head -30
echo "-- my queue now"; squeue -u $USER 2>&1 | head
echo "-- pending GPU jobs in cbcb: $(squeue -p cbcb -t PD -h 2>/dev/null | wc -l)   running: $(squeue -p cbcb -t R -h 2>/dev/null | wc -l)"

echo; echo "################ 2. disk"
df -h $ME 2>&1 | tail -1; du -sh $ME/unimeth_bench $ME/spo1_typing $ME/orca_feat 2>/dev/null

echo; echo "################ 3. benchmark samples (pod5 size, my move-table BAM)"
for d in $BENCH/bacteria/*/; do S=$(basename $d)
    P=$(ls $d/pod5/*.pod5 2>/dev/null | head -1); B=$ME/unimeth_bench/bam/$S.moves.bam
    printf "  %-24s pod5=%-7s movesBAM=%s\n" $S "$( [ -n "$P" ] && du -h $P | cut -f1 || echo none)" "$( [ -s $B ] && echo "yes($(du -h $B | cut -f1))" || echo no)"
done
ls $BENCH/ 2>/dev/null | tr '\n' ' '; echo

echo; echo "################ 4. ONT synthetic (pure-chemistry samples = exact per-read labels)"
ls -la $ONT/ 2>&1 | head -12
echo "-- subset pod5:";  ls -la $ONT/subset/ 2>&1 | head -14
echo "-- basecalls:";    ls -la $ONT/basecalls/ 2>&1 | head -14
echo "-- references:";   ls -la $ONT/references/ 2>&1 | head -14
for b in $ONT/basecalls/*.bam; do [ -r "$b" ] && echo "  $(basename $b): mv tags in first 200 records = $(has_mv $b), header: $($SAM view -H $b 2>/dev/null | grep -m1 '@PG' | cut -c1-160)"; done 2>/dev/null | head -10
for f in $ONT/references/all_5mers_*_sites.bed; do [ -r "$f" ] && echo "  $(basename $f): $(wc -l < $f) sites, first: $(head -1 $f | tr '\t' ' ')"; done
grep -c ">" $ONT/references/all_5mers.fa 2>/dev/null | sed 's/^/  all_5mers.fa contigs: /'

echo; echo "################ 5. SPO1 re-basecalled BAMs (per-read 5mC/5hmC/6mA calls)"
for b in $ME/spo1_typing/*.mod.sorted.bam; do [ -r "$b" ] && echo "  $(basename $b): $(du -h $b | cut -f1), mv=$(has_mv $b), MM=$($SAM view $b | head -200 | grep -c 'MM:Z')"; done

echo; echo "################ 6. H. pylori differential GT + existing site-level features"
for f in $ME/hp_labels/gt_*.bed; do echo "  $(basename $f): $(wc -l < $f) rows, first: $(head -1 $f | tr '\t' ' ')"; done
ls -d $ME/orca_feat/*/ $ME/orca_feat_spo1/run1_jan31/single_end/*/ $ME/orca_feat_bench/*/ 2>/dev/null | head -20

echo; echo "################ 7. python environment"
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null
python - <<'PY'
import importlib
for m in ("torch", "numpy", "sklearn", "pysam", "pod5", "h5py", "pandas", "scipy"):
    try:
        x = importlib.import_module(m); print(f"  {m:8s} {getattr(x, '__version__', '?')}")
    except Exception as e: print(f"  {m:8s} MISSING ({e.__class__.__name__})")
PY
/fs/nexus-scratch/bds062/envs/mod/bin/python -c "import torch, h5py; print('  bhargav mod env: torch', torch.__version__, 'h5py', h5py.__version__)" 2>&1 | tail -1
