#!/bin/bash
# Follow-up to the overnight per-read typing run (2026-09-21). Reuses the
# extracted windows of that run; nothing is basecalled or extracted again.
#
#   seeds     E01 / E04 / E07 / E10 with two more seeds (is a difference real?)
#   N-series  oligo <-> genome transfer with each window re-normalized by its own
#             median / MAD, and with a narrower window (is the transfer failure a
#             read-level normalization artefact? does less context help?)
#   M-series  E. coli with ONLY motif sites, so 'none' comes from the matched
#             dam-/dcm- controls alone (removes the class-prior bias that made
#             unmodified reads at motif sites look modified)
#   P         the same typing question asked of RawMod's frozen encoder
#             (rawmod_embedding_probe.py on Bhargav's strand-resolved features)
#
# MODE=preflight (login node, nothing queued) | submit | status
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
ME=/fs/nexus-scratch/vgandhi
SRC=${SRC:-$ME/typing_overnight_scav}          # npz from the overnight run
FU=${FU:-$ME/typing_followup}
N=$SRC/npz
FEAT=/fs/cbcb-lab/storm/bds062/rawmod_strand_resolved/features
ENVACT="source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh; conda activate $ME/envs/unimeth"
SBASE=${SBASE:-"--account=scavenger --partition=scavenger --qos=scavenger --requeue"}
GRES=${GRES:-gpu:rtxa5000:1}
MODE=${MODE:-status}
mkdir -p $FU/{runs,logs,code,status} 2>/dev/null

ECO="$N/ecoli_wt.npz,$N/ecoli_dm.npz,$N/ecoli_msssi.npz"; HP="$N/hp_wt.npz,$N/hp_wga.npz"
SYN=$(for R in rep1 rep2; do for C in control 5mC 5hmC 6mA; do printf "$N/syn_${C}_$R.npz,"; done; done | sed 's/,$//')
NOVEL="anabaena=$N/anabaena.npz,tdent=$N/tdent.npz,j99=$N/j99.npz"
# name | train | classes | inputs | extra | note
experiments() {
for S in 1 2; do cat <<EOF
E01_syn_type4_all_s$S|$SYN|none,5mC,5hmC,6mA|sig,dwell,seq|--seed $S|seed replicate of E01
E04_syn_5mC_vs_5hmC_s$S|$SYN|none,5mC,5hmC|sig,dwell,seq|--groups synC --seed $S|seed replicate of E04
E07_ecoli_type3_all_s$S|$ECO|none,6mA,5mC|sig,dwell,seq|--seed $S|seed replicate of E07
E10_bact_type4_all_s$S|$ECO,$HP|none,6mA,5mC,4mC|sig,dwell,seq|--eval $NOVEL --seed $S|seed replicate of E10 (incl. unseen organisms)
EOF
done
cat <<EOF
N01_syn_type4_winnorm|$SYN|none,5mC,5hmC,6mA|sig,dwell,seq|--renorm window|does per-window normalization cost anything in-domain
N01_syn_type4_crop5|$SYN|none,5mC,5hmC,6mA|sig,dwell,seq|--crop 5|11-base window: less context to memorise, better on unseen 5-mers?
N13_syn_to_bact_winnorm|$SYN|none,5mC,5hmC,6mA|sig,dwell,seq|--renorm window --eval ecoli=${ECO//,/+},hp=${HP//,/+},$NOVEL|oligo -> genome with window normalization
N13_syn_to_bact_winnorm_crop5|$SYN|none,5mC,5hmC,6mA|sig,dwell|--renorm window --crop 5 --eval ecoli=${ECO//,/+},hp=${HP//,/+},$NOVEL|same, 11-base window, signal only
N14_bact_to_syn_winnorm|$ECO,$HP|none,6mA,5mC,4mC|sig,dwell,seq|--renorm window --eval syn=${SYN//,/+}|genome -> oligo with window normalization
N14_bact_to_syn_winnorm_crop5|$ECO,$HP|none,6mA,5mC,4mC|sig,dwell|--renorm window --crop 5 --eval syn=${SYN//,/+}|same, 11-base window, signal only
M07_ecoli_motifonly_all|$ECO|none,6mA,5mC|sig,dwell,seq|--groups dam,dcm,cpg|'none' only from matched controls at the same motif sites
M08_ecoli_motifonly_signal|$ECO|none,6mA,5mC|sig,dwell|--groups dam,dcm,cpg|same without sequence input
M10_bact_motifonly_signal|$ECO,$HP|none,6mA,5mC,4mC|sig,dwell|--groups dam,dcm,cpg,hp6mA,hp5mC,hp4mC,d6mA,dneg --eval $NOVEL|all bacteria, matched controls only, signal only
EOF
}
train_cmd() { IFS='|' read -r NAME TRAIN CLASSES INPUTS EXTRA NOTE <<< "$1"
    echo "python $FU/code/scripts/typing/train_read_typing.py --name $NAME --out $2 --train $TRAIN --classes $CLASSES --inputs $INPUTS $EXTRA --note \"$NOTE\" ${3:-}"; }
CK=$FU/code/checkpoints/results20_sad_dim16
probe_cmd() { echo "\$PYP $FU/code/scripts/typing/rawmod_embedding_probe.py --repo $FU/code --features $FEAT --checkpoints mixed=$CK/mixed.pt,loco_5hmC=$CK/loco_5hmC.pt,loco_5mC=$CK/loco_5mC.pt,loco_6mA=$CK/loco_6mA.pt --out $1 ${2:-}"; }
# the probe imports RawMod's training code (needs matplotlib); fall back to Bhargav's env if mine lacks it
PICKPY="if python -c 'import matplotlib, h5py, sklearn, torch' 2>/dev/null; then PYP=python; else PYP=/fs/nexus-scratch/bds062/envs/mod/bin/python; fi"
snapshot() { rm -rf $FU/code && mkdir -p $FU/code && git -C $REPO archive HEAD scripts rawmod checkpoints | tar -x -C $FU/code; }

case $MODE in
preflight)
    FAILS=0; ok() { echo "  PASS  $1"; }; bad() { echo "  FAIL  $1"; FAILS=$((FAILS+1)); }
    snapshot && ok "code + checkpoints snapshot ($(ls $CK | wc -l) checkpoints)" || bad "snapshot"
    FREE=$(df --output=avail -BG $ME | tail -1 | tr -dc 0-9); [ "${FREE:-0}" -ge 5 ] && ok "free space ${FREE}G" || bad "only ${FREE}G free"
    for f in $(echo $SYN,$ECO,$HP | tr ',' ' ') $N/anabaena.npz $N/tdent.npz $N/j99.npz; do [ -s $f ] || bad "missing $f"; done; ok "checked 16 input npz"
    eval "$ENVACT"; eval "$PICKPY"; echo "  info  probe python: $PYP"
    mkdir -p $FU/pf
    echo "--- trainer with the new options, one epoch on a slice of real data (login node CPU)"
    python - <<PY
import numpy as np
for s in ("syn_control_rep1", "syn_5mC_rep1", "ecoli_wt", "ecoli_dm"):
    z = np.load("$N/%s.npz" % s); n = len(z["label"]); keep = np.sort(np.random.default_rng(0).choice(n, min(n, 6000), replace=False))
    out = {k: (z[k][keep] if k in ("sig", "dwell", "base", "label", "group", "site", "read") else z[k]) for k in z.files}
    np.savez_compressed("$FU/pf/%s.npz" % s, **out)
PY
    P=$FU/pf
    pf() { local nm=$1; shift; eval "python $FU/code/scripts/typing/train_read_typing.py --epochs 1 --batch 256 --name pf_$nm --out $P/runs/$nm $*" > $P/$nm.log 2>&1
           local st=$(cut -f2 $P/runs/$nm/status.txt 2>/dev/null); case "$st" in OK*) ok "train $nm: $st $(grep -m1 '\[test\]' $P/$nm.log | cut -c1-100)";; *) bad "train $nm: ${st:-no status}"; tail -5 $P/$nm.log;; esac; }
    pf winnorm_crop --train $P/syn_control_rep1.npz,$P/syn_5mC_rep1.npz --classes none,5mC --inputs sig,dwell --renorm window --crop 5 --seed 1 --eval eco=$P/ecoli_wt.npz
    pf motifonly    --train $P/ecoli_wt.npz,$P/ecoli_dm.npz --classes none,6mA,5mC --groups dam,dcm,cpg
    echo "--- RawMod embedding probe on the real strand-resolved features (40 images per class, CPU)"
    for c in 5mC 5hmC 6mA control; do [ -r $FEAT/ONT_${c}_plus.h5 ] && ok "readable ONT_${c}_plus.h5" || bad "cannot read $FEAT/ONT_${c}_plus.h5"; done
    eval "$(probe_cmd $P/probe "--max-per-class 40")" > $P/probe.log 2>&1
    st=$(cut -f2 $P/probe/status.txt 2>/dev/null); case "$st" in OK*) ok "probe: $st"; grep -E "images at|acc=" $P/probe.log | head -12 | sed 's/^/        /';; *) bad "probe: ${st:-no status}"; tail -8 $P/probe.log;; esac
    sbatch --test-only $SBASE --gres=$GRES --cpus-per-task=6 --mem=64G --time=03:00:00 --wrap=true > $P/sb.txt 2>&1 && ok "job template: $(tr '\n' ' ' < $P/sb.txt | cut -c1-80)" || bad "job template: $(cat $P/sb.txt)"
    echo; [ $FAILS -eq 0 ] && echo "PREFLIGHT: ALL PASS - safe to run MODE=submit" || echo "PREFLIGHT: $FAILS FAILURE(S) - paste this output"
    exit $FAILS ;;
submit)
    snapshot; : > $FU/jobs.tsv
    PRE="$ENVACT; set -uo pipefail"
    while IFS= read -r ROW; do NAME=${ROW%%|*}
        J=$(sbatch --parsable $SBASE --gres=$GRES --cpus-per-task=6 --mem=64G --time=03:00:00 --job-name=fu_$NAME --output=$FU/logs/${NAME}_%j.log --wrap="$PRE; $(train_cmd "$ROW" $FU/runs/$NAME)")
        echo -e "train\t$NAME\t${J%%;*}" >> $FU/jobs.tsv
    done < <(experiments)
    J=$(sbatch --parsable $SBASE --gres=$GRES --cpus-per-task=6 --mem=64G --time=03:00:00 --job-name=fu_P01_rawmod_probe --output=$FU/logs/P01_rawmod_probe_%j.log --wrap="$PRE; $PICKPY; $(probe_cmd $FU/runs/P01_rawmod_probe)")
    echo -e "probe\tP01_rawmod_probe\t${J%%;*}" >> $FU/jobs.tsv
    echo "submitted $(grep -c . $FU/jobs.tsv) jobs:"; column -t -s$'\t' $FU/jobs.tsv; squeue -u $USER -o "%.9i %.34j %.3t %.8M %R" | head -30 ;;
status)
    eval "$ENVACT"; squeue -u $USER -o "%.9i %.34j %.3t %.8M %R" | head -30
    mkdir -p $REPO/benchmark_results/typing_followup
    python - $FU/runs $REPO/benchmark_results/typing_followup <<'PY'
import glob, json, os, shutil, sys
runs, dest = sys.argv[1:3]
print(f"\n{'experiment':36s} {'status':18s} {'macroF1':>8s} {'modAUROC':>9s} {'site_acc':>9s}  held-out sets (macroF1 / modAUROC)")
for d in sorted(glob.glob(runs + "/*")):
    nm = os.path.basename(d); f = d + "/metrics.json"
    if not os.path.exists(f): print(f"{nm:36s} not finished"); continue
    m = json.load(open(f)); shutil.copy(f, f"{dest}/{nm}.json")
    if nm.startswith("P01"):
        print(f"{nm:36s} {m['status'][:18]:18s}")
        for ck, r in m.get("checkpoints", {}).items():
            print(f"    {ck:12s} site_rep F1={r['site_rep_96d']['macro_f1']}  per-read post-TF F1={r['read_post_transformer']['macro_f1']} (pooled acc {r['read_post_transformer'].get('pooled_by_image_acc')})"
                  f"  pre-TF F1={r['read_pre_transformer']['macro_f1']}  raw baseline F1={r['read_raw_signal_baseline']['macro_f1']}"
                  + (f"  | 5mC-vs-5hmC per-read acc={r['pair_5mC_vs_5hmC']['read_post_transformer']['acc']} (raw {r['pair_5mC_vs_5hmC']['read_raw_signal_baseline']['acc']})" if "pair_5mC_vs_5hmC" in r else ""))
        continue
    t = m["sets"].get("test", {})
    ev = "  ".join(f"{k}={v.get('read_macro_f1')}/{v.get('read_auroc_mod_vs_none')}" for k, v in m["sets"].items() if not k.startswith("test") and isinstance(v, dict) and "read_macro_f1" in v)
    print(f"{nm:36s} {m['status'][:18]:18s} {str(t.get('read_macro_f1')):>8s} {str(t.get('read_auroc_mod_vs_none')):>9s} {str(t.get('site_acc')):>9s}  {ev}")
print(f"\nmetrics copied to {dest}")
PY
    ;;
print) while IFS= read -r ROW; do train_cmd "$ROW" $FU/runs/${ROW%%|*}; done < <(experiments); probe_cmd $FU/runs/P01_rawmod_probe ;;
*) echo "MODE must be preflight | submit | status | print";;
esac
