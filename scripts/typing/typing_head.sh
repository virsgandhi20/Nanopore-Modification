#!/bin/bash
# Per-read typing on RawMod's encoder (rawmod_typing_head.py): the experiment grid, as SLURM jobs on the cbcb partition.
#
#   bash typing_head.sh preflight   # self-test on the login node with the test environment (CPU, ~1 min)
#   bash typing_head.sh data        # which feature files exist (Bhargav's strand-resolved set, the benchmark pipeline's)
#   bash typing_head.sh submit [T01 T02 ...]   # sbatch the grid (default: every experiment whose inputs exist)
#                                              # T09 to T15 (Oct 31): the trained 'other' class and the Mahalanobis baseline
#   bash typing_head.sh status      # one line per experiment; copies metrics.json into the fork's benchmark_results/typing_head/
#   bash typing_head.sh print       # the commands, nothing submitted
#
# Reuses the Oct 28 repository-test setup: the released RawMod clone, its conda env (torch 2.14.1+cu126, which runs on
# both GPU generations of the partition) and the oligo features that the benchmark pipeline wrote.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
PY_SCRIPT=$HERE/rawmod_typing_head.py
ROOT=/fs/cbcb-lab/storm/vgandhi/rawmod-test
REPO=$ROOT/RawMod
ENV=$ROOT/env
CONDA_SH=$HOME/miniconda3/etc/profile.d/conda.sh
RUN=${RUN:-/fs/cbcb-lab/storm/vgandhi/typing_head}       # outputs (lab share: scratch is full)
FORK=${FORK:-$HOME/Nanopore-Modification}                 # cluster checkout of the fork, for the results copy
SRC_A=/fs/cbcb-lab/storm/bds062/rawmod_strand_resolved/features   # Bhargav's strand-resolved features (Sep 21 probe)
SRC_B=$ROOT/work                                                   # the benchmark pipeline's own features (Oct 28)
TAG=own1_k10_seed0_q10_r15
SLURM="--account=cbcb --partition=cbcb --qos=high --gres=gpu:1 --cpus-per-task=8 --mem=32G --time=03:00:00 --exclude=legacy43"
unset TMPDIR

fa() { local s=$1; [ "$s" = none ] && s=control; echo "$SRC_A/ONT_${s}_plus.h5"; }   # chem -> file, source A (control = none)
fb() { local s=$1; [ "$s" = none ] && s=control; echo "$SRC_B/syn_$s/rawmod_features/$TAG/features_plus.h5"; }   # source B
fb2() { local s=$1; [ "$s" = none ] && s=control; echo "$SRC_B/syn_${s}_rep2/rawmod_features/$TAG/features_plus.h5"; }
files() { local src=$1; shift; local out=""; for c in "$@"; do out+="${out:+,}$c=$($src "$c")"; done; echo "$out"; }
CK=$REPO/checkpoints

# name | checkpoint | mode | classes | files | eval-files | extra
declare -A EXP NOTE
EXP[T01]="$CK/region_holdout/mixed.pt|frozen|none,5mC,5hmC,6mA|$(files fa none 5mC 5hmC 6mA)||"
EXP[T02]="$CK/region_holdout/mixed.pt|finetune|none,5mC,5hmC,6mA|$(files fa none 5mC 5hmC 6mA)||"
EXP[T03]="$CK/region_holdout/mixed.pt|finetune|none,5mC,5hmC,6mA|$(files fb none 5mC 5hmC 6mA)|$(files fb2 none 5mC 5hmC 6mA)|"
EXP[T04]="$CK/leaveout/ont5hmC.pt|finetune|none,5mC,6mA|$(files fa none 5mC 5hmC 6mA)||"
EXP[T05]="$CK/leaveout/ont6mA.pt|finetune|none,5mC,5hmC|$(files fa none 5mC 5hmC 6mA)||"
EXP[T06]="$CK/region_holdout/loco_5hmC.pt|finetune|none,5mC,6mA|$(files fa none 5mC 5hmC 6mA)||"
EXP[T07]="$CK/rawmod.pt|frozen|none,5mC,5hmC,6mA|$(files fa none 5mC 5hmC 6mA)||"
EXP[T08]="$CK/region_holdout/mixed.pt|finetune|none,5mC,5hmC,6mA|$(files fa none 5mC 5hmC 6mA)||--detect-weight 0"
# Oct 31: a trained 'other' class (--other) after the confidence-based open set of T04 to T06 failed (AUROC 0.30 to 0.39).
# The other class is fed by one foreign chemistry and the unknown is a different one; the encoder never saw the unknown.
EXP[T09]="$CK/leaveout/ont5hmC.pt|finetune|none,5mC|$(files fa none 5mC 5hmC 6mA)||--other 6mA"
EXP[T10]="$CK/leaveout/ont6mA.pt|finetune|none,5mC|$(files fa none 5mC 5hmC 6mA)||--other 5hmC"
EXP[T11]="$CK/region_holdout/loco_5hmC.pt|finetune|none,6mA|$(files fa none 5mC 5hmC 6mA)||--other 5mC"
EXP[T12]="$CK/leaveout/ont5hmC.pt|frozen|none,5mC|$(files fa none 5mC 5hmC 6mA)||--other 6mA"
EXP[T13]="$CK/leaveout/ont5hmC.pt|finetune|none,5mC,6mA|$(files fa none 5mC 5hmC 6mA)||"
EXP[T14]="$CK/leaveout/ont6mA.pt|finetune|none,5mC,5hmC|$(files fa none 5mC 5hmC 6mA)||"
EXP[T15]="$CK/region_holdout/loco_5hmC.pt|finetune|none,6mA|$(files fa none 5mC 5hmC 6mA)||"
NOTE[T01]="frozen encoder, mixed.pt: the non-linear probe (compare with the Sep 21 linear probe)"
NOTE[T02]="fine-tuned encoder + head, mixed.pt, detection kept by BCE"
NOTE[T03]="as T02 on the benchmark pipeline's own features; replicate-2 libraries scored unseen"
NOTE[T04]="open set: 5hmC unknown, encoder never saw the 5hmC oligo library (leaveout/ont5hmC.pt)"
NOTE[T05]="open set: 6mA unknown, encoder never saw the 6mA oligo library (leaveout/ont6mA.pt)"
NOTE[T06]="open set: 5hmC unknown, encoder never saw 5hmC at all (region_holdout/loco_5hmC.pt)"
NOTE[T07]="frozen encoder, rawmod.pt (all data)"
NOTE[T08]="as T02 without the detection term: does detection degrade?"
NOTE[T09]="trained other: other=6mA, unknown=5hmC, encoder never saw the 5hmC library (ont5hmC.pt)"
NOTE[T10]="trained other: other=5hmC, unknown=6mA, encoder never saw the 6mA library (ont6mA.pt)"
NOTE[T11]="trained other: other=5mC, unknown=5hmC (the hard pair as other), encoder never saw 5hmC at all (loco_5hmC.pt)"
NOTE[T12]="as T09 with the encoder frozen: does the other class need fine-tuning?"
NOTE[T13]="T04 rerun for the Mahalanobis score: 5hmC unknown, no other class (ont5hmC.pt)"
NOTE[T14]="T05 rerun for the Mahalanobis score: 6mA unknown, no other class (ont6mA.pt)"
NOTE[T15]="T06 rerun for the Mahalanobis score, two known classes (none, 6mA) so it pairs with T11"
ORDER="T01 T02 T03 T04 T05 T06 T07 T08 T09 T10 T11 T12 T13 T14 T15"

inputs_exist() {   # every path of the files and eval-files fields
  IFS='|' read -r ck mode classes fl ev extra <<<"${EXP[$1]}"
  [ -f "$ck" ] || { echo "  $1: missing checkpoint $ck"; return 1; }
  for spec in ${fl//,/ } ${ev//,/ }; do [ -f "${spec#*=}" ] || { echo "  $1: missing ${spec#*=}"; return 1; }; done
}
cmd() {
  IFS='|' read -r ck mode classes fl ev extra <<<"${EXP[$1]}"
  echo "$ENV/bin/python $PY_SCRIPT --repo $REPO --checkpoint $ck --mode $mode --classes $classes --files $fl ${ev:+--eval-files $ev} $extra --name $1 --out $RUN/$1 --note '${NOTE[$1]}'"
}

case "${1:-}" in
  preflight)
    source "$CONDA_SH"; conda activate "$ENV"; mkdir -p "$RUN"
    python "$PY_SCRIPT" --selftest --repo "$REPO" --checkpoint "$CK/region_holdout/mixed.pt" --mode finetune --classes none,5mC,6mA --name preflight --out "$RUN/preflight" 2>&1 | grep -v "UserWarning\|self.encoder = " | tail -12
    echo "STAGE preflight OK";;
  data)
    echo "--- source A (Bhargav's strand-resolved features)"; ls -la "$SRC_A"/ONT_*_plus.h5 2>&1 | head
    echo "--- source B (benchmark pipeline features)"; ls -la "$SRC_B"/syn_*/rawmod_features/$TAG/features_plus.h5 2>&1 | head -12
    echo "--- experiments whose inputs exist"; for e in $ORDER; do inputs_exist "$e" && echo "  $e ok: ${NOTE[$e]}"; done; echo "STAGE data OK";;
  print)
    for e in ${2:-$ORDER}; do echo "== $e"; cmd "$e"; done;;
  submit)
    shift; mkdir -p "$RUN/logs"; sel=${*:-$ORDER}
    for e in $sel; do
      inputs_exist "$e" || continue
      jid=$(sbatch --parsable --export=ALL --job-name="ty_$e" --output="$RUN/logs/${e}_%j.out" $SLURM --wrap "$(cmd "$e")")
      echo "$jid  $e  ${NOTE[$e]}"
    done
    squeue -u "$USER" -o "%.10i %.12j %.9T %.8M %R" | grep -E 'JOBID|ty_' || true
    echo "STAGE submit OK";;
  status)
    printf '%-5s %-10s %-7s %-7s %-7s %-7s %-7s %-7s %-7s  %s\n' exp status readF1 imgAcc detBef detAft softAUC mahAUC pOthAUC note
    for e in $ORDER; do
      m="$RUN/$e/metrics.json"
      if [ -f "$m" ]; then
        python3 - "$m" "$e" "${NOTE[$e]}" <<'PY'
import json, sys
m = json.load(open(sys.argv[1])); t = m.get('test', {}); o = m.get('open_set', {})
mah = (o.get('mahalanobis') or {}).get('auroc_known_vs_unknown_reads', '-'); po = (o.get('p_other') or {}).get('auroc_known_vs_unknown_reads', '-')
print(f"{sys.argv[2]:<5} {m.get('status','?')[:10]:<10} {t.get('read_macro_f1','-')!s:<7} {t.get('image_acc','-')!s:<7} "
      f"{m.get('detection_auroc_test_before','-')!s:<7} {m.get('detection_auroc_test_after','-')!s:<7} {o.get('auroc_known_vs_unknown_reads','-')!s:<7} {mah!s:<7} {po!s:<7}  {sys.argv[3]}")
PY
      else
        st=$(squeue -h -u "$USER" -n "ty_$e" -o %T | head -1); printf '%-5s %-10s %s\n' "$e" "${st:-not run}" "${NOTE[$e]}"
      fi
    done
    if [ -d "$FORK/benchmark_results" ]; then
      mkdir -p "$FORK/benchmark_results/typing_head"
      for e in $ORDER; do [ -f "$RUN/$e/metrics.json" ] && cp "$RUN/$e/metrics.json" "$FORK/benchmark_results/typing_head/$e.json"; done
      echo "(metrics copied to $FORK/benchmark_results/typing_head/)"
    fi
    echo "--- recent failures"; grep -l "FAILED\|Traceback" "$RUN"/logs/*.out 2>/dev/null | tail -3 || true
    echo "STAGE status OK";;
  *) sed -n '2,12p' "$0"; exit 1;;
esac
