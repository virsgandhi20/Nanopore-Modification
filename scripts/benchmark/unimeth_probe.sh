#!/bin/bash
# Why does UniMeth's output not depend on its input? (login-node safe, ~2 min)
# All six diagnostic configs were blind, and switching normalization modes
# moved mean P(mod) only in the 4th decimal, so the model is not responding to
# the signal. Distinguish "checkpoint not loaded" from "signal not reaching the
# model", and look for UniMeth's own demo data to use as a positive control.
set -uo pipefail
D=/fs/nexus-scratch/vgandhi/unimeth_bench/diag
U=/fs/nexus-scratch/vgandhi/Unimeth
CK=/fs/nexus-scratch/vgandhi/unimeth_models/checkpoints
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
R=$REPO/benchmark_results/unimeth/probe.txt
source /nfshomes/vgandhi/miniconda3/etc/profile.d/conda.sh 2>/dev/null
conda activate /fs/nexus-scratch/vgandhi/envs/unimeth 2>/dev/null
q() { sort -g | awk '{a[NR]=$1} END{n=NR; m=split("0.01 0.5 0.9 0.99 0.999 0.9999",p," "); for(i=1;i<=m;i++){k=int(n*p[i]); if(k<1)k=1; printf "  q%s=%s", p[i], a[k]} printf "  max=%s  n=%d\n", a[n], n}'; }
{
echo "################ 1. runlog A, everything except progress bars"
grep -avE "record/s|s/ record|it/s\]" $D/A_d140_auto_6mA.runlog | cut -c1-220 | head -70
echo
echo "################ 2. load/key/warning lines across all runlogs"
grep -aihE "missing|unexpected|strict|warn|not found|mismatch|fallback|randomly|initializ|skip" $D/*.runlog | cut -c1-200 | sort | uniq -c | sort -rn | head -25
echo
echo "################ 3. is P(mod) a constant? quantiles"
echo "6mA (A):";  awk '$2>=0{print $9}' $D/A_d140_auto_6mA.txt | q
echo "5mC (E):";  awk '$2>=0{print $9}' $D/E_d140_auto_5mC.txt | q
echo
echo "################ 4. does the 5mC output depend on SEQUENCE context only?"
awk '$2>=0{s[$7]+=$9; n[$7]++} END{for(k in n) printf "  %-8s n=%-8d meanP=%.4f\n", k, n[k], s[k]/n[k]}' $D/E_d140_auto_5mC.txt
echo "   (different constants per context + no GT separation = model sees tokens, not signal)"
echo
echo "################ 5. same read, same position: A (auto norm) vs B (legacy norm)"
paste <(awk '$2>=0' $D/A_d140_auto_6mA.txt | head -200000 | cut -f1,2,5,9) <(awk '$2>=0' $D/B_d140_legacy_6mA.txt | head -200000 | cut -f9) \
  | awk '{d=$4-$5; if(d<0)d=-d; s+=d; if(d>m)m=d; n++} END{printf "  mean |P_auto - P_legacy| = %.6f   max = %.4f   over %d calls\n", s/n, m, n}'
echo "   (~0 means the signal is NOT influencing the output)"
echo
echo "################ 6. do BAM read ids exist in the pod5?"
SAM=/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools
B=/fs/nexus-scratch/vgandhi/unimeth_bench/bam/Ecoli_WT_5kHz.moves.bam
P=/fs/cbcb-lab/storm/bds062/data/benchmark/bacteria/Ecoli_WT_5kHz/pod5/Ecoli_WT_5kHz.pod5
$SAM view $B | head -300 | cut -f1 | sort -u > /tmp/um_ids_$$.txt
python - <<PY
import pod5
ids = [l.strip() for l in open("/tmp/um_ids_$$.txt")]
with pod5.Reader("$P") as r:
    have = {str(x) for x in r.read_ids}
hit = sum(1 for i in ids if i in have)
print(f"  {hit} of {len(ids)} BAM read ids found in the pod5 ({len(have):,} reads in pod5)")
PY
rm -f /tmp/um_ids_$$.txt
echo "  tags on first BAM record:"; $SAM view $B | head -1 | cut -f12- | tr '\t' '\n' | cut -c1-40 | grep -E "^(mv|ts|ns|sm|sd|pi|sp|st|MN):" | sed 's/^/    /'
echo
echo "################ 7. checkpoint structure"
python - <<PY
import torch
ck = torch.load("$CK/unimeth_r10.4.1_5kHz_6mA.pt", map_location="cpu", weights_only=False)
print("  type:", type(ck).__name__)
if isinstance(ck, dict):
    ks = list(ck.keys()); print("  top-level keys:", ks[:12], "... total", len(ks))
    sd = ck.get("state_dict") or ck.get("model") or ck.get("model_state_dict") or ck
    if isinstance(sd, dict):
        names = [k for k in sd.keys() if hasattr(sd[k], "shape")]
        print("  tensors:", len(names)); [print("   ", n, tuple(sd[n].shape)) for n in names[:6]]; print("    ..."); [print("   ", n, tuple(sd[n].shape)) for n in names[-4:]]
    for k in ("config", "args", "model_config", "hparams", "version"):
        if k in ck: print(f"  {k}:", str(ck[k])[:400])
PY
echo
echo "################ 8. demo / test data shipped with UniMeth (positive control?)"
find $U -maxdepth 4 \( -iname "*demo*" -o -name "*.pod5" -o -name "*.bam" -o -iname "*example*" -o -iname "test*" \) -not -path "*/.git/*" | head -20
grep -n -iE "demo|subset_18|example data|test data|zenodo|figshare|drive.google" $U/README.md | cut -c1-200 | head -15
echo
echo "################ 9. flags we have not used"
unimeth-infer --help 2>&1 | grep -A3 -E -- "--signal_index|--bam_mode|--model_type|--keep_mv" | cut -c1-160
} > $R 2>&1
echo "wrote $R ($(wc -l < $R) lines)"
