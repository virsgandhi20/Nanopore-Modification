#!/bin/bash
# Differential ground-truth labels for H. pylori 26695 from the benchmark
# modbams: a site is a positive for modification M when dorado calls it highly
# modified in the native (WT) sample AND essentially unmodified in the WGA
# sample. The WGA condition removes dorado's systematic false positives:
# WGA DNA carries no modifications, so anything called there is basecaller
# artifact, not biology.
#
# Output: $OUT/gt_{4mC,6mA,5mC}.bed  (chrom, 0-based pos, pos+1, mod, frac_wt)
set -euo pipefail

MODKIT=${MODKIT:-/fs/nexus-scratch/vgandhi/dist_modkit_v0.6.4_cd85862/modkit}
SAM=${SAM:-/fs/cbcb-software/RedHat-8-x86_64/local/samtools/1.16/bin/samtools}
BENCH=${BENCH:-/fs/cbcb-lab/storm/bds062/data/benchmark}
OUT=${OUT:-/fs/nexus-scratch/vgandhi/hp_labels}
MODS=${MODS:-"4mC 6mA 5mC"}
T=${T:-4}

# differential thresholds
MIN_COV=${MIN_COV:-10}
WT_MIN_FRAC=${WT_MIN_FRAC:-50}     # percent modified required in WT
WGA_MAX_FRAC=${WGA_MAX_FRAC:-10}   # percent allowed in WGA

mkdir -p $OUT

for MOD in $MODS; do
    for S in HP26695_WT_5kHz HP26695_WGA_5kHz; do
        BED=$OUT/${S}_${MOD}.bed
        [ -s $BED ] && { echo "reusing $BED"; continue; }
        SRC=$BENCH/bacteria/$S/modbam/${S}_sup_v5r3_${MOD}.bam
        # the share is read-only, so index via a workspace symlink
        ln -sf $SRC $OUT/${S}_${MOD}.bam
        [ -s $OUT/${S}_${MOD}.bam.bai ] || $SAM index -@ $T $OUT/${S}_${MOD}.bam
        $MODKIT pileup $OUT/${S}_${MOD}.bam $BED --threads $T 2> $OUT/${S}_${MOD}.log
        echo "$S $MOD: $(wc -l < $BED) pileup rows"
    done

    # differential: high in WT, low in WGA, covered in both.
    # bedMethyl columns: 1 chrom, 2 start, 3 end, 4 mod code, 6 strand,
    # 10 valid coverage, 11 percent modified.
    python3 - "$OUT" "$MOD" "$MIN_COV" "$WT_MIN_FRAC" "$WGA_MAX_FRAC" <<'PY'
import sys
out, mod, min_cov, wt_min, wga_max = sys.argv[1], sys.argv[2], int(sys.argv[3]), float(sys.argv[4]), float(sys.argv[5])

def load(path):
    d = {}
    with open(path) as f:
        for line in f:
            p = line.split("\t")
            cov, frac = int(p[9]), float(p[10])
            d[(p[0], int(p[1]), p[5])] = (cov, frac)
    return d

wt = load(f"{out}/HP26695_WT_5kHz_{mod}.bed")
wga = load(f"{out}/HP26695_WGA_5kHz_{mod}.bed")
n = 0
with open(f"{out}/gt_{mod}.bed", "w") as fo:
    for key, (cov_wt, frac_wt) in wt.items():
        if cov_wt < min_cov or frac_wt < wt_min:
            continue
        cov_wga, frac_wga = wga.get(key, (0, 0.0))
        if cov_wga < min_cov or frac_wga > wga_max:
            continue
        chrom, pos, strand = key
        fo.write(f"{chrom}\t{pos}\t{pos+1}\t{mod}\t{frac_wt:.1f}\t{strand}\n")
        n += 1
print(f"gt_{mod}.bed: {n} differential sites "
      f"(WT frac>={wt_min}, WGA frac<={wga_max}, cov>={min_cov} both)")
PY
done
echo "=== HP LABELS DONE: $(date) ==="
