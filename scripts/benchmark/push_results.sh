#!/bin/bash
# Commit and push the small result tables under benchmark_results/ so they can
# be read from the laptop side without pasting terminal output.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."
B=/fs/nexus-scratch/vgandhi/unimeth_bench
mkdir -p benchmark_results/unimeth
[ -s $B/table1_unimeth.tsv ] && cp $B/table1_unimeth.tsv benchmark_results/unimeth/
[ -s $B/table1_sanity.tsv ] && cp $B/table1_sanity.tsv benchmark_results/unimeth/
git add benchmark_results
git commit -q -m "${1:-Benchmark results update} ($(date '+%Y-%m-%d %H:%M'))" || { echo "nothing new to commit"; exit 0; }
git push -q && echo "pushed: $(git log --oneline -1)"
