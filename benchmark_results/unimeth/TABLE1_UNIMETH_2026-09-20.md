# UniMeth column for Table 1, bacteria (2026-09-20)

UniMeth v0.3.1 from source, `unimeth_r10.4.1_5kHz_{6mA,5mC}.pt`, Dorado 1.4.0 `sup@v5.0.0 --emit-moves`,
run with `--frequency 4khz` (see `UPSTREAM_ISSUE_DRAFT.md`: the 5 kHz path applies pA-unit `sm`/`sd` tags to
raw DAC signal and is blind on current-Dorado BAMs). First 15,000 reads of each coordinate-sorted BAM (a
contiguous region at full depth), coverage floor 10, strands collapsed per position, negatives = every other
covered A (6mA) or C (5mC) position. Raw per-read calls are kept under
`/fs/nexus-scratch/vgandhi/unimeth_bench/<sample>.<mod>/calls.txt`.

| dataset | GT | sites | positives | AUROC (call freq) | AUPRC | AUROC (mean P) | AUPRC | per-read AUROC |
|---|---|---|---|---|---|---|---|---|
| E. coli WT 6mA | Dam GATC motif | 681,286 | 11,439 | 0.9998 | 0.978 | 0.9998 | 0.977 | 0.985 |
| Anabaena 6mA | GATC motif | 1,244,766 | 10,336 | 0.9998 | 0.982 | 0.9999 | 0.984 | 0.977 |
| H. pylori 26695 6mA | WT vs WGA differential | 379,814 | 18,726 | 0.9988 | 0.965 | 0.9989 | 0.967 | 0.967 |
| E. coli M.SssI 5mC CpG | CG motif | 549,098 | 163,704 | 1.0000 | 1.000 | 1.0000 | 1.000 | 0.997 |
| E. coli WT 5mC non-CpG | Dcm CCWGG motif | 465,699 | 4,738 | 0.7244 | 0.212 | 0.9992 | 0.914 | 0.979 |
| H. pylori 26695 5mC | WT vs WGA differential | 237,194 | 6,516 | 0.8460 | 0.662 | 0.9290 | 0.692 | 0.915 |
| T. denticola 6mA | `tdenticola` preset | 1,115,505 | 9,768 | 0.4905 | 0.008 | 0.5540 | 0.009 | 0.621 |
| H. pylori J99 6mA | GTNNNNNNAC, A on both strands | 1,004,601 | 4,699 | 0.4993 | 0.005 | 0.5224 | 0.005 | 0.472 |

Call frequency = UniMeth's own site output (fraction of reads with P > 0.5). Mean P = prob_1_sum / coverage.
They differ only where the model is under-confident (Dcm: mean P 0.24 at true sites, so most true sites have
call frequency 0 and tie with the negatives).

## The last two rows measure the ground truth, not the caller

De novo motif enrichment (`scripts/benchmark/motif_enrichment.py`) on UniMeth's calls and, independently, on
Dorado's 6mA calls from the benchmark modbams (`*_sup_v5r3_6mA.bam` via modkit pileup). Methylated = site
frequency >= 0.7 at coverage >= 10.

| organism | preset sites methylated, UniMeth | preset sites methylated, Dorado | adenines methylated genome-wide | motifs both callers find (modified base in brackets) |
|---|---|---|---|---|
| E. coli (control) | 10,642 / 11,422 (93.2%) | 37,516 / 37,974 (98.8%) | 1.6% / 1.7% | G[A]TC, 97% / 96% of methylated sites, specificity 100% |
| Anabaena | 6,261 / 10,323 (60.7%) | not run | 0.5% | G[A]TC, 99.5% of methylated sites |
| T. denticola | 1 / 9,766 (0.0%) | 3 / 15,368 (0.0%) | 1.3% / 1.6% | AA[A]TTT, AA[A]TTC, GA[A]TTT, GA[A]TTC (= RA[A]TTY, 77% / 72%), CTA[A]T (9.5% / 8.5%), GAAG[A]G (4.2% / 5.4%) |
| H. pylori J99 | 114 / 4,699 (2.4%) | 127 / 4,689 (2.7%) | 4.7% / 5.3% | C[A]TG (27% / 26%), G[A]TC (21% / 20%), G[A]GG (9% / 9%), G[A]NTC, GCCT[A], GTC[A]T, GAC[A]T, GAC[A]C, ATTA[A]T |

Two unrelated callers return the same motifs in the same proportions and agree that the preset sites of
T. denticola (GATC, TATAC) and H. pylori J99 (GTNNNNNNAC) are not methylated in these samples. Any tool scored
against those two presets is being scored on unmethylated positions.
