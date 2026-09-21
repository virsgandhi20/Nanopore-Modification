# Dorado column for Table 1, bacteria (2026-09-21, job 7603725)

Dorado sup v5 modbams shipped with the benchmark (`*_sup_v5r3_{4mC,5mC,5mCG,6mA}.bam`) -> `modkit pileup`
(v0.6.4, all reads, whole genome) -> `scripts/benchmark/score_sites.py` with the same coverage floor (10),
strand collapsing and ground truth as the UniMeth column. Score = modkit percent modified for the row's mod code.

| row | modbam / code | sites | positives | pos rate | AUROC | AUPRC |
|---|---|---|---|---|---|---|
| E. coli M.SssI 5mC CpG, negatives = same CpGs in the DM strain | 5mCG / m | 1,357,324 | 687,670 | 0.5066 | 0.9999 | 0.9999 |
| E. coli M.SssI 5mC CpG, negatives = DM strain, all-context model | 5mC / m | 4,632,546 | 687,714 | 0.1485 | 0.9980 | 0.9793 |
| E. coli M.SssI 5mC CpG, within sample, all-context model | 5mC / m | 2,332,281 | 687,714 | 0.2949 | 0.9952 | 0.9794 |
| E. coli WT 5mC non-CpG (Dcm CCWGG) | 5mC / m | 2,303,208 | 23,960 | 0.0104 | 0.9999 | 0.9979 |
| E. coli WT 6mA (Dam GATC) | 6mA / a | 2,266,982 | 37,974 | 0.0168 | 0.9997 | 0.9750 |
| Anabaena 6mA (GATC) | 6mA / a | 3,753,487 | 30,920 | 0.0082 | 0.9997 | 0.9728 |
| T. denticola 6mA, `tdenticola` preset | 6mA / a | 1,765,332 | 15,368 | 0.0087 | 0.4972 | 0.0083 |
| H. pylori J99 6mA, GTNNNNNNAC preset | 6mA / a | 1,003,517 | 4,689 | 0.0047 | 0.4903 | 0.0044 |
| H. pylori J99 4mC, TCNNNNNNNGC preset | 4mC / 21839 | 644,535 | 9,683 | 0.0150 | 0.5264 | 0.0155 |

Notes
- A CpG-context model reports only CpG positions and every CpG in the M.SssI sample is methylated, so that
  row takes its negatives from the same CpG positions in the unmethylated dam-/dcm- strain (design agreed
  with Bhargav on 2026-09-21).
- T. denticola and J99 6mA: Dorado reproduces UniMeth's chance-level result on the motif presets (0.497 and
  0.490 vs 0.49 and 0.50), consistent with the de novo motif check: the preset positions are not methylated in
  these samples. Bhargav's decision (2026-09-21): drop the two rows from the table for now.
- J99 4mC preset also scores at chance (0.526) with Dorado's 4mC model. One caller only, so not a conclusion,
  but it is the ground truth behind the "4mC, bacterial motifs" row of Table 1 and should be checked the same way.
- Unlike the UniMeth column (first 15,000 reads = one region at full depth), these rows use every read.
