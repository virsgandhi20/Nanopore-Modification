# DeepMod2, CpG 5mC row (2026-09-21, job 7604052, RTX A5000, 7 minutes)

DeepMod2 (WGLab, main branch), model `bilstm_r10.4.1_5khz_v5.0`, `--seq_type dna --file_type pod5`, Dorado 1.4.0
sup@v5.0.0 move-table BAMs aligned to E. coli K-12, first 15,000 reads of each coordinate-sorted BAM.
Positives = CpGs in E. coli DM + M.SssI; negatives = the same CpGs in the DM strain without M.SssI.
Score = DeepMod2 per-site `mod_fraction`, coverage floor 10, `scripts/benchmark/score_sites.py`.

| row | sites | positives | pos rate | AUROC | AUPRC | F1 at 0.5 |
|---|---|---|---|---|---|---|
| E. coli M.SssI vs DM, 5mC CpG | 304,368 | 154,167 | 0.5065 | 1.0000 | 1.0000 | 0.9971 |

Per-site rows: 158,503 (M.SssI) and 155,599 (DM). DeepMod2 calls CpG only, so it is N/A for every other row of Table 1.
