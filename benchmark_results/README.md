# Results index

Small result tables and metric files only (large outputs stay on cluster scratch). Figures are in `../figures/`.

| folder | line of work | contents |
|---|---|---|
| `unimeth/` | RawMod paper, Table 1 | UniMeth column (`TABLE1_UNIMETH_2026-09-20.md`), the diagnostic chain that found the normalization bug (`diag_*`, `probe.txt`, `sigstats.*`), the motif ground-truth check (`gt_check_*`), the draft upstream issue |
| `dorado/` | RawMod paper, Table 1 | Dorado column from the benchmark modbams (`TABLE1_DORADO_2026-09-21.md`, raw table, status) |
| `deepmod2/` | RawMod paper, Table 1 | DeepMod2 CpG row |
| `strand/` | RawMod paper, code to-do | strand orientation measurement and ground-truth strand audit |
| `typing_overnight/` | typing project | 18 per-read typing experiments, one metrics JSON each, `summary.tsv` |
| `typing_followup/` | typing project | seed replicates, normalization and window tests, motif-matched controls, RawMod embedding probe |
