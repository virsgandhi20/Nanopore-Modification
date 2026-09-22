**Title:** 5 kHz inference is silently blind on BAMs from Dorado > 0.7.1: pA-unit `sm`/`sd` tags are applied to the raw DAC signal

**Version:** UniMeth v0.3.1 from source (11215d4), `unimeth_r10.4.1_5kHz_6mA.pt` and `unimeth_r10.4.1_5kHz_5mC.pt`, R10.4.1 5 kHz, Dorado 1.4.0 and 0.9.2 with `dna_r10.4.1_e8.2_400bps_sup@v5.0.0 --emit-moves`.

**Symptom.** Inference runs without any warning, but the calls carry no information. On E. coli K-12 (Dam GATC is ~fully 6mA-methylated), mean P(6mA) is 0.0009 at GATC adenines and 0.0010 at all other adenines; per-read AUROC 0.48. The 5mC model gives P ~ 0.31 everywhere (Dcm CCWGG vs background, per-read AUROC 0.48). Same result with Dorado 0.9.2 and with `--dorado_version 0.7.1` forced.

**Cause.** In `unimeth/data/pipeline.py::get_norm_params`, the `5khz` branch for Dorado > 0.7.1 returns `shift = sm, scale = sd`, and `patch_sequence` then computes `(raw_dac - shift) / scale` on `pod5_read.signal`. Current Dorado writes `sm`/`sd` in pA (with the v5 models they are the fixed standardisation constants, `sm=93.6924`, `sd=23.5067` on every read), so the pod5 calibration `pA = (dac + offset) * scale` has to be applied first. Only the `4khz` branch does that conversion.

Median over 40 reads of the value the network receives:

| data | Dorado | `sm` / `sd` | 5khz new branch | 5khz legacy branch | 4khz branch (calibrated) |
|---|---|---|---|---|---|
| your demo (`demo.bam`, `subset_18.pod5`) | 0.7.1 | -776.84 / 0.0080 | 192746 +/- 14992 | **-0.12 +/- 0.95** | 108819 +/- 2805 |
| our E. coli | 1.4.0 | 93.69 / 23.51 | 13.27 +/- 4.26 | 11714 +/- 2353 | **-0.11 +/- 0.94** |

So the demo works because its BAM predates the tag change and goes through the legacy branch; a BAM from a current Dorado goes through the new branch and the model gets unstandardised input.

**Check.** `--frequency` is read only in `get_norm_params`, so `--frequency 4khz` with the 5 kHz checkpoints is a one-flag test of calibrated input. Same reads, same checkpoints:

| model, sites | flag | mean P at GT | mean P background | per-read AUROC |
|---|---|---|---|---|
| 6mA, Dam GATC | `--frequency 5khz` | 0.0009 | 0.0010 | 0.478 |
| 6mA, Dam GATC | `--frequency 4khz` | 0.862 | 0.011 | **0.985** |
| 5mC, Dcm CCWGG | `--frequency 5khz` | 0.313 | 0.318 | 0.476 |
| 5mC, Dcm CCWGG | `--frequency 4khz` | 0.234 | 0.005 | **0.979** |

Your demo with README flags is healthy (CpG calls: 39% below 0.1, 40% above 0.9), and goes flat if forced through the 4khz branch, as expected for the old tag convention.

**Suggested fix.** For Dorado > 0.7.1, use the calibrated conversion regardless of sampling rate:
`shift = sm / scale_dacs_to_pa - shift_dacs_to_pa; scale = sd / scale_dacs_to_pa`. A cheap guard would also help: warn when the normalised signal of the first reads is far from mean 0 / sd 1.

Two smaller things found on the way: the PyPI wheel is still v0.1.0 and ships without `configs/`, and the README uses `unimeth infer` where v0.3.1 installs `unimeth-infer`.

---
Posted 2026-09-22 as https://github.com/sekeyWang/Unimeth/issues/21 (account virsgandhi20).
