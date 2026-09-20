#!/usr/bin/env python
"""What does UniMeth's network actually see? (login-node safe, ~1 min)

UniMeth's 5 kHz path normalizes the RAW DAC signal with the BAM's sm/sd tags,
which Dorado writes in pA, without applying the pod5 calibration (its 4 kHz
path does apply it). Print, for a handful of reads, the calibration and the
mean/sd of the model input under each of UniMeth's three formulas, so two
datasets (ours vs UniMeth's demo) can be compared on the same scale.
"""
import argparse, statistics as st
import numpy as np, pod5, pysam

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True); ap.add_argument("--pod5", required=True)
ap.add_argument("--label", required=True); ap.add_argument("--n", type=int, default=40)
ap.add_argument("--out", help="append one TSV row here")
a = ap.parse_args()

tags = {}
with pysam.AlignmentFile(a.bam, check_sq=False) as b:
    hdr = [pg for pg in b.header.to_dict().get("PG", []) if "dorado" in pg.get("PN", pg.get("ID", "")).lower()]
    for r in b:
        if r.is_secondary or r.is_supplementary or not r.has_tag("sm"): continue
        rid = r.get_tag("pi") if r.has_tag("pi") else r.query_name
        tags[rid] = (r.get_tag("sm"), r.get_tag("sd"), r.get_tag("ts"))
        if len(tags) >= a.n: break
print(f"[{a.label}] dorado in BAM header: {[(p.get('PN'), p.get('VN')) for p in hdr][:2]}")

rows = []
with pod5.Reader(a.pod5) as rd:
    for rec in rd.reads(selection=list(tags), missing_ok=True):
        sm, sd, ts = tags[str(rec.read_id)]
        o, s = rec.calibration.offset, rec.calibration.scale
        dac = rec.signal[ts:].astype(np.float64)
        x_auto = (dac - sm) / sd                          # 5khz, dorado > 0.7.1
        x_leg = (dac - (1 - sm)) / (1 / sd)               # 5khz, dorado <= 0.7.1
        x_pa = ((dac + o) * s - sm) / sd                  # 4khz path = calibrated pA
        rows.append((o, s, sm, sd, dac.mean(), dac.std(), x_auto.mean(), x_auto.std(),
                     x_leg.mean(), x_leg.std(), x_pa.mean(), x_pa.std()))
if not rows: raise SystemExit(f"[{a.label}] no reads matched between BAM and pod5")
m = [st.median(c) for c in zip(*rows)]
names = "offset scale sm sd dac_mean dac_sd auto_mean auto_sd legacy_mean legacy_sd pA_mean pA_sd".split()
print(f"[{a.label}] medians over {len(rows)} reads")
for k, v in zip(names, m): print(f"    {k:12s} {v:12.4f}")
print(f"    device guess: {'PromethION-like' if m[1] < 0.16 else 'MinION/GridION-like'} (scale {m[1]:.4f})")
if a.out:
    import os
    new = not os.path.exists(a.out)
    with open(a.out, "a") as f:
        if new: f.write("label\tn\t" + "\t".join(names) + "\n")
        f.write(f"{a.label}\t{len(rows)}\t" + "\t".join(f"{v:.4f}" for v in m) + "\n")
