# -*- coding: utf-8 -*-
"""step1_6d_far_region.py — 重建 §2.3 的 6D 远区（|Z|>400 pc）逐星表，并核对该节数字。

输入：../../数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv
      ../../数据/原始输入/gaia/gaia_dr3_mstars.csv
输出：../../数据/6D远区_逐星_316颗.csv
      输出/step1_6d_far.json
"""
import os, json
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"),
                  low_memory=False, usecols=["uid","subclass","teff","logg","ewha","ewha_err",
                                             "snrr","ra","dec","gaia_source_id","gaia_g_mean_mag"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")

g = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                usecols=["source_id","parallax","parallax_error","ruwe","radial_velocity"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)].copy()
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()

d = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=d*un.pc, frame="icrs")
u["absz_pc"] = np.abs(d*np.sin(np.deg2rad(sc.galactic.b.deg)))
u["z_gc_pc"] = sc.transform_to("galactocentric").z.to(un.pc).value
u["dist_pc"] = d
u["MG"] = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(d, 1, None)) + 5
u = u[u.absz_pc <= 3000].copy()

rv = u.radial_velocity.values
sel6 = np.isfinite(rv) & (np.abs(rv) < 500) & (np.abs(u["z_gc_pc"].values) < 3000)
six = u[sel6].copy()
far = six[six.absz_pc > 400].copy()
print("6D sample:", len(six), " | |Z|>400 pc:", len(far))

act = far[far.active]
brt = far[far.MG < 6]
res = dict(
    n=int(len(far)), n_act=int(far.active.sum()),
    n_snr_lt3_act=int((act.snrr < 3).sum()),
    n_logg_sentinel_act=int((act.logg == -9999).sum()),
    n_logg_lt35_real_act=int(((act.logg != -9999) & (act.logg < 3.5)).sum()),
    g_med_act=float(act.gaia_g_mean_mag.median()),
    g_med_act_snr_lt3=float(act[act.snrr < 3].gaia_g_mean_mag.median()),
    g_med_act_sentinel=float(act[act.logg == -9999].gaia_g_mean_mag.median()),
    mg_med_bin=float(far.MG.median()),
    n_mg_lt6=int(len(brt)), n_mg_lt6_act=int(brt.active.sum()),
    mg_med_bright=float(brt.MG.median()),
    n_inact=int((~far.active).sum()), n_mg_lt6_inact=int((~brt.active).sum()),
    dist_med_pc=float(far.dist_pc.median()),
    n_postcut=int((far.MG >= 6).sum()), n_postcut_act=int(far[far.MG >= 6].active.sum()),
)
for k, v in res.items(): print("  %-26s %s" % (k, round(v,3) if isinstance(v,float) else v))

far.rename(columns={"gaia_g_mean_mag":"G"}).loc[:, ["uid","subclass","G","MG","snrr","logg","teff",
    "dist_pc","absz_pc","z_gc_pc","active"]].to_csv(
    os.path.join(D, "6D远区_逐星_316颗.csv"), index=False, encoding="utf-8-sig")
json.dump(res, open(os.path.join(OUT, "step1_6d_far.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=2)
print("  -> 数据/6D远区_逐星_316颗.csv, 输出/step1_6d_far.json")
print("STEP1_OK")
