# -*- coding: utf-8 -*-
"""step1_build.py — 重建 §3.4 使用的 2D 逐星样本（637,973 颗，含 snrr）。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv
输出：输出/step1_2d.csv（uid, subclass, active, ewha_err, snrr）
"""
import os
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
CAT = os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")
GAI = os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv")

cat = pd.read_csv(CAT, low_memory=False,
                  usecols=["uid","subclass","ewha","ewha_err","snrr","ra","dec",
                           "gaia_source_id","gaia_g_mean_mag"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50)]
cat = cat[(cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(GAI, low_memory=False, usecols=["source_id","parallax","parallax_error","ruwe"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)].copy()
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
dist = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
u["absz"] = np.abs(dist*np.sin(np.deg2rad(sc.galactic.b.deg)))
u = u[u.absz <= 3000].copy()
MG = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(1000.0/u.parallax.values, 1, None)) + 5
u = u[MG >= 6].copy()
print("  2D sample %d (expect 637,973)" % len(u))
u[["uid","subclass","active","ewha_err","snrr"]].to_csv(
    os.path.join(OUT, "step1_2d.csv"), index=False, encoding="utf-8-sig")
print("  -> 输出/step1_2d.csv")
