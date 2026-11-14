# -*- coding: utf-8 -*-
"""step1_sample.py — §3.5 的周期样本：3,340 颗的构成与"可用天体测量"子样本（含银心 |z|）。

输入：数据/原始输入/processed/dr13/fig5_sample_strict_v4.csv
      数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv（uid -> gaia_source_id）
      数据/原始输入/gaia/gaia_dr3_mstars.csv（视差、ruwe）
输出：输出/step1_sample.csv、输出/step1_counts.csv
"""
import os
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

f5 = pd.read_csv(os.path.join(D, "原始输入", "processed", "dr13", "fig5_sample_strict_v4.csv"),
                 low_memory=False, encoding="utf-8-sig", comment="#")
spt = f5.spt.astype(str).str.strip()
f5["panel"] = np.where(spt.isin(["dM0","dM1","dM2","dM3"]), "dM0-dM3",
              np.where(spt.isin(["dM4","dM5","dM6"]), "dM4-dM6", "outside"))
cnt = f5.panel.value_counts().to_dict()
print("  periodic sample %d : dM0-dM3 %d, dM4-dM6 %d, outside %d"
      % (len(f5), cnt.get("dM0-dM3",0), cnt.get("dM4-dM6",0), cnt.get("outside",0)))

cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"),
                  low_memory=False, usecols=["uid","ra","dec","gaia_source_id"]).drop_duplicates("uid")
cat["gid"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                usecols=["source_id","parallax","parallax_error","ruwe"])
m = cat.merge(g, left_on="gid", right_on="source_id", how="inner")
m = m[np.isfinite(m.parallax) & (m.parallax > 0)].copy()
dist = 1000.0/m.parallax.values
sc = SkyCoord(ra=m.ra.values*un.deg, dec=m.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
m["z_gc"] = np.abs(sc.transform_to("galactocentric").z.to(un.pc).value)
J = f5.merge(m[["uid","parallax","parallax_error","ruwe","z_gc"]], on="uid", how="left",
             suffixes=("", "_gaia"))
J["usable"] = ((J.parallax > 0.2) & (J.parallax/J.parallax_error > 5)
               & (pd.to_numeric(J.ruwe_gaia, errors="coerce") < 1.4))
print("  usable astrometry: %d（正文 3,142）" % int(J.usable.sum()))
J["tau_c"] = pd.to_numeric(J.P, errors="coerce") / pd.to_numeric(J.Ro, errors="coerce")
J[["uid","spt","panel","teff","M","P","Ro","sRo","RHa","tau_c","z_gc","usable"]].to_csv(
    os.path.join(OUT, "step1_sample.csv"), index=False, encoding="utf-8-sig")
U = J[J.usable]
print("  Ro 范围 %.4f-%.3f（中位 %.3f）" % (U.Ro.min(), U.Ro.max(), U.Ro.median()))
print("  tau_c = P/Ro：中位 %.1f d，范围 %.1f-%.1f d" % (U.tau_c.median(), U.tau_c.min(), U.tau_c.max()))
print("  sRo/Ro 中位 %.3f ; tau_c 与 M 的 Spearman %.3f"
      % ((U.sRo/U.Ro).median(), pd.Series(U.tau_c).corr(pd.Series(U.M), method="spearman")))
pd.DataFrame([dict(n_total=len(J), n_dM0_dM3=cnt.get("dM0-dM3",0), n_dM4_dM6=cnt.get("dM4-dM6",0),
                   n_outside=cnt.get("outside",0), n_usable=int(J.usable.sum()))]).to_csv(
    os.path.join(OUT, "step1_counts.csv"), index=False, encoding="utf-8-sig")
print("  -> 输出/step1_sample.csv、输出/step1_counts.csv")