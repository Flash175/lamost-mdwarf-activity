# -*- coding: utf-8 -*-
"""step2_mside.py — §3.6 的 M 侧：逐谱活动率、逐星可探测比例与 S/N 6-20 计数。

两个逐星样本都给出（正文在该段分别用到）：
  · 去巨星前 639,711（视差 SNR>5、ruwe<1.4、每星一条最大 EW 谱）
  · 去巨星后 637,973（再加 M_G>=6 与 |Z|<=3 kpc；可探测比例 40.2/30.6 与 89/7 用的就是它）

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv
输出：输出/step2_mside.csv、step2_matched.csv（639,711）、step2_matched_2d.csv（637,973）
"""
import os
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出")

cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                  usecols=["uid","subclass","ewha","ewha_err","snrr","gaia_source_id",
                           "gaia_g_mean_mag","ra","dec"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)

rows = []
for s in ["dM0","dM1","dM2","dM3","dM4","dM5","dM6","dM7","dM8"]:
    x = cat[cat.subclass == s]
    rows.append(dict(subtype=s, n_spectra=len(x), act_pct_spectrum=100*x.active.mean(),
                     detect_pct_spectrum=100*(3*x.ewha_err <= 0.75).mean()))
M = pd.DataFrame(rows).round(3)
M.to_csv(os.path.join(OUT, "step2_mside.csv"), index=False, encoding="utf-8-sig")
print(M.to_string(index=False))

cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                usecols=["source_id","parallax","parallax_error","ruwe"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)]
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
print("\n  去巨星前逐星样本 %d（正文 639,711）" % len(u))

def tally(df, label, cols=("dM6","dM7","dM8")):
    out = []
    for s in cols:
        x = df[df.subclass == s]; mm = x[(x.snrr >= 6) & (x.snrr <= 20)]
        out.append(dict(subtype=s, n_stars=len(x), n_matched=len(mm),
                        act_pct=100*x.active.mean(),
                        act_pct_matched=100*mm.active.mean() if len(mm) else np.nan,
                        detect_pct=100*(3*x.ewha_err <= 0.75).mean()))
    T = pd.DataFrame(out).round(3)
    print("\n  %s" % label); print(T.to_string(index=False))
    return T

Y = tally(u, "639,711（去巨星前）")
Y.to_csv(os.path.join(OUT, "step2_matched.csv"), index=False, encoding="utf-8-sig")

# 去巨星后 637,973
dist = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
u["absz"] = np.abs(dist*np.sin(np.deg2rad(sc.galactic.b.deg)))
MG = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(1000.0/u.parallax.values, 1, None)) + 5
u2 = u[(MG >= 6) & (u.absz <= 3000)].copy()
print("\n  去巨星后逐星样本 %d（正文 637,973）" % len(u2))
Z = tally(u2, "637,973（去巨星 + |Z|<=3 kpc）", cols=("dM7","dM8"))
Z.to_csv(os.path.join(OUT, "step2_matched_2d.csv"), index=False, encoding="utf-8-sig")
print("  -> 输出/step2_mside.csv、step2_matched.csv、step2_matched_2d.csv")
