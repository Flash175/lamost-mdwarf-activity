# -*- coding: utf-8 -*-
"""step5_crosscheck.py — §3.3 第二段：同时有两种 RV 的星上做交叉核对。

对同一批星分别用 Gaia RV 与 LAMOST RV 算 v_z，比较每个 25 pc 片的 Δσ；
并给出两种 RV 相差 >50 km/s 的比例（正文 2.2 per cent）。

输入：输出/step1_lamost.csv
输出：输出/step5_crosscheck.csv、输出/step5_summary.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_lamost.csv"), low_memory=False, encoding="utf-8-sig")
b = df[np.isfinite(df.radial_velocity) & (df.radial_velocity.abs() < 500)
       & np.isfinite(df.vz_gaia)].copy()
print("  both velocities (|RV|<500 on both) %d (body 81,719)" % len(b))

def mad(x):
    x = np.asarray(x, float)
    return 1.4826*np.median(np.abs(x-np.median(x)))

rows = []
for lo in range(0, 300, 25):
    s = b[(b.absz >= lo) & (b.absz < lo+25)]
    a, i = s[s.active], s[~s.active]
    if len(a) < 50 or len(i) < 50: continue
    rows.append(dict(z_lo=lo, n=len(s), n_act=len(a), n_inact=len(i),
                     d_std_gaia=a.vz_gaia.std()-i.vz_gaia.std(),
                     d_mad_gaia=mad(a.vz_gaia)-mad(i.vz_gaia),
                     d_std_lam=a.vz.std()-i.vz.std(),
                     d_mad_lam=mad(a.vz)-mad(i.vz)))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step5_crosscheck.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))

dv = (b.rv_lam - b.radial_velocity).abs()
S = pd.DataFrame([dict(n_both=len(b), pct_gt50=100*float((dv > 50).mean()),
                       median_abs_drv=float(dv.median()))]).round(3)
S.to_csv(os.path.join(OUT, "step5_summary.csv"), index=False, encoding="utf-8-sig")
print(); print(S.to_string(index=False))
print("  -> 输出/step5_crosscheck.csv、输出/step5_summary.csv")
