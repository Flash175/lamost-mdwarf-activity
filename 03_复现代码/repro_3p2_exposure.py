# -*- coding: utf-8 -*-
"""step7_exposure.py — §3.2 末段引用的 Figure~\ref{fig:exposures}（曝光采样）。

按 uid 统计曝光组数、相邻观测时间间隔（100 天一格），并给出图注用到的三个比例。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv（designation, uid, lmjd）
输出：输出/step7_exposure.csv、输出/step7_intervals.csv、输出/step7_sampling.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据"); OUT = os.path.join(HERE, "输出")
CAT = os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")

d = pd.read_csv(CAT, low_memory=False, usecols=["designation","uid","lmjd"])
d = d.dropna(subset=["designation","lmjd"]).sort_values(["uid","lmjd"])
n_exp = d.groupby("uid").size()
ge2 = n_exp[n_exp >= 2]
print("  DR13 dM 目录：%d 星 / %d 条谱（有 lmjd）" % (len(n_exp), len(d)))

dt = np.diff(d.lmjd.to_numpy(float))[(d.uid.to_numpy()[1:] == d.uid.to_numpy()[:-1])]
dt = dt[np.isfinite(dt) & (dt > 0)]
bins = np.arange(0, 2900, 100)
h, edges = np.histogram(dt, bins=bins)
pd.DataFrame(dict(bin_lo=edges[:-1], bin_hi=edges[1:], n_pairs=h)).to_csv(
    os.path.join(OUT, "step7_intervals.csv"), index=False, encoding="utf-8-sig")

nmax = int(min(n_exp.max(), 35))
counts = np.array([int((n_exp == k).sum()) for k in range(1, nmax+1)])
cum = np.cumsum(counts[1:])/counts[1:].sum()*100
pd.DataFrame(dict(n_exposures=np.arange(1, nmax+1), n_stars=counts,
                  cum_pct_ge2=np.r_[np.nan, cum])).to_csv(
    os.path.join(OUT, "step7_exposure.csv"), index=False, encoding="utf-8-sig")

S = pd.DataFrame([dict(n_stars=int(len(n_exp)), n_spectra=int(len(d)),
                       n_ge2=int(len(ge2)),
                       pct_exactly2_of_ge2=100.0*(ge2 == 2).sum()/len(ge2),
                       pct_lt16_of_ge2=100.0*(ge2 < 16).sum()/len(ge2),
                       n_ge16=int((n_exp >= 16).sum()),
                       n_pairs=int(len(dt)), median_dt=float(np.median(dt)),
                       first_bin_pairs=int(h[0]))]).round(3)
S.to_csv(os.path.join(OUT, "step7_sampling.csv"), index=False, encoding="utf-8-sig")

print("  >=2 次曝光的 %d 颗：恰好 2 次 %.1f%%、<16 次 %.2f%%（>=16 次 %d 颗）"
      % (S.n_ge2[0], S.pct_exactly2_of_ge2[0], S.pct_lt16_of_ge2[0], S.n_ge16[0]))
print("  相邻观测间隔 %d 对，中位 %.0f 天；首个 100 天箱 %d 对"
      % (S.n_pairs[0], S.median_dt[0], S.first_bin_pairs[0]))
print("  累计（分母= >=2 次的 %d 颗）：2 次 %.1f%%、16 次 %.1f%%"
      % (int(counts[1:].sum()), cum[0], cum[16-2]))
print("  -> 输出/step7_exposure.csv、输出/step7_intervals.csv、输出/step7_sampling.csv")
