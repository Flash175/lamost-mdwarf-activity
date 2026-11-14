# -*- coding: utf-8 -*-
"""step2_bins.py — §3.3 第一段：50 pc 档的 sigma(vz) 与 |Z| 中位；25 pc 细片的两种估计量。

输入：输出/step1_6d.csv
输出：输出/step2_bins.csv、输出/step2_slices.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_6d.csv"), low_memory=False, encoding="utf-8-sig")

def mad(x):
    x = np.asarray(x, float)
    return 1.4826*np.median(np.abs(x-np.median(x)))

rows = []
for lo, hi in ((0,50),(50,100),(100,150),(150,200)):
    s = df[(df.zgc >= lo) & (df.zgc < hi)]
    a, i = s[s.active], s[~s.active]
    rows.append(dict(zbin="%d-%d" % (lo,hi), n_act=len(a), n_inact=len(i),
                     sd_act=a.vz.std(), sd_inact=i.vz.std(),
                     zmed_act=float(a.zgc.median()), zmed_inact=float(i.zgc.median())))
B = pd.DataFrame(rows).round(3)
B.to_csv(os.path.join(OUT, "step2_bins.csv"), index=False, encoding="utf-8-sig")
print(B.to_string(index=False))

rows = []
for lo, hi in ((0,25),(25,50),(50,75),(75,100),(100,150),(150,200)):
    s = df[(df.zgc >= lo) & (df.zgc < hi)]
    a, i = s[s.active], s[~s.active]
    rows.append(dict(zslice="%d-%d" % (lo,hi), n_act=len(a), n_inact=len(i),
                     zmed_act=float(a.zgc.median()), zmed_inact=float(i.zgc.median()),
                     d_std=a.vz.std()-i.vz.std(), d_mad=mad(a.vz)-mad(i.vz)))
S = pd.DataFrame(rows).round(3)
S.to_csv(os.path.join(OUT, "step2_slices.csv"), index=False, encoding="utf-8-sig")
print()
print(S.to_string(index=False))
w = (S.n_act + S.n_inact).astype(float)
print("  count-weighted mean d_std = %.3f (body -1.56)" % float((S.d_std*w).sum()/w.sum()))
print("  -> 输出/step2_bins.csv、输出/step2_slices.csv")
