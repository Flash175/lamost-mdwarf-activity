# -*- coding: utf-8 -*-
"""step3_errors.py — §3.3 第一段的标准误、显著性与亚型分层（表 5 口径）。

SE(std) = sqrt(sa^2/(2(na-1)) + si^2/(2(ni-1)))；SE(MAD) = sqrt((1.858 sa)^2/na + (1.858 si)^2/ni)

输入：输出/step1_6d.csv
输出：输出/step3_errors.csv、输出/step3_window.csv、输出/step3_strata.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_6d.csv"), low_memory=False, encoding="utf-8-sig")

def mad(x):
    x = np.asarray(x, float)
    return 1.4826*np.median(np.abs(x-np.median(x)))

def stats(a, i):
    sa, si = a.vz.std(), i.vz.std()
    ma, mi = mad(a.vz), mad(i.vz)
    se_s = np.sqrt(sa**2/(2*(len(a)-1)) + si**2/(2*(len(i)-1)))
    se_m = np.sqrt((1.858*sa)**2/len(a) + (1.858*si)**2/len(i))
    return sa, si, ma, mi, se_s, se_m

rows = []
for lo, hi in ((0,25),(25,50),(50,75),(75,100),(100,150),(150,200),(40,80)):
    s = df[(df.zgc >= lo) & (df.zgc < hi)]
    a, i = s[s.active], s[~s.active]
    sa, si, ma, mi, se_s, se_m = stats(a, i)
    rows.append(dict(slice="%d-%d" % (lo,hi), n_act=len(a), n_inact=len(i),
                     d_std=sa-si, se_std=se_s, z_std=(sa-si)/se_s,
                     d_mad=ma-mi, se_mad=se_m, z_mad=(ma-mi)/se_m))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step3_errors.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))

s = df[(df.zgc >= 40) & (df.zgc < 80)]
a, i = s[s.active], s[~s.active]
sa, si, ma, mi, se_s, se_m = stats(a, i)
W = pd.DataFrame([dict(n_act=len(a), n_inact=len(i),
                       zmed_act=float(a.zgc.median()), zmed_inact=float(i.zgc.median()),
                       sd_act=sa, sd_inact=si, d_std=sa-si, se_std=se_s, z_std=(sa-si)/se_s,
                       mad_act=ma, mad_inact=mi, d_mad=ma-mi, se_mad=se_m, z_mad=(ma-mi)/se_m)]).round(3)
W.to_csv(os.path.join(OUT, "step3_window.csv"), index=False, encoding="utf-8-sig")
print("\n  40-80 pc window:"); print(W.to_string(index=False))

rows = []
for grp, subs in (("dM0-dM2", ["dM0","dM1","dM2"]), ("dM3-dM5", ["dM3","dM4","dM5"])):
    for lo, hi in ((0,25),(25,50),(50,75),(40,80)):
        s = df[(df.zgc >= lo) & (df.zgc < hi) & df.subclass.isin(subs)]
        a, i = s[s.active], s[~s.active]
        sa, si, ma, mi, se_s, se_m = stats(a, i)
        rows.append(dict(group=grp, slice="%d-%d" % (lo,hi), n_act=len(a), n_inact=len(i),
                         d_std=sa-si, se_std=se_s, z_std=(sa-si)/se_s, d_mad=ma-mi))
G = pd.DataFrame(rows).round(3)
G.to_csv(os.path.join(OUT, "step3_strata.csv"), index=False, encoding="utf-8-sig")
print(); print(G.to_string(index=False))
print("  -> 输出/step3_errors.csv、输出/step3_window.csv、输出/step3_strata.csv")
