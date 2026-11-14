# -*- coding: utf-8 -*-
"""step4_lamost.py — §3.3 第二段：LAMOST-RV 样本的分片比较（表 5 的 SE 口径）。

输入：输出/step1_lamost.csv
输出：输出/step4_lamost.csv、输出/step4_counts.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_lamost.csv"), low_memory=False, encoding="utf-8-sig")

def mad(x):
    x = np.asarray(x, float)
    return 1.4826*np.median(np.abs(x-np.median(x)))

rows = []
for lo in range(0, 625, 25):
    s = df[(df.absz >= lo) & (df.absz < lo+25)]
    a, i = s[s.active], s[~s.active]
    if len(a) < 50 or len(i) < 50: continue
    sa, si = a.vz.std(), i.vz.std()
    ma, mi = mad(a.vz), mad(i.vz)
    se_s = np.sqrt(sa**2/(2*(len(a)-1)) + si**2/(2*(len(i)-1)))
    se_m = np.sqrt((1.858*sa)**2/len(a) + (1.858*si)**2/len(i))
    rows.append(dict(z_lo=lo, z_hi=lo+25, n=len(s), n_act=len(a), n_inact=len(i),
                     sig_std_act=sa, sig_std_inact=si, sig_mad_act=ma, sig_mad_inact=mi,
                     d_std=sa-si, se_std=se_s, z_std=(sa-si)/se_s,
                     d_mad=ma-mi, se_mad=se_m, z_mad=(ma-mi)/se_m))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step4_lamost.csv"), index=False, encoding="utf-8-sig")
print(T[["z_lo","z_hi","n_act","n_inact","d_std","z_std","d_mad","z_mad"]].to_string(index=False))

C = pd.DataFrame([dict(
    n_stars=len(df), median_sig_rv=float(np.median(df.sig_rv)),
    n_400_600=int(((df.absz >= 400) & (df.absz < 600)).sum()),
    n_600_1000=int(((df.absz >= 600) & (df.absz < 1000)).sum()),
    n_both_rv=int(np.isfinite(df.radial_velocity).sum()))]).round(3)
C.to_csv(os.path.join(OUT, "step4_counts.csv"), index=False, encoding="utf-8-sig")
print(); print(C.to_string(index=False))
sub = T[T.z_lo >= 200]
print("  beyond 200 pc: d_std %.2f - %.2f ; z_std %+.1f - %+.1f"
      % (sub.d_std.min(), sub.d_std.max(), sub.z_std.min(), sub.z_std.max()))
print("  -> 输出/step4_lamost.csv、输出/step4_counts.csv")
