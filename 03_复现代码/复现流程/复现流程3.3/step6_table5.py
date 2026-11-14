# -*- coding: utf-8 -*-
"""step6_table5.py — 附录 Table 5（tab:kin6d）：6D 子样本六个 |Z| 档的三分量速度弥散。

输入：输出/step1_6d.csv
输出：输出/step6_table5.csv、输出/step6_summary.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_6d.csv"), low_memory=False, encoding="utf-8-sig")

BINS = [(0,50),(50,100),(100,150),(150,200),(200,400),(400,3000)]
rows = []
for lo, hi in BINS:
    s = df[(df.zgc >= lo) & (df.zgc < hi)]
    a, i = s[s.active], s[~s.active]
    rows.append(dict(zbin="%d-%d" % (lo,hi), n_act=len(a), n_inact=len(i),
                     svx_act=a.vx.std(), svy_act=a.vy.std(), svz_act=a.vz.std(),
                     svx_inact=i.vx.std(), svy_inact=i.vy.std(), svz_inact=i.vz.std(),
                     zmed_act=float(a.zgc.median()), zmed_inact=float(i.zgc.median())))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step6_table5.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))

S = pd.DataFrame([dict(n_6d=len(df), n_beyond400=int((df.zgc >= 400).sum()),
                       pct_beyond400=100.0*(df.zgc >= 400).mean())]).round(3)
S.to_csv(os.path.join(OUT, "step6_summary.csv"), index=False, encoding="utf-8-sig")
print(); print(S.to_string(index=False))
print("  -> 输出/step6_table5.csv、输出/step6_summary.csv")