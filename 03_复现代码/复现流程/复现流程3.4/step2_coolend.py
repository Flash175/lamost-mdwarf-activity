# -*- coding: utf-8 -*-
"""step2_coolend.py — §3.4 第一句：dM5-dM8 的星数、活动率、Wilson 区间与可探测比例。

输入：输出/step1_2d.csv
输出：输出/step2_coolend.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_2d.csv"), low_memory=False, encoding="utf-8-sig")

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return (100*(c-h), 100*(c+h))

rows = []
for sp in ["dM5","dM6","dM7","dM8"]:
    s = df[df.subclass == sp]
    k, n = int(s.active.sum()), len(s)
    lo, hi = wilson(k, n)
    rows.append(dict(subtype=sp, n=n, n_act=k, frac=100*k/n, lo=lo, hi=hi,
                     detect_pct=100*(3*s.ewha_err <= 0.75).mean()))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step2_coolend.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))
print("  -> 输出/step2_coolend.csv")
