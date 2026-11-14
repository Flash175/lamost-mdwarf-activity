
# -*- coding: utf-8 -*-
"""step1b_distance_precision.py — 从原始星表重建逐星视差相对误差，并给出分位数。

用途：守护 §4.12 的 "a median of 2.3 per cent (20 per cent at the selection cut)"。
输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv、
      数据/金属度_清洁样本_637973.csv（只取 uid 名单）
输出：数据/表4.12_距离精度_逐星.csv、输出/step1b_distance.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); RAW = os.path.join(D, "原始输入")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

cl = pd.read_csv(os.path.join(D, "金属度_清洁样本_637973.csv"), low_memory=False,
                 encoding="utf-8-sig", usecols=["uid"])
S = set(cl.uid)
d = pd.read_csv(os.path.join(RAW, "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                encoding="utf-8-sig", usecols=["uid", "gaia_source_id"])
d = d[d.uid.isin(S)].drop_duplicates("uid")
g = pd.read_csv(os.path.join(RAW, "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                encoding="utf-8-sig", usecols=["source_id", "parallax", "parallax_error"])
m = d.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m["frac_err_pct"] = 100.0 * m.parallax_error / m.parallax
m = m.replace([np.inf, -np.inf], np.nan).dropna(subset=["frac_err_pct"])
out = m[["uid", "frac_err_pct"]].copy()
out["frac_err_pct"] = out.frac_err_pct.round(4)
out.to_csv(os.path.join(D, "表4.12_距离精度_逐星.csv"), index=False, encoding="utf-8-sig")

v = out.frac_err_pct.values
rows = [dict(item="n stars with a usable parallax", value=len(v)),
        dict(item="median relative distance error (per cent)", value=float(np.median(v))),
        dict(item="16th percentile (per cent)", value=float(np.percentile(v, 16))),
        dict(item="84th percentile (per cent)", value=float(np.percentile(v, 84))),
        dict(item="95th percentile (per cent)", value=float(np.percentile(v, 95))),
        dict(item="worst case at the plx/err>5 cut (per cent)", value=20.0),
        dict(item="fraction above 10 per cent", value=float((v > 10).mean()))]
pd.DataFrame(rows).to_csv(os.path.join(OUT, "step1b_distance.csv"), index=False, encoding="utf-8-sig")
print(pd.DataFrame(rows).to_string(index=False))
