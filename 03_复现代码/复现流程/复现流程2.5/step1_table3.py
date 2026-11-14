# -*- coding: utf-8 -*-
"""step1_table3.py — 复现 Table 3（各亚型 EW 分位与 0.75 Å 所处位置）。

输入：../../数据/subtype_ew_percentiles.csv
输出：输出/step1_table3.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

df = pd.read_csv(os.path.join(D, "subtype_ew_percentiles.csv"), **RD)
subs = ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]
df = df.set_index("subtype").loc[subs]
df["pct_rank"] = (100 - df["frac_gt_075"]).round().astype(int)
out = df.reset_index()[["subtype","n","ew_med","ew_p90","frac_gt_075","pct_rank"]]
out.to_csv(os.path.join(OUT, "step1_table3.csv"), index=False, encoding="utf-8-sig")
print(out.round(2).to_string(index=False))
print("  sum n = %d ; 另有 dM7–dM9 621 颗未列 (总计 637,973)" % int(df.n.sum()))
print("  -> 输出/step1_table3.csv")
