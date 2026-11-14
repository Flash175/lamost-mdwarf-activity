# -*- coding: utf-8 -*-
"""step1_bins.py — §3.1 的六档活动率与 Wilson 区间、dM0 占比。

输入：../../数据/新图1_数据_2D分箱.csv
输出：输出/step1_bins.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

df = pd.read_csv(os.path.join(D, "新图1_数据_2D分箱.csv"), **RD)
order = ["<100", "100-200", "200-400", "400-600", "600-1000", ">1000"]
df = df.set_index("zbin").loc[order].reset_index()
df[["zbin","n","n_act","frac","lo","hi","absz_med","dm0_share"]].to_csv(
    os.path.join(OUT, "step1_bins.csv"), index=False, encoding="utf-8-sig")
print(df[["zbin","n","frac","lo","hi","dm0_share"]].round(2).to_string(index=False))
print("  sum n =", int(df.n.sum()))
print("  -> 输出/step1_bins.csv")
