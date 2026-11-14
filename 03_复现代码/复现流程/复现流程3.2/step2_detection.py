# -*- coding: utf-8 -*-
"""step2_detection.py — §3.2 第一段的探测能力：3sigma(EW) <= 0.75 A 的星占比，
以及各 |Z| 档的目录 S/N 中位与合并活动率（供 step4 使用）。

输入：输出/step1_2d.csv
输出：输出/step2_detection.csv、输出/step2_bins.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_2d.csv"), low_memory=False, encoding="utf-8-sig")

df["det"] = 3*df.ewha_err <= 0.75
tab = df.groupby("subclass").agg(n=("det","size"), n_det=("det","sum"),
                                 frac=("det","mean")).reset_index()
tab["detect_pct"] = 100*tab.frac
tab = tab[["subclass","n","n_det","detect_pct"]].round({"detect_pct":3})
tab.to_csv(os.path.join(OUT, "step2_detection.csv"), index=False, encoding="utf-8-sig")
print(tab.to_string(index=False))

EDGES = [0,100,200,400,600,1000,1e9]
LAB = ["<100","100-200","200-400","400-600","600-1000",">1000"]
df["zb"] = pd.cut(df.absz, EDGES, labels=LAB)
b = df.groupby("zb", observed=False).agg(n=("active","size"), n_act=("active","sum"),
                                         snr_med=("snrr","median")).reset_index()
b["frac"] = 100*b.n_act/b.n
b = b[["zb","n","n_act","frac","snr_med"]].round({"frac":3,"snr_med":2})
b.to_csv(os.path.join(OUT, "step2_bins.csv"), index=False, encoding="utf-8-sig")
print()
print(b.to_string(index=False))
print("  sum n =", int(b.n.sum()), " active =", int(b.n_act.sum()))
print("  -> 输出/step2_detection.csv、输出/step2_bins.csv")
