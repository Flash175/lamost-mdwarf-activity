# -*- coding: utf-8 -*-
"""step1_hbeta_tracks.py — §2.6 第一条：Hβ 与 Hα 的三条轨迹（按 |Z|）。

输入：../../数据/hbeta_activity_by_zbin.csv
输出：输出/step1_hbeta_tracks.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

df = pd.read_csv(os.path.join(D, "hbeta_activity_by_zbin.csv"), **RD)
order = ["<100", "100-200", "200-400", "400-600", "600-1000", ">1000"]
df = df.set_index("zbin").loc[order].reset_index()
df["ratio_400600_over_6001000"] = df["f_cat"] / df["f_cat"].shift(0)
out = df[["zbin","n","f_cat","f_ha","f_hb"]].copy()
out.to_csv(os.path.join(OUT, "step1_hbeta_tracks.csv"), index=False, encoding="utf-8-sig")

print(out.round(3).to_string(index=False))
d = df.set_index("zbin")
for col, name in (("f_cat","catalogue Ha"), ("f_ha","local Ha"), ("f_hb","local Hb")):
    print("  %-14s 400-600/600-1000 = %.3f" % (name, d.loc["400-600", col] / d.loc["600-1000", col]))
print("  Hbeta 相对 catalogue Ha 的偏低: 近平面 %.2f pp, 最外档 %.2f pp" % (
    d.loc["<100","f_cat"] - d.loc["<100","f_hb"], d.loc[">1000","f_cat"] - d.loc[">1000","f_hb"]))
print("  -> 输出/step1_hbeta_tracks.csv")
