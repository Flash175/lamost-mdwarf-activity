# -*- coding: utf-8 -*-
"""step5_epoch.py — §3.2 第三段：只看单历元星（去掉"最大 EW"选择偏差）后的近/远比例。

历元数 = 同一 uid 在"EW 门 + Gaia 交叉 + 视差/ruwe 门"之后的谱数。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv、
      输出/step1_2d.csv
输出：输出/step5_epoch.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据"); OUT = os.path.join(HERE, "输出")
CAT = os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")
GAI = os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv")

cat = pd.read_csv(CAT, low_memory=False,
                  usecols=["uid","ewha","ewha_err","gaia_source_id"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50)]
cat = cat[(cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(GAI, low_memory=False, usecols=["source_id","parallax","parallax_error","ruwe"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)]
ep = m.groupby("uid").size().rename("n_ep")
print("  有历元计数的星 %d（中位历元数 %.0f）" % (len(ep), ep.median()))

df = pd.read_csv(os.path.join(OUT, "step1_2d.csv"), low_memory=False, encoding="utf-8-sig")
j = df.join(ep, on="uid")
one = j[j.n_ep == 1]
print("  2D 样本 %d 颗；单历元 %d 颗 = %.2f%%" % (len(j), len(one), 100*len(one)/len(j)))
print("  >=2 历元占 %.2f%%" % (100*(j.n_ep >= 2).mean()))

rows = []
for sp in ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]:
    a = one[(one.subclass == sp) & (one.absz < 100)]
    b = one[(one.subclass == sp) & (one.absz >= 200) & (one.absz < 400)]
    if len(a) < 30 or len(b) < 30: continue
    fa, fb = 100*a.active.mean(), 100*b.active.mean()
    rows.append(dict(subtype=sp, near_pct=fa, n_near=len(a), far_pct=fb, n_far=len(b),
                     ratio=fa/fb))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step5_epoch.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))
print("  倍数范围 %.2f - %.2f（%s）" % (T.ratio.min(), T.ratio.max(),
      "、".join(T.subtype)))
print("  -> 输出/step5_epoch.csv")
