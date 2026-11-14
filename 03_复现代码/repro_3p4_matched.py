# -*- coding: utf-8 -*-
"""step3_matched.py — §3.4 的 S/N 匹配：dM6 与 dM7 在三档 S/N 内的活动率，及 dM8 剩余谱数。

输入：输出/step1_2d.csv
输出：输出/step3_matched.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_2d.csv"), low_memory=False, encoding="utf-8-sig")

# 与 s551 一致：前两档半开区间（>=lo, <hi），6-20 档闭区间（>=6, <=20）
CUTS = [("3-6", 3, 6, False), ("6-12", 6, 12, False), ("6-20", 6, 20, True)]
rows = []
for lab, lo, hi, closed in CUTS:
    for sp in ["dM6","dM7","dM8"]:
        s = df[df.subclass == sp]
        s = s[(s.snrr >= lo) & ((s.snrr <= hi) if closed else (s.snrr < hi))]
        rows.append(dict(cut=lab, subtype=sp, n=len(s), n_act=int(s.active.sum()),
                         frac=100*s.active.mean() if len(s) else float("nan")))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step3_matched.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))
d8 = T[T.subtype == "dM8"]
print("  dM8 剩余谱数: %s" % list(d8.n))
print("  -> 输出/step3_matched.csv")
