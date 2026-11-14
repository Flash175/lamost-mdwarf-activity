# -*- coding: utf-8 -*-
"""step2_quantile.py — 复现“改用各亚型 80 分位判据后”的亚型差异与下降倍数。

输入：../../数据/subtype_quantile_by_Z.csv、../../数据/subtype_ew_percentiles.csv
输出：输出/step2_quantile.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

q = pd.read_csv(os.path.join(D, "subtype_quantile_by_Z.csv"), **RD).set_index("subtype")
subs = ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]
q = q.loc[subs]
q["ratio_000_to_200_400"] = q["q80_<100"] / q["q80_200-400"]
out = q.reset_index()[["subtype","thr_p80","q80_<100","q80_200-400","ratio_000_to_200_400"]]
out.to_csv(os.path.join(OUT, "step2_quantile.csv"), index=False, encoding="utf-8-sig")

p = pd.read_csv(os.path.join(D, "subtype_ew_percentiles.csv"), **RD).set_index("subtype").loc[subs]
print(out.round(2).to_string(index=False))
print("\n  固定阈值判据：<100 pc 的亚型跨度 %.1f–%.1f pp（差 %.1f）" %
      (p.frac_gt_075.min(), p.frac_gt_075.max(), p.frac_gt_075.max()-p.frac_gt_075.min()))
print("  分位判据    ：<100 pc 的亚型跨度 %.1f–%.1f pp（差 %.1f）" %
      (q["q80_<100"].min(), q["q80_<100"].max(), q["q80_<100"].max()-q["q80_<100"].min()))
print("  下降倍数    ：%.2f–%.2f（q80(<100 pc) / q80(200–400 pc)）" %
      (q["ratio_000_to_200_400"].min(), q["ratio_000_to_200_400"].max()))
print("  -> 输出/step2_quantile.csv")
