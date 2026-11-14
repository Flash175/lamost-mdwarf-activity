# -*- coding: utf-8 -*-
"""step2_ratios.py — §3.1 的三个倍数（3.06 / 3.58 / 3.35）与 >400 pc 合并值。

输入：../../数据/新图1_数据_2D分箱.csv、../../数据/giant_cut_by_MG.csv
输出：输出/step2_ratios.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

b = pd.read_csv(os.path.join(D, "新图1_数据_2D分箱.csv"), **RD).set_index("zbin")
far = pd.read_csv(os.path.join(D, "giant_cut_by_MG.csv"), **RD).set_index("zbin")
r1 = b.loc["<100","frac"] / b.loc["200-400","frac"]
r2 = b.loc["<100","frac"] / b.loc["400-600","frac"]
r3 = b.loc["<100","frac"] / far.loc[">400","frac_all"]
out = pd.DataFrame([
    dict(quantity="<100 / 200-400", value=round(float(r1),2), note="merged decline, 3.06"),
    dict(quantity="<100 / 400-600", value=round(float(r2),2), note="minimum, 3.58"),
    dict(quantity="<100 / >400 combined", value=round(float(r3),2), note="3.35, uses 7.11 per cent"),
    dict(quantity=">400 combined fraction", value=round(float(far.loc[">400","frac_all"]),2), note="7.11 per cent"),
])
out.to_csv(os.path.join(OUT, "step2_ratios.csv"), index=False, encoding="utf-8-sig")
print(out.to_string(index=False))
print("  -> 输出/step2_ratios.csv")
