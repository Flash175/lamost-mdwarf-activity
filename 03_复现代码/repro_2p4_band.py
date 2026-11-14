# -*- coding: utf-8 -*-
"""step1_band_fractions.py — §2.4 压线带比例与影响上限（Table 2 的四列）。

输入：../../数据/threshold_band_by_Z.csv、../../数据/ew_threshold_band.csv
输出：输出/step1_band.csv（逐档比例、上限、差分）
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

bd = pd.read_csv(os.path.join(D, "threshold_band_by_Z.csv"), **RD)
glob = pd.read_csv(os.path.join(D, "ew_threshold_band.csv"), **RD)
g = glob.set_index("band")

order = ["<100","100-200","200-400","400-600","600-1000",">1000"]
bd = bd.set_index("zbin").loc[order]
trend = ["200-400","400-600","600-1000"]          # 参与趋势的三档

out = pd.DataFrame({
    "zbin": order,
    "n": bd["n"].astype(int).values,
    "band_070_080_pct": bd["band_070_080"].round(3).values,
    "halfkill_pp": bd["bias_max_pp"].round(3).values,
    "band_050_100_pct": bd["band_050_100"].round(3).values,
    "wide_pp": bd["bias_wide_pp"].round(3).values,
})
out.to_csv(os.path.join(OUT, "step1_band.csv"), index=False, encoding="utf-8-sig")

print("  全局 0.70-0.80 A 带宽（目录级）: %.2f %%" % float(g.loc["0.7-0.8","pct"]))
print("  <100 pc 带宽            : %.2f %%  → 上限 %.3f pp" % (bd.loc["<100","band_070_080"], bd.loc["<100","bias_max_pp"]))
print("  200-1000 pc 三档范围     : %.2f - %.2f %%" % (bd.loc[trend,"band_070_080"].min(), bd.loc[trend,"band_070_080"].max()))
h = bd["bias_max_pp"]; w = bd["bias_wide_pp"]
print("  单档最大上限            : %.2f pp" % h.max())
print("  逐档差分（含最内档）    : %.2f pp" % (h.max()-h.min()))
print("  三档内差分              : %.2f pp" % (h.loc[trend].max()-h.loc[trend].min()))
print("  宽带差分（含最内档）    : %.2f pp" % (w.max()-w.min()))
print("  -> 输出/step1_band.csv")
