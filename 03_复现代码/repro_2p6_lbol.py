# -*- coding: utf-8 -*-
"""step3_lbol.py — §2.6 第二条：L_Hα/L_bol 判据。

输入：../../数据/chi_lbol_by_zbin.csv、../../数据/chi_lbol_by_subtype.csv
输出：输出/step3_lbol_bins.csv、输出/step3_lbol_subtype.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
os.makedirs(OUT, exist_ok=True)
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

b = pd.read_csv(os.path.join(D, "chi_lbol_by_zbin.csv"), **RD)
order = ["<100", "100-200", "200-400", "400-600", "600-1000", ">1000"]
b = b.set_index("zbin").loc[order].reset_index()
b["ratio_100_to_600_1000"] = b["pct_lB"].iloc[0] / b["pct_lB"].loc[b.zbin == "600-1000"].iloc[0]
b[["zbin","n","pct_EW","pct_lA","pct_lB"]].to_csv(os.path.join(OUT, "step3_lbol_bins.csv"),
                                                    index=False, encoding="utf-8-sig")

s = pd.read_csv(os.path.join(D, "chi_lbol_by_subtype.csv"), **RD)
s["diff_EW_minus_lB"] = s["pct_EW"] - s["pct_lB"]
s[["spt","n","pct_EW","pct_lA","pct_lB","logchi","ew_at_1e4","diff_EW_minus_lB"]].to_csv(
    os.path.join(OUT, "step3_lbol_subtype.csv"), index=False, encoding="utf-8-sig")

print(b[["zbin","pct_EW","pct_lB"]].round(2).to_string(index=False))
ew = b.set_index("zbin")
print("\n  L/L_bol 判据的降幅: %.3f  (EW 判据: %.3f)" % (
    b["pct_lB"].iloc[0] / float(b.loc[b.zbin == "600-1000","pct_lB"].iloc[0]),
    float(ew.loc["<100","pct_EW"]) / float(ew.loc["600-1000","pct_EW"])))
print("  EW 等价 1e-4: dM0 %.2f A -> dM7 %.2f A" % (
    s.loc[s.spt=="dM0","ew_at_1e4"].iloc[0], s.loc[s.spt=="dM7","ew_at_1e4"].iloc[0]))
print("  两判据差: dM0 %.2f pp -> dM6 %.2f pp" % (
    s.loc[s.spt=="dM0","diff_EW_minus_lB"].iloc[0], s.loc[s.spt=="dM6","diff_EW_minus_lB"].iloc[0]))
print("  -> 输出/step3_lbol_bins.csv, step3_lbol_subtype.csv")
