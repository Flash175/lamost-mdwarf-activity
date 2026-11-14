# -*- coding: utf-8 -*-
"""step2_bins_and_band.py — 复现 Table 2 的六档 n / active% / 压线带四列。"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

MANUSCRIPT = {
    "<100":     (177128, 23.81, 1.53, 0.77, 8.22, 4.11),
    "100-200":  (170974, 12.90, 1.34, 0.67, 7.60, 3.80),
    "200-400":  (191179,  7.78, 1.13, 0.57, 6.24, 3.12),
    "400-600":  ( 71688,  6.65, 0.95, 0.47, 5.07, 2.53),
    "600-1000": ( 25773,  7.76, 1.05, 0.53, 5.01, 2.51),
    ">1000":    (  1231, 13.24, 1.30, 0.65, 6.42, 3.21),
}
ORDER = ["<100","100-200","200-400","400-600","600-1000",">1000"]
st = pd.read_csv(os.path.join(D, "金属度_清洁样本_637973.csv"), usecols=["zb","active"], **RD)
g = st.groupby("zb")["active"].agg(["count","sum"])
bd = pd.read_csv(os.path.join(D, "threshold_band_by_Z.csv"), **RD).set_index("zbin")
rows = []
for b in ORDER:
    n_m, a_m, w_m, h_m, ww_m, wi_m = MANUSCRIPT[b]
    v = dict(zbin=b, n_src=int(g.loc[b,"count"]), n_paper=n_m,
             act_src=round(100.0*g.loc[b,"sum"]/g.loc[b,"count"],2), act_paper=a_m,
             band_src=round(float(bd.loc[b,"band_070_080"]),2), band_paper=w_m,
             half_src=round(float(bd.loc[b,"bias_max_pp"]),2), half_paper=h_m,
             wide_src=round(float(bd.loc[b,"band_050_100"]),2), wide_paper=ww_m,
             wb_src=round(float(bd.loc[b,"bias_wide_pp"]),2), wb_paper=wi_m)
    v["fail"] = ";".join([k for k, s, p in (
        ("n", v["n_src"], n_m), ("act", v["act_src"], a_m), ("band", v["band_src"], w_m),
        ("half", v["half_src"], h_m), ("wide", v["wide_src"], ww_m), ("wb", v["wb_src"], wi_m)) if s != p])
    rows.append(v)
pd.DataFrame(rows).to_csv(os.path.join(OUT, "step2_table2.csv"), index=False, encoding="utf-8-sig")
print("  2D n =", int(g["count"].sum()), "| -> 输出/step2_table2.csv")
