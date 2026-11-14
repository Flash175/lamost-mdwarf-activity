# -*- coding: utf-8 -*-
"""step3_decompositions.py — Table 4（六参量）与 Table 9（移位—份额）。"""
import os, json
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")
six = pd.read_csv(os.path.join(D, "多因素_六因素贡献_snrr版_降幅口径.csv"), **RD)
res = {}
for b in ("100-200", "200-400"):
    s = six[six["bin"] == b].set_index("factor")
    res[b] = dict(observed=float(s.loc["OBSERVED_DECLINE","contribution_to_decline_pp"]),
                  explained_pp=float(s.loc["EXPLAINED_SUM","contribution_to_decline_pp"]),
                  explained_pct=float(s.loc["EXPLAINED_SUM","share_of_observed_pct"]),
                  unexplained_pp=float(s.loc["UNEXPLAINED","contribution_to_decline_pp"]),
                  unexplained_pct=float(s.loc["UNEXPLAINED","share_of_observed_pct"]))
ss = pd.read_csv(os.path.join(D, "shift_share.csv"), **RD)
tr = "200-400>400-600"; sub = ss[ss["transition"] == tr]
res["shift_share_" + tr] = dict(rate=float(sub["rate_pp"].sum()), comp=float(sub["comp_pp"].sum()),
                                inter=float(sub["inter_pp"].sum()),
                                total=float(sub["rate_pp"].sum()+sub["comp_pp"].sum()+sub["inter_pp"].sum()))
json.dump(res, open(os.path.join(OUT, "step3_decomp.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=2)
print("  -> 输出/step3_decomp.json")
