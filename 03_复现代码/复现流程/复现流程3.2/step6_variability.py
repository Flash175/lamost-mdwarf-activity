# -*- coding: utf-8 -*-
"""step6_variability.py — §3.2 末段引用的 Figure~\ref{fig:variability}（历元间 Hα 变幅）。

逐星统计 EW 的峰峰值（max-min），并定位图注里的三颗示例星：
5 / 8 / 22 个历元，峰峰值 4.95 / 4.95 / 4.88 A（容差 0.02 A）。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv
输出：输出/step6_examples.csv、输出/step6_series.csv、输出/step6_variability.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据"); OUT = os.path.join(HERE, "输出")
CAT = os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")

d = pd.read_csv(CAT, low_memory=False,
                usecols=["uid","subclass","ewha","ewha_err","lmjd"])
w = d[(d.ewha > -20) & (d.ewha <= 50) & (d.ewha_err > 0) & (d.ewha_err <= 5)].copy()
g = w.groupby("uid").ewha.agg(n_ep="size", ew_max="max", ew_min="min")
g["pp"] = g.ew_max - g.ew_min
print("  EW 门内星 %d；多历元（>=2）%d 颗；峰峰值中位 %.2f A"
      % (len(g), int((g.n_ep >= 2).sum()), float(g.loc[g.n_ep >= 2, "pp"].median())))

TARGET = [(5, 4.95), (8, 4.95), (22, 4.88)]
rows, picks = [], []
for ne, pp in TARGET:
    m = g[(g.n_ep == ne) & (np.abs(g.pp - pp) < 0.02)]
    cand = sorted(m.index.tolist())
    rows.append(dict(n_ep_target=ne, pp_target=pp, n_candidates=len(cand), uid=cand[0]))
    picks.append(cand[0])
    print("  目标 n_ep=%2d 峰峰=%.2f A -> 命中 %d 颗（选用 %s）" % (ne, pp, len(cand), cand[0]))
E = pd.DataFrame(rows)
E.to_csv(os.path.join(OUT, "step6_examples.csv"), index=False, encoding="utf-8-sig")

ser = w[w.uid.isin(picks)][["uid","subclass","lmjd","ewha","ewha_err"]].sort_values(["uid","lmjd"])
ser.to_csv(os.path.join(OUT, "step6_series.csv"), index=False, encoding="utf-8-sig")
print(ser.groupby("uid").agg(subclass=("subclass","first"), n_ep=("ewha","size"),
                             ew_min=("ewha","min"), ew_max=("ewha","max")).round(3).to_string())

S = pd.DataFrame([dict(n_stars_with_ew=int(len(g)),
                       n_multi_epoch=int((g.n_ep >= 2).sum()),
                       frac_pp_ge1A=100*float((g.pp >= 1).mean()),
                       frac_pp_ge3A=100*float((g.pp >= 3).mean()),
                       median_pp_multi=float(g.loc[g.n_ep >= 2, "pp"].median()),
                       median_pp_all=float(g.pp.median()))]).round(3)
S.to_csv(os.path.join(OUT, "step6_variability.csv"), index=False, encoding="utf-8-sig")
print("  峰峰 >=1 A 占 %.2f%%，>=3 A 占 %.2f%%（全部 EW 门内星）"
      % (float(S.frac_pp_ge1A.iloc[0]), float(S.frac_pp_ge3A.iloc[0])))
print("  -> 输出/step6_examples.csv、输出/step6_series.csv、输出/step6_variability.csv")
