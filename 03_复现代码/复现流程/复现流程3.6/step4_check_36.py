# -*- coding: utf-8 -*-
"""step4_check_36.py — §3.6（含 §3.7 的 Teff 轴）逐项比对。

输入：输出/step1_ktype.csv、step1_kstats.csv、step1_teff.csv、step2_mside.csv、step2_matched.csv
输出：输出/step4_check.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")
K = pd.read_csv(os.path.join(OUT, "step1_ktype.csv"), **RD).set_index("subtype")
S = pd.read_csv(os.path.join(OUT, "step1_kstats.csv"), **RD).set_index("sample")
T = pd.read_csv(os.path.join(OUT, "step1_teff.csv"), **RD).set_index("teff_ctr")
M = pd.read_csv(os.path.join(OUT, "step2_mside.csv"), **RD).set_index("subtype")
Y = pd.read_csv(os.path.join(OUT, "step2_matched.csv"), **RD).set_index("subtype")
Z = pd.read_csv(os.path.join(OUT, "step2_matched_2d.csv"), **RD).set_index("subtype")

rows = []
def ck(item, body, got, tol=0.06, hard=0.5):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 3),
                     diff=round(float(got)-float(body), 3), status=st))
def ck_int(item, body, got):
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got)-int(body),
                     status="PASS" if int(got) == int(body) else "STALE"))

# ---- 未清理 6,111 谱：逐型比例与计数 ----
NB = {"K0": 282, "K1": 566, "K2": 72, "K3": 1317, "K4": 411, "K5": 2368, "K7": 1095}
PB = {"K0": 1.4, "K1": 1.1, "K2": 0.0, "K3": 1.4, "K4": 2.9, "K5": 1.7, "K7": 3.9}
for s in NB:
    ck_int("n %s" % s, NB[s], K.loc[s, "n_raw"])
    ck("act %% %s" % s, PB[s], K.loc[s, "pct_raw"])
ck_int("total K spectra (6,111)", 6111, int(S.loc["raw", "n"]))
ck("overall act raw (2.0)", 2.0, float(S.loc["raw", "act_pct"]))
ck("chi2 (29.8)", 29.8, float(S.loc["raw", "chi2"]), tol=0.05)
ck("dof (6)", 6.0, float(S.loc["raw", "dof"]), tol=0.0)
ck("K1-K4 (1.86)", 1.86, float(S.loc["raw", "d_k1k4"]))
ck("err K1-K4 (0.94)", 0.94, float(S.loc["raw", "se_d"]), tol=0.01)

# ---- 清理后 5,047 谱：Teff 轴与边界 ----
ck_int("cleaned spectra (5,047)", 5047, int(S.loc["clean", "n"]))
ck("cleaned overall (1.74)", 1.74, float(S.loc["clean", "act_pct"]), tol=0.01)
ck("Teff axis min (0.0)", 0.0, float(T.act_pct.min()))
ck("Teff axis max (4.34 -> 4.3)", 4.34, float(T.act_pct.max()), tol=0.01)
ck_int("boundary n at 3,950 K (392)", 392, T.loc[3950, "n"])
ck("boundary act at 3,950 K (4.34)", 4.34, float(T.loc[3950, "act_pct"]), tol=0.01)
ck_int("boundary n at 4,050 K (235)", 235, T.loc[4050, "n"])
ck("boundary act at 4,050 K (1.70)", 1.70, float(T.loc[4050, "act_pct"]), tol=0.01)
ck_int("boundary n at 4,150 K (283)", 283, T.loc[4150, "n"])
ck("boundary act at 4,150 K (2.12)", 2.12, float(T.loc[4150, "act_pct"]), tol=0.01)
ck("cleaned overall act (1.7)", 1.7, float(S.loc["clean", "act_pct"]), tol=0.05)
# S/N 6-20 匹配位移：单侧界（计算值不超过正文给的界即为 PASS）
shift_max = float(S.loc["clean", "matched_max_shift_pp"])
rows.append(dict(item="S/N matched max shift (bound 1.9 pp)", body="<=1.9",
                 computed=round(shift_max, 2), diff="",
                 status="PASS" if shift_max <= 1.9 else "STALE"))

# ---- M 侧逐谱 ----
MP = {"dM0": 7.0, "dM1": 8.1, "dM2": 13.6, "dM3": 16.4, "dM4": 20.9,
      "dM5": 58.0, "dM6": 71.3, "dM7": 69.9, "dM8": 57.3}
for s, v in MP.items():
    ck("M spectrum %s" % s, v, M.loc[s, "act_pct_spectrum"])
ck("detectable dM7 (40.2)", 40.2, float(Z.loc["dM7", "detect_pct"]))
ck("detectable dM8 (30.6)", 30.6, float(Z.loc["dM8", "detect_pct"]))
ck_int("S/N 6-20 keeps dM7 (89)", 89, Z.loc["dM7", "n_matched"])
ck_int("S/N 6-20 keeps dM8 (7)", 7, Z.loc["dM8", "n_matched"])

from scipy.stats import spearmanr
rho, pval = spearmanr(T.index.values.astype(float), T.act_pct.values.astype(float))
ck("K bins Spearman rho (-0.11)", -0.11, float(rho), tol=0.02)
ck("K bins Spearman p (0.69)", 0.69, float(pval), tol=0.05)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step4_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step4_check.csv")