# -*- coding: utf-8 -*-
"""step2_check_37.py — §3.7（清理后的 Teff 轴）逐项比对，38 项。

输入：输出/step1_teffaxis.csv、输出/step1_mclean.csv
输出：输出/step2_check.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")
TA = pd.read_csv(os.path.join(OUT, "step1_teffaxis.csv"), **RD).set_index("teff_ctr")
MC = pd.read_csv(os.path.join(OUT, "step1_mclean.csv"), **RD).iloc[0]

rows = []
def ck(item, body, got, tol=0.06, hard=0.5):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 3),
                     diff=round(float(got)-float(body), 3), status=st))
def ck_int(item, body, got):
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got)-int(body),
                     status="PASS" if int(got) == int(body) else "STALE"))

# ---- 清理样本链 ----
ck_int("M dedup 759,607", 759607, MC.n_dup)
ck("M dedup active (14.35)", 14.35, float(MC.pct_dup), tol=0.01)
ck_int("M cleaned 642,596", 642596, MC.n_clean)
ck("M cleaned active (13.53)", 13.53, float(MC.pct_clean), tol=0.01)
ck_int("K cleaned spectra in bins (5,047)", 5047, MC.k_total)
ck_int("K bins between 3,950 and 5,450 (16)", 16, MC.n_bins_K)
ck_int("M bins between 3,150 and 4,150 (11)", 11, MC.n_bins_M)
ck_int("M stars inside bins (633,213)", 633213, MC.n_M_in_bins)
ck_int("M stars below 3,100 K (9,383)", 9383, MC.n_M_below3100)
ck_int("M stars above 4,150 K (0)", 0, MC.n_M_above4150)

ck("K cleaned active (1.74)", 1.74, float(MC.k_pct), tol=0.01)
AB = TA[TA.index >= 4250]
ck("above 4,250 K min (0.00)", 0.0, float(AB.pct_K.min()))
ck("above 4,250 K max (3.45)", 3.45, float(AB.pct_K.max()), tol=0.01)
ck_int("above 4,250 K smallest bin (32)", 32, int(AB.n_K.min()))
ck_int("above 4,250 K largest bin (587)", 587, int(AB.n_K.max()))

# ---- M 侧 11 个箱 ----
MB = {3150: (11867, 70.25), 3250: (32033, 39.73), 3350: (56156, 16.75), 3450: (69341, 17.30),
      3550: (87852, 15.13), 3650: (111735, 10.12), 3750: (106608, 6.62), 3850: (94986, 6.56),
      3950: (53609, 5.84), 4050: (8800, 7.72), 4150: (226, 10.18)}
for c_, (n_, f_) in MB.items():
    ck_int("M bin n at %d K" % c_, n_, TA.loc[c_, "n_M"])
    ck("M bin act %% at %d K" % c_, f_, float(TA.loc[c_, "pct_M"]), tol=0.01)

# ---- K 侧三个边界箱 ----
KB = {3950: (392, 4.34), 4050: (235, 1.70), 4150: (283, 2.12)}
for c_, (n_, f_) in KB.items():
    ck_int("K bin n at %d K" % c_, n_, TA.loc[c_, "n_K"])
    ck("K bin act %% at %d K" % c_, f_, float(TA.loc[c_, "pct_K"]), tol=0.01)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step2_check.csv")
