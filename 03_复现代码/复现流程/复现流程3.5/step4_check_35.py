# -*- coding: utf-8 -*-
"""step4_check_35.py — §3.5 正文与 Figure 2 图注的数字逐项比对。

状态：PASS = 复现；ROUND = 仅 bootstrap 抖动；STALE = 与本样本不符（需改正文/图注）。

输入：输出/step1_counts.csv、step2_panels.csv、step2_halves.csv、step2_summary.csv
输出：输出/step4_check.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")
cnt = pd.read_csv(os.path.join(OUT, "step1_counts.csv"), **RD).iloc[0]
P = pd.read_csv(os.path.join(OUT, "step2_panels.csv"), **RD).set_index("panel")
H = pd.read_csv(os.path.join(OUT, "step2_halves.csv"), **RD).set_index("half")
S = pd.read_csv(os.path.join(OUT, "step2_summary.csv"), **RD).iloc[0]
F = pd.read_csv(os.path.join(OUT, "step3_figure.csv"), **RD).set_index("panel")
sam = pd.read_csv(os.path.join(OUT, "step1_sample.csv"), **RD)
TV = pd.read_csv(os.path.join(OUT, "step2b_tau_variant.csv"), **RD).set_index("scheme")
SU = sam[sam.usable]

rows = []
def ck(item, body, got, tol=0.006, hard=0.02):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 3),
                     diff=round(float(got)-float(body), 3), status=st))
def ck_int(item, body, got, note=""):
    ok = int(got) == int(body)
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got)-int(body),
                     status="PASS" if ok else "STALE"))

# ---------- 样本构成 ----------
ck_int("periodic sample", 3340, cnt.n_total)
ck_int("panel dM0-dM3", 2083, cnt.n_dM0_dM3)
ck_int("panel dM4-dM6", 1245, cnt.n_dM4_dM6)
ck_int("outside either panel", 12, cnt.n_outside)
ck_int("usable astrometry", 3150, cnt.n_usable)

# ---------- |Z| 分层（正文） ----------
ck("median |z| (72 pc)", 72.0, float(S.z_median), tol=0.6)
ck("half median |z| low (33)", 33.0, float(H.loc["low","zmed"]), tol=0.6)
ck("half median |z| high (156)", 156.0, float(H.loc["high","zmed"]), tol=0.6)
ck("beta low (-0.901)", -0.901, float(H.loc["low","beta"]))
ck("err beta low (0.117)", 0.117, float(H.loc["low","err"]), tol=0.01)
ck("beta high (-0.864)", -0.864, float(H.loc["high","beta"]))
ck("err beta high (0.133)", 0.133, float(H.loc["high","err"]), tol=0.01)
ck("delta beta, high - low (+0.037)", 0.037, float(S.dbeta))
ck("err delta beta (0.177)", 0.177, float(S.err_dbeta), tol=0.01)
ck("delta beta significance (0.21 sigma)", 0.21, float(S.sigma), tol=0.02)

# ---------- Figure 2 图注（面板拟合） ----------
ck("caption beta dM0-dM3 (-0.890)", -0.890, float(P.loc["dM0-dM3","beta"]))
ck("caption err dM0-dM3 (0.093)", 0.093, float(P.loc["dM0-dM3","err"]), tol=0.01)
ck("caption beta dM4-dM6 (-0.822)", -0.822, float(P.loc["dM4-dM6","beta"]))
ck("caption err dM4-dM6 (0.178)", 0.178, float(P.loc["dM4-dM6","err"]), tol=0.01)
ck("caption beta all (-0.857)", -0.857, float(P.loc["all","beta"]))
ck("caption err all (0.083)", 0.083, float(P.loc["all","err"]), tol=0.01)
db = float(P.loc["dM4-dM6","beta"]) - float(P.loc["dM0-dM3","beta"])
ed = float(np.hypot(P.loc["dM0-dM3","err"], P.loc["dM4-dM6","err"]))
ck("caption delta beta (0.068)", 0.068, db)
ck("caption err delta beta (0.201)", 0.201, ed, tol=0.01)
ck("caption delta beta significance (0.34 sigma)", 0.34, abs(db)/ed, tol=0.02)

# ---------- Ro = P/tau_c 的构造与图点 ----------
ck("tau_c median (d)", 33.0, float(SU.tau_c.median()), tol=1.5)
ck("tau_c min (d)", 15.5, float(SU.tau_c.min()), tol=0.6)
ck("tau_c max (d)", 78.1, float(SU.tau_c.max()), tol=0.6)
ck("Ro min", 0.002, float(SU.Ro.min()), tol=0.0006, hard=0.01)
ck("Ro max", 2.405, float(SU.Ro.max()), tol=0.01, hard=0.05)
ck_int("figure panel dM0-dM3 N (2,083)", 2083, F.loc["dM0-dM3", "n_total"])
ck_int("figure panel dM4-dM6 N (1,245)", 1245, F.loc["dM4-dM6", "n_total"])
ck_int("figure unsat points dM0-dM3 (819)", 819, F.loc["dM0-dM3", "n_unsat"])
ck_int("figure unsat points dM4-dM6 (187)", 187, F.loc["dM4-dM6", "n_unsat"])
ck("figure fit beta dM0-dM3", -0.890, float(F.loc["dM0-dM3", "beta"]))
ck("figure fit beta dM4-dM6", -0.822, float(F.loc["dM4-dM6", "beta"]))

ck("95 per cent upper limit on |dbeta| (0.39)", 0.39, float(S["upper95"]), tol=0.01)
ck("2 sigma detection threshold (0.35)", 0.35, float(S["thr2"]), tol=0.006)

TVl = TV.loc["Li (-0.30)"]; TVw = TV.loc["Wright (+0.31)"]
ck("tau variant: break with +0.31 (0.069)", 0.069, float(TVw.ro_sat_argmin), tol=0.002)
ck("tau variant: break with -0.30 (0.135)", 0.135, float(TVl.ro_sat_argmin), tol=0.002)
ck("tau variant: median tau ratio (1.40)", 1.40, float(TVw.tau_med)/float(TVl.tau_med), tol=0.01)
ck("tau variant: dbeta with +0.31 (-0.22)", -0.22, float(TVw.dbeta), tol=0.01)
ck("tau variant: sigma with +0.31 (0.9)", 0.93, float(TVw.sigma), tol=0.01)
ck("tau variant: dbeta with -0.30 (+0.038)", 0.038, float(TVl.dbeta), tol=0.006)
ck("tau variant: panel dbeta with +0.31 (0.34)", 0.34, float(TVw.dbeta_panel), tol=0.006)
mx = max(float(TVl.sigma), float(TVw.sigma))
rows.append(dict(item="tau variant: both nulls below 1 sigma", body=1.0, computed=round(mx, 4),
                 diff=np.nan, status="PASS" if mx < 1.0 else "STALE"))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step4_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step4_check.csv")