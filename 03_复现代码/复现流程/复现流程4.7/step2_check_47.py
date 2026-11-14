
# -*- coding: utf-8 -*-
"""step2_check_47.py — §4.7 的逐项比对。"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
V = pd.read_csv(os.path.join(OUT, "step1_values.csv"), encoding="utf-8-sig").set_index("item")["value"]
rows = []
def ck(item, body, got, tol=0.006, hard=0.05):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 4),
                     diff=round(float(got) - float(body), 4), status=st))
def ck_int(item, body, got, tol=0):
    ok = abs(int(got) - int(body)) <= tol
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got) - int(body),
                     status="PASS" if ok else "STALE"))

# |Z| 两半
ck("median |z| low (33 pc)", 33.0, float(V["zmed low"]), tol=1.0)
ck("median |z| high (156 pc)", 156.0, float(V["zmed high"]), tol=1.0)
ck("beta low (-0.901)", -0.901, float(V["beta low"]), tol=0.006)
ck("err beta low (0.117)", 0.117, float(V["err low"]), tol=0.006)
ck("beta high (-0.864)", -0.864, float(V["beta high"]), tol=0.006)
ck("err beta high (0.133)", 0.133, float(V["err high"]), tol=0.006)
ck("delta beta high - low (+0.037)", 0.037, float(V["dbeta halves"]), tol=0.006)
ck("err delta beta (0.177)", 0.177, float(V["err dbeta"]), tol=0.01)
ck("significance (0.21 sigma)", 0.21, float(V["sigma halves"]), tol=0.02)
ck("two-sigma bound (0.39)", 0.39, float(V["two sigma bound"]), tol=0.006)
ck("95 per cent upper limit (0.38)", 0.38, float(V["upper95 bound"]), tol=0.006)

# 断点扫描
ck("argmin Ro_sat (0.132)", 0.132, float(V["ro argmin"]), tol=0.002)
ck("1 sigma range lower (0.092)", 0.092, float(V["range 1sig lo"]), tol=0.002)
ck("1 sigma range upper, contiguous (0.154)", 0.154, float(V["1sig contiguous hi"]), tol=0.002)
ck("isolated grid point at 1 sigma (0.160)", 0.160, float(V["1sig isolated hi"]), tol=0.002)
ck("1 per cent range lower (0.032)", 0.032, float(V["range 1pct lo"]), tol=0.002)
ck("1 per cent range upper (0.252)", 0.252, float(V["range 1pct hi"]), tol=0.002)
ck("0.5 per cent range lower (0.065)", 0.065, float(V["range 0.5pct lo"]), tol=0.002)
ck("0.5 per cent range upper (0.212)", 0.212, float(V["range 0.5pct hi"]), tol=0.002)
ck("gap depth in DeltaRSS (1.24 archival, 1.21 on this grid)", 1.241, float(V["gap max drss"]), tol=0.05)
ck("gap deepest break (0.1577)", 0.1577, float(V["gap argmax ro"]), tol=0.001)
ck_int("gap stars in (0.155, 0.160] (29)", 29, V["gap n stars"])
ck_int("gap stars dM0 (9)", 9, V["gap n dM0"])
ck_int("gap stars dM1 (9)", 9, V["gap n dM1"])
ck_int("gap stars dM2 (5)", 5, V["gap n dM2"])
ck_int("gap stars dM4 (6)", 6, V["gap n dM4"])

# 面板与四估计量
ck("beta M0-M3 (-0.890)", -0.890, float(V["beta M0M3"]), tol=0.006)
ck("err M0-M3 (0.093)", 0.093, float(V["err M0M3"]), tol=0.006)
ck("beta M4-M6 (-0.822)", -0.822, float(V["beta M4M6"]), tol=0.006)
ck("err M4-M6 (0.178)", 0.178, float(V["err M4M6"]), tol=0.006)
ck("panel dbeta (0.068)", 0.068, float(V["2x2 dbeta fi_fix"]), tol=0.006)
ck("panel dbeta significance (0.34 sigma)", 0.34, float(V["2x2 sig fi_fix"]), tol=0.02)
ck("2x2 free intercept, free break (+0.05)", 0.05, float(V["2x2 dbeta fi_free"]), tol=0.01)
ck("2x2 free intercept, free break (0.2 sigma)", 0.2, float(V["2x2 sig fi_free"]), tol=0.05)
ck("2x2 anchored, fixed break (-0.588)", -0.588, float(V["2x2 dbeta ac_fix"]), tol=0.006)
ck("2x2 anchored, fixed break (4.1 sigma)", 4.1, float(V["2x2 sig ac_fix"]), tol=0.05)
ck("2x2 anchored, free break (+0.267)", 0.267, float(V["2x2 dbeta ac_free"]), tol=0.006)
ck("2x2 anchored, free break (1.3 sigma)", 1.3, float(V["2x2 sig ac_free"]), tol=0.05)

# tau_c 约定
ck("tau_c ratio Wright/Li (1.4)", 1.4, float(V["tau ratio"]), tol=0.02)
ck("Ro_sat under +0.31 (0.069)", 0.069, float(V["tau ro sat wright"]), tol=0.002)
ck("halves dbeta under +0.31 (-0.22)", -0.22, float(V["tau dbeta wright"]), tol=0.01)
ck("halves significance under +0.31 (0.9 sigma)", 0.93, float(V["tau sigma wright"]), tol=0.02)
ck("panel dbeta under +0.31 (0.34)", 0.34, float(V["tau panel wright"]), tol=0.006)
ck("halves dbeta under -0.30 (+0.038)", 0.038, float(V["tau dbeta li"]), tol=0.006)

# 适用域
ck("periodic subsample below 0.75 A (21.8 %)", 21.8, float(V["frac below 0.75"]), tol=0.06)
ck("EW median offset vs parent (2.5 A)", 2.5, float(V["ew offset"]), tol=0.05)
ck("G offset vs parent (-1.4 mag, i.e. brighter)", -1.4, float(V["g offset"]), tol=0.05)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
