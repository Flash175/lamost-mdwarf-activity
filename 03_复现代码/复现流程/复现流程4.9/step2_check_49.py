
# -*- coding: utf-8 -*-
"""step2_check_49.py — §4.9 自有条目的比对（其余条目见 README 的交叉引用表）。"""
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

for lab, v in (("0-100", 0.0043), ("100-200", 0.0063), ("200-400", 0.0079)):
    ck("binomial error in ln f %s (%s)" % (lab, v), v, float(V["sig lnf %s" % lab]), tol=0.0001)
ck("unweighted a (-1.2262)", -1.2262, float(V["a unw"]), tol=0.0002)
ck("unweighted b (-0.004920 pc^-1)", -0.004920, float(V["b unw"]), tol=0.000002)
ck("unweighted H (203 pc)", 203.0, float(V["H unw"]), tol=0.5)
ck("binomial a (-1.1783)", -1.1783, float(V["a bin"]), tol=0.0002)
ck("binomial b (-0.005198 pc^-1)", -0.005198, float(V["b bin"]), tol=0.000002)
ck("binomial nominal H (192 pc)", 192.0, float(V["H bin"]), tol=0.5)
ck("binomial H error (1.5 pc)", 1.45, float(V["H se bin"]), tol=0.06)
ck("binomial stored chi2 (470)", 470.0, float(V["chi2 stored bin"]), tol=0.5)
for tag, v in (("0-100", 0.031), ("100-200", -0.114), ("200-400", 0.072)):
    ck("binomial residual %s (%s)" % (tag, v), v, float(V["resid bin %s" % tag]), tol=0.002)
for tag, v in (("0-100", 0.064), ("100-200", -0.107), ("200-400", 0.043)):
    ck("unweighted residual %s (%s)" % (tag, v), v, float(V["resid unw %s" % tag]), tol=0.002)
ck("chi2 recomputed unweighted (537)", 537.0, float(V["chi2 recomputed unw"]), tol=2.0)
ck("chi2 recomputed binomial (470 at the optimum, 466.5 at rounded coefficients)", 470.0,
   float(V["chi2 recomputed bin"]), tol=5.0)
ck("residual / error lower (7)", 7.0, min(float(V["resid over sig bin %s" % t]) for t in ("0-100", "100-200", "200-400")), tol=0.5)
ck("residual / error upper (18)", 18.0, max(float(V["resid over sig bin %s" % t]) for t in ("0-100", "100-200", "200-400")), tol=0.5)
ck("per-star b (-0.004961 pc^-1)", -0.004961, float(V["perstar b"]), tol=0.000004)
ck("per-star b error (0.000040)", 0.000040, float(V["perstar b se"]), tol=0.000002)
ck("per-star H (202 pc)", 202.0, float(V["perstar H"]), tol=0.5)
ck("per-star H error (2 pc)", 2.0, float(V["perstar H se"]), tol=0.45)

# 逐亚型阈值换算（§4.9 第 4 段）
ck("equivalent EW threshold at dM0 (0.86 A)", 0.86, float(V["dM0 ew at 1e4"]), tol=0.006)
ck("equivalent EW threshold at dM6 (5.70 A)", 5.70, float(V["dM6 ew at 1e4"]), tol=0.006)
ck("EW pass / L fail at dM0 (0.6 pp)", 0.6, float(V["dM0 gap EW-L"]), tol=0.02)
ck("EW pass / L fail at dM6 (36.5 pp)", 36.5, float(V["dM6 gap EW-L"]), tol=0.06)
ck("EW pass / L fail minimum over dM0-dM6 (0.6)", 0.6, float(V["gap min dM0-dM6"]), tol=0.02)
ck("EW pass / L fail maximum over dM0-dM6 (36.5)", 36.5, float(V["gap max dM0-dM6"]), tol=0.06)
rows.append(dict(item="EW pass / L fail increases monotonically from dM0 to dM6", body="True",
                 computed=bool(round(float(V["gap monotone dM0-dM6"]))), diff=np.nan,
                 status="PASS" if round(float(V["gap monotone dM0-dM6"])) == 1 else "STALE"))
ck("same stars: n dM0 (201,622)", 201622.0, float(V["dM0 n"]), tol=0.5)
ck("same stars: n dM6 (3,671)", 3671.0, float(V["dM6 n"]), tol=0.5)
ck("indicative dM7 L fraction (3.1 %)", 3.1, float(V["dM7 pct L"]), tol=0.06)
ck("indicative dM8 L fraction (4.7 %)", 4.7, float(V["dM8 pct L"]), tol=0.06)
ck("EW criterion dM0 (7.1 %)", 7.1, float(V["dM0 pct EW"]), tol=0.06)
ck("EW criterion dM6 (72.9 %)", 72.9, float(V["dM6 pct EW"]), tol=0.06)
ck("dropping the 3-sigma cut inflates dM0 (0.75 pp)", 0.75, float(V["L drop dM0"]), tol=0.006)
ck("dropping the 3-sigma cut inflates dM1 (0.99 pp)", 0.99, float(V["L drop dM1"]), tol=0.006)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
