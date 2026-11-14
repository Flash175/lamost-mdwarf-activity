
# -*- coding: utf-8 -*-
"""step2_check_48.py — §4.8 的逐项比对。"""
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

# 原始 6,111 条
ck_int("raw pilot spectra (6,111)", 6111, V["raw n"])
ck("raw chi-square (29.8)", 29.8, float(V["raw chi2"]), tol=0.1)
ck_int("raw dof (6)", 6, V["raw dof"])
ck("raw p value (4e-5)", 4.0, float(V["raw p"])*1e5, tol=0.6)
ck("K7 raw (3.9 %)", 3.9, float(V["raw pct K7"]), tol=0.06)
ck("K4 raw (2.9 %)", 2.9, float(V["raw pct K4"]), tol=0.06)
ck("K1 raw (1.1 %)", 1.1, float(V["raw pct K1"]), tol=0.06)

# 清理后 5,047 条
ck_int("cleaned pilot spectra (5,047)", 5047, V["clean n"])
ck_int("cleaned Teff bins (16)", 16, V["clean bins"])
ck_int("cleaned bins cover all spectra (5,047)", 5047, V["clean bins n"])
ck("cleaned max bin fraction (4.3 %)", 4.3, float(V["clean pct max"]), tol=0.06)
ck("cleaned min bin fraction (0.0 %)", 0.0, float(V["clean pct min"]), tol=0.01)
ck("no monotonic trend: Spearman rho (-0.11)", -0.11, float(V["clean spearman rho"]), tol=0.02)
ck("no monotonic trend: Spearman p (0.69)", 0.69, float(V["clean spearman p"]), tol=0.05)

# EW 分位（原始试点）
ck("median EW K0 (-1.01 A)", -1.01, float(V["ew med K0"]), tol=0.006)
ck("median EW K7 (-0.27 A)", -0.27, float(V["ew med K7"]), tol=0.006)
rows.append(dict(item="median EW negative for every subtype", body="<0", computed=round(float(V["ew med max (all negative?)"]), 4),
                 diff=np.nan, status="PASS" if float(V["ew med max (all negative?)"]) < 0 else "STALE"))
ck("90th percentile below +0.09 A", 0.09, float(V["ew p90 max"]), tol=0.006)
ck("0.75 A criterion passes 0 per cent (minimum)", 0.0, float(V["frac075 min"]), tol=0.01)
ck("0.75 A criterion passes 3.9 per cent (maximum)", 3.9, float(V["frac075 max"]), tol=0.06)

# 三段分解
for s, a, b, c in (("dM5", 58.0, 57.8, 56.8), ("dM7", 69.9, 70.9, 73.3)):
    ck("%s per spectrum (%s)" % (s, a), a, float(V["%s spec" % s]), tol=0.06)
    ck("%s per star (%s)" % (s, b), b, float(V["%s star" % s]), tol=0.06)
    ck("%s per star with 2D cuts (%s)" % (s, c), c, float(V["%s 2d" % s]), tol=0.06)

# 试点选择范围
ck("selection G minimum (14)", 14.0, float(V["sel G min"]), tol=0.05)
ck("selection G maximum (18)", 18.0, float(V["sel G max"]), tol=0.05)
ck("selection S/N maximum (30)", 30.0, float(V["sel snr max"]), tol=0.05)
ck("selection Teff minimum (3,900 K)", 3900.0, float(V["sel teff min"]), tol=2.0)
ck("selection Teff maximum (5,500 K)", 5500.0, float(V["sel teff max"]), tol=10.0)
rows.append(dict(item="selection log g above 3.8", body=">3.8", computed=round(float(V["sel logg min"]), 3),
                 diff=np.nan, status="PASS" if float(V["sel logg min"]) > 3.8 else "STALE"))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
