# -*- coding: utf-8 -*-
"""step4_check_31.py — 把 §3.1 的数字逐项与计算值比对。"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
b = pd.read_csv(os.path.join(OUT, "step1_bins.csv")).set_index("zbin")
r = pd.read_csv(os.path.join(OUT, "step2_ratios.csv")).set_index("quantity")
s = pd.read_csv(os.path.join(OUT, "step3_subtypes.csv")).set_index("subtype")

rows = []
def ck(item, exp, got, tol=0.006):
    ok = abs(float(got) - float(exp)) <= tol
    rows.append(dict(item=item, expected=exp, computed=got, status="PASS" if ok else "FAIL"))

order = ["<100","100-200","200-400","400-600","600-1000",">1000"]
EXP = {
    "<100":     (177128, 23.81, 23.61, 24.01, 19.8),
    "100-200":  (170974, 12.90, 12.74, 13.06, 25.0),
    "200-400":  (191179,  7.78,  7.66,  7.90, 34.3),
    "400-600":  ( 71688,  6.65,  6.47,  6.84, 54.0),
    "600-1000": ( 25773,  7.76,  7.44,  8.10, 71.8),
    ">1000":    (  1231, 13.24, 11.46, 15.25, 81.1),
}
for k in order:
    n, f, lo, hi, d = EXP[k]
    ck("%s n" % k, n, int(b.loc[k,"n"]), 0)
    ck("%s active" % k, f, float(b.loc[k,"frac"]), 0.006)
    ck("%s Wilson lo" % k, lo, float(b.loc[k,"lo"]), 0.006)
    ck("%s Wilson hi" % k, hi, float(b.loc[k,"hi"]), 0.006)
    ck("%s dM0 share" % k, d, float(b.loc[k,"dm0_share"]), 0.05)

ck("sum n", 637973, int(b.n.sum()), 0)
ck("merged decline 3.06", 3.06, float(r.loc["<100 / 200-400","value"]), 0.006)
ck("minimum factor 3.58", 3.58, float(r.loc["<100 / 400-600","value"]), 0.006)
ck(">400 combined 3.35", 3.35, float(r.loc["<100 / >400 combined","value"]), 0.006)
ck(">400 combined fraction 7.11", 7.11, float(r.loc[">400 combined fraction","value"]), 0.006)

SUBT = {
    "dM0": (13.0,  5.0, 2.57), "dM1": (14.4,  5.6, 2.57), "dM2": (20.8, 9.2, 2.25),
    "dM3": (23.6, 11.6, 2.03), "dM4": (25.1, 14.5, 1.73), "dM5": (64.6, 36.4, 1.77),
    "dM6": (74.8, 52.6, 1.42),
}
for k, (near, far, ratio) in SUBT.items():
    ck("%s near" % k, near, float(s.loc[k,"near_cut_pct"]), 0.05)
    ck("%s far" % k, far, float(s.loc[k,"far_cut_pct"]), 0.05)
    ck("%s ratio" % k, ratio, float(s.loc[k,"ratio_cut"]), 0.006)

ck("factor range lower", 1.42, round(float(s.ratio_cut.min()),2), 0.006)
ck("factor range upper", 2.57, round(float(s.ratio_cut.max()),2), 0.006)
ck("dM0 n inner", 35091, int(s.loc["dM0","n_inner"]), 0)
ck("dM0 n far", 65497, int(s.loc["dM0","n_far"]), 0)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step4_check.csv"), index=False, encoding="utf-8-sig")
bad = res[res.status == "FAIL"]
print(bad.to_string(index=False) if len(bad) else "  (no mismatch)")
print("  PASS %d / %d" % (int((res.status=="PASS").sum()), len(res)))
raise SystemExit(0 if len(bad) == 0 else 2)
