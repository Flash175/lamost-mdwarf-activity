# -*- coding: utf-8 -*-
"""step3_check_25.py — 把 §2.5 / Table 3 与 §4.2 引用的数字逐项比对。"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
t3 = pd.read_csv(os.path.join(OUT, "step1_table3.csv")).set_index("subtype")
q = pd.read_csv(os.path.join(OUT, "step2_quantile.csv")).set_index("subtype")

EXPECT = {   # 正文 Table 3
    "dM0": (201622, -0.12, 0.55,  8.0, 92),
    "dM1": (157205, -0.01, 0.72,  9.6, 90),
    "dM2": (109103,  0.11, 1.44, 15.4, 85),
    "dM3": ( 35568,  0.18, 2.10, 19.8, 80),
    "dM4": (107251,  0.24, 3.00, 23.4, 77),
    "dM5": ( 22932,  1.77, 7.28, 61.2, 39),
    "dM6": (  3671,  4.38, 10.41, 77.0, 23),
}
rows = []
def ck(item, exp, got, tol=0.005):
    ok = abs(float(got) - float(exp)) <= tol
    rows.append(dict(item=item, expected=exp, computed=got, status="PASS" if ok else "FAIL"))

for s, (n, med, p90, frac, rank) in EXPECT.items():
    ck("%s n" % s, n, int(t3.loc[s, "n"]), 0)
    ck("%s median" % s, med, round(float(t3.loc[s, "ew_med"]), 2), 0.006)
    ck("%s p90" % s, p90, round(float(t3.loc[s, "ew_p90"]), 2), 0.006)
    ck("%s frac>0.75" % s, frac, round(float(t3.loc[s, "frac_gt_075"]), 1), 0.05)
    ck("%s percentile" % s, rank, int(t3.loc[s, "pct_rank"]), 0)

ck("sum n (dM0-dM6)", 637352, int(t3["n"].sum()), 0)
ck("dM7-dM9 residual", 621, 637973 - int(t3["n"].sum()), 0)

# §4.2 引用的分位判据结论
sub = ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]
r = q.loc[sub, "ratio_000_to_200_400"]
ck("下降倍数下限 (≥1.4)", 1.4, round(float(r.min()), 2), 0.05)
ck("下降倍数上限 (≤3.7)", 3.7, round(float(r.max()), 2), 0.05)
ck("分位判据下亚型跨度 (pp)", 7.7, round(float(q.loc[sub, "q80_<100"].max() - q.loc[sub, "q80_<100"].min()), 1), 0.1)
ck("固定阈值下亚型跨度 (pp)", 69.0, round(float(t3["frac_gt_075"].max() - t3["frac_gt_075"].min()), 1), 0.1)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step3_check.csv"), index=False, encoding="utf-8-sig")
bad = res[res.status == "FAIL"]
print(bad.to_string(index=False) if len(bad) else "  (no mismatch)")
print("  PASS %d / %d" % (int((res.status == "PASS").sum()), len(res)))
raise SystemExit(0 if len(bad) == 0 else 2)
