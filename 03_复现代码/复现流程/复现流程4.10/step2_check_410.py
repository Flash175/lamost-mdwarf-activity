
# -*- coding: utf-8 -*-
"""step2_check_410.py — §4.10 的逐项比对。"""
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

# 覆盖
ck_int("usable [Fe/H] sample (584,029)", 584029, V["n usable"])
ck("coverage (91.5 per cent)", 91.5, float(V["cov overall"]), tol=0.05)
ck("median [Fe/H] uncertainty (0.12 dex)", 0.12, float(V["median mh err"]), tol=0.006)

# 可用子集 vs 全样本的活动率
for zb, mh_v, all_v in (("<100", 23.29, 23.81), ("100-200", 12.57, 12.90), ("200-400", 7.51, 7.78),
                        ("400-600", 6.42, 6.65), ("600-1000", 7.39, 7.76), (">1000", 12.34, 13.24)):
    ck("active fraction, [Fe/H] subset %s (%s)" % (zb, mh_v), mh_v, float(V["act mh %s" % zb]), tol=0.02)
    ck("active fraction, full sample %s (%s)" % (zb, all_v), all_v, float(V["act all %s" % zb]), tol=0.02)
ck("decline in the [Fe/H] subset (15.8 pp)", 15.8, float(V["decline mh subset"]), tol=0.06)

# 垂直梯度
for zb, v in (("<100", -0.085), ("100-200", -0.114), ("200-400", -0.154), ("400-600", -0.201), ("600-1000", -0.203)):
    ck("median [Fe/H] %s (%s dex)" % (zb, v), v, float(V["mh med %s" % zb]), tol=0.006)
ck("gradient over the first kpc (0.12 dex)", 0.12, float(V["gradient first kpc"]), tol=0.006)
ck("linear slope (-2.3e-4 dex/pc)", -2.3e-4, float(V["slope dex per pc"]), tol=0.05e-4)
ck("Pearson r (-0.195)", -0.195, float(V["pearson r"]), tol=0.002)

# 固定 |Z| 的四分位
for zb, v in (("<100", 5.7), ("100-200", 3.9), ("200-400", 2.7), ("400-600", 2.2), ("600-1000", 0.9)):
    ck("top minus bottom [Fe/H] quartile %s (%s pp)" % (zb, v), v, float(V["q4-q1 %s" % zb]), tol=0.06)

# 通道预算
ck("global logistic slope (+0.588)", 0.588, float(V["beta global"]), tol=0.002)
ck("global slope error (0.020)", 0.020, float(V["beta global se"]), tol=0.002)
ck("global contribution (-0.5 pp)", -0.5, float(V["delta pp global"]), tol=0.02)
ck("Delta[Fe/H] innermost to 200-400 pc (-0.069 dex)", -0.069, float(V["dmh"]), tol=0.002)
ck("fixed-effect slope (-1.51)", -1.51, float(V["beta FE"]), tol=0.006)
ck("fixed-effect slope error (0.024)", 0.024, float(V["beta FE se"]), tol=0.002)
ck("fixed-effect contribution (+1.3 pp)", 1.3, float(V["delta pp FE"]), tol=0.02)

# 亚型退化
ck("median [Fe/H] at dM0 (-0.21 dex)", -0.21, float(V["mh med dM0"]), tol=0.006)
ck("median [Fe/H] at dM6 (+0.23 dex)", 0.23, float(V["mh med dM6"]), tol=0.006)
ck("correlation with Teff at dM0 (+0.42)", 0.42, float(V["rho mh teff dM0"]), tol=0.006)
ck("correlation with Teff at dM5 (-0.43)", -0.43, float(V["rho mh teff dM5"]), tol=0.006)
for sp, v in (("dM0", 0.55), ("dM1", -1.52), ("dM2", -0.83), ("dM3", -1.34),
              ("dM4", -2.07), ("dM5", -2.21), ("dM6", -4.14)):
    ck("per-subtype slope %s (%s dex^-1)" % (sp, v), v, float(V["slope %s" % sp]), tol=0.006)
ck("per-subtype slope error at dM0 (0.07)", 0.07, float(V["slope dM0 se"]), tol=0.006)

# 四分位组成（<100 pc）
ck("late-type share, metal-poor quartile (15 per cent, text rounds to integers)", 15.0, float(V["late share q1"]), tol=0.5)
ck("late-type share, metal-rich quartile (69 per cent, text rounds to integers)", 69.0, float(V["late share q4"]), tol=0.5)
ck("dM4 share, metal-poor quartile (13 per cent, text rounds to integers)", 13.0, float(V["dM4 share q1"]), tol=0.5)
ck("dM4 share, metal-rich quartile (46 per cent, text rounds to integers)", 46.0, float(dM4q4) if False else float(V["dM4 share q4"]), tol=0.5)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
