# -*- coding: utf-8 -*-
"""step2_check_42.py — §4.2 的逐项比对（38 项）。

输入：输出/step1_values.csv
输出：输出/step2_check.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
V = pd.read_csv(os.path.join(OUT, "step1_values.csv"), encoding="utf-8-sig").set_index("item")["value"]

rows = []
def ck(item, body, got, tol=0.06, hard=0.5):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 3),
                     diff=round(float(got)-float(body), 3), status=st))
def ck_int(item, body, got):
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got)-int(body),
                     status="PASS" if int(got) == int(body) else "STALE"))

# 1) dM0 占比
for lab, v in (("in100", 19.8), ("m200", 34.3), ("m400", 54.0), ("m600", 71.8)):
    ck("dM0 share %s" % lab, v, float(V["dM0 share %s" % lab]), tol=0.05)

# 2) 逐型倍数
for sp, v in (("dM0", 2.57), ("dM1", 2.57), ("dM2", 2.25), ("dM3", 2.03),
              ("dM4", 1.73), ("dM5", 1.77), ("dM6", 1.42)):
    ck("factor %s" % sp, v, float(V["factor %s" % sp]), tol=0.006)
ck("Wilson intervals disjoint (all 7)", 1.0, float(V["wilson disjoint all"]), tol=0.0)

# 3) 单历元
ck_int("single-epoch stars (502,986)", 502986, V["single epoch n"])
ck("single-epoch pct (~79)", 78.84, float(V["single epoch pct"]), tol=0.2)
ck("single-epoch factor min (1.46)", 1.46, float(V["single-epoch factor min"]), tol=0.006)
ck("single-epoch factor max (2.82)", 2.82, float(V["single-epoch factor max"]), tol=0.006)
ck_int("single-epoch dM6 far (99)", 99, V["single-epoch nfar dM6"])

# 4) |Z| 与 d sin b 的相关（银心 z）
ck("corr(|z_gc|, d sin b) 0.994", 0.994, float(V["corr zgc dsinb"]), tol=0.0006)

# 5) 分类检验 C1 / C2
for sp, v in (("dM0", 2.43), ("dM1", 2.75), ("dM2", 2.27), ("dM3", 1.83),
              ("dM4", 1.58), ("dM5", 1.13)):
    ck("C1 ratio %s" % sp, v, float(V["C1 ratio %s" % sp]), tol=0.006)
ck_int("C1 dM5 far stars (672)", 672, V["C1 dM5 nfar"])
for sp, v in (("dM0", 2.55), ("dM1", 2.32), ("dM4", 1.39), ("dM5", 1.08)):
    ck("C2 ratio %s" % sp, v, float(V["C2 ratio %s" % sp]), tol=0.006)

# 6) L 判据
for sp, v in (("dM5", 36.4), ("dM6", 36.3), ("dM7", 3.1), ("dM8", 4.7),
              ("dM0", 6.5), ("dM4", 14.3)):
    ck("L pct %s" % sp, v, float(V["L pct %s" % sp]), tol=0.05)
ck("spread EW dM0-dM6 (65.8)", 65.8, float(V["spread EW dM0-dM6"]), tol=0.05)
ck("spread L dM0-dM6 (29.9)", 29.9, float(V["spread L dM0-dM6"]), tol=0.05)
ck("spread EW incl dM7 (67.7)", 67.7, float(V["spread EW incl dM7"]), tol=0.05)
ck("spread L incl dM7 (33.2)", 33.2, float(V["spread L incl dM7"]), tol=0.05)
ck("L factor in100/m200 (3.21)", 3.21, float(V["L factor in100 over m200"]), tol=0.01)

# 现象数（分箱与倍数）
for zb, v in (("<100", 23.81), ("100-200", 12.90), ("200-400", 7.78), ("400-600", 6.65), ("600-1000", 7.76)):
    ck("bin f %s" % zb, v, float(V["bin f %s" % zb]), tol=0.006)
ck("factor 3.06", 3.06, float(V["factor 3.06"]), tol=0.006)
ck("factor 3.58", 3.58, float(V["factor 3.58"]), tol=0.006)

# 0.75 A 的分位
ck("0.75 A percentile dM0 (92)", 92.0, float(V["pct075 dM0"]), tol=0.1)
ck("0.75 A percentile dM6 (23)", 23.0, float(V["pct075 dM6"]), tol=0.1)

# 分位判据倍数范围
ck("quantile ratio lower (1.4)", 1.4, float(V["quantile ratio min"]), tol=0.05)
ck("quantile ratio upper (3.7)", 3.7, float(V["quantile ratio max"]), tol=0.05)

# 16 点分解
ck("16.0 points drop", 16.0, float(V["drop 200-400"]), tol=0.1)
ck("10.7 within subtypes", 10.7, float(V["within 200-400"]), tol=0.05)
ck("5.4 from composition", 5.4, float(V["comp 200-400"]), tol=0.05)
ck("components sum 16.09", 16.09, float(V["sum components"]), tol=0.01)
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step2_check.csv")
