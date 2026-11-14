
# -*- coding: utf-8 -*-
"""step2_check_45.py — §4.5 的逐项比对。"""
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

# 三个内盘切片：std / MAD 的差、误差与显著性
for lab, d, se, z, dm, sem, zm in (("0-25", 6.08, 0.23, -26.0, 7.82, 0.62, -12.7),
                                   ("25-50", 4.03, 0.28, -14.2, 6.43, 0.74, -8.6),
                                   ("50-75", 4.39, 0.29, -15.3, 6.49, 0.75, -8.6)):
    ck("std diff %s (%s)" % (lab, d), d, float(V["std diff %s" % lab]), tol=0.006)
    ck("std se %s (%s)" % (lab, se), se, float(V["std se %s" % lab]), tol=0.006)
    ck("std z %s (%s)" % (lab, z), z, float(V["std z %s" % lab]), tol=0.06)
    ck("mad diff %s (%s)" % (lab, dm), dm, float(V["mad diff %s" % lab]), tol=0.006)
    ck("mad se %s (%s)" % (lab, sem), sem, float(V["mad se %s" % lab]), tol=0.006)
    ck("mad z %s (%s)" % (lab, zm), zm, float(V["mad z %s" % lab]), tol=0.06)

# 档内 |Z| 对齐
for lab, a, b, dd in (("0-25", 12.1, 12.4, 0.3), ("25-50", 36.6, 37.7, 1.1), ("50-75", 61.7, 62.5, 0.9)):
    ck("zmed active %s (%s)" % (lab, a), a, float(V["zmed act %s" % lab]), tol=0.06)
    ck("zmed inactive %s (%s)" % (lab, b), b, float(V["zmed inact %s" % lab]), tol=0.06)
    ck("zmed difference %s (%s)" % (lab, dd), dd, float(V["zmed diff %s" % lab]), tol=0.06)

# 亚型分层
ck("strata dM0-dM2 lower z (7.7)", 7.7, float(V["strata dM0-dM2 zmin"]), tol=0.06)
ck("strata dM0-dM2 upper z (14.2)", 14.2, float(V["strata dM0-dM2 zmax"]), tol=0.06)
ck("strata dM3-dM5 lower z (9.1)", 9.1, float(V["strata dM3-dM5 zmin"]), tol=0.06)
ck("strata dM3-dM5 upper z (19.7)", 19.7, float(V["strata dM3-dM5 zmax"]), tol=0.06)

# 75 pc 以外两种估计量的分歧
ck("std z 75-100 (+1.9)", 1.9, float(V["std z far 75-100"]), tol=0.06)
ck("std z 100-150 (+0.7)", 0.7, float(V["std z far 100-150"]), tol=0.06)
ck("std z 150-200 (+5.9)", 5.9, float(V["std z far 150-200"]), tol=0.06)
ck("mad z 75-100 (-4.7)", -4.7, float(V["mad z far 75-100"]), tol=0.06)
ck("mad z 100-150 (-4.4)", -4.4, float(V["mad z far 100-150"]), tol=0.06)

# 内盘非活动 σ(vz) 随亚型
ck("inactive svz lower (19.66)", 19.66, float(V["inactive svz min"]), tol=0.006)
ck("inactive svz upper (21.14)", 21.14, float(V["inactive svz max"]), tol=0.006)
ck("inactive svz range (1.5)", 1.5, float(V["inactive svz range"]), tol=0.05)
ck("active fraction dM0 (9.1)", 9.1, float(V["active frac dM0"]), tol=0.06)
ck("active fraction dM5 (70.4)", 70.4, float(V["active frac dM5"]), tol=0.06)
ck("active fraction factor (7.8)", 7.8, float(V["active frac factor"]), tol=0.06)

# 6D 子样本
ck_int("6D subsample (82,518)", 82518, V["6D n"])
ck("6D inside 400 pc (99.9 %)", 99.9, float(V["6D frac inside 400"]), tol=0.05)
ck_int("6D beyond 400 pc (47)", 47, V["6D n beyond 400"])

# 亮端偏置：Gaia RV 覆盖率（冻结 2D 样本重算；正文原写 95 / 55 / 0.02）
ck("RV coverage G 12-14 (97)", 97.0, float(V["rv cover 12-14"]), tol=0.5)
ck("RV coverage G 14-16 (53)", 53.0, float(V["rv cover 14-16"]), tol=0.5)
ck("RV coverage G 16-18 (0.003)", 0.003, float(V["rv cover 16-18"]), tol=0.001)

# 两类各自的绝对弥散与样本量（zmatch 表，独立于 diff 表）
ck_int("n active 0-25 (2,589)", 2589, V["n act 0-25"])
ck_int("n inactive 0-25 (9,425)", 9425, V["n inact 0-25"])
ck_int("n active 25-50 (2,020)", 2020, V["n act 25-50"])
ck_int("n inactive 50-75 (10,326)", 10326, V["n inact 50-75"])
for sl, a, b in (("0-25", 13.40, 19.49), ("25-50", 15.65, 19.68), ("50-75", 15.39, 19.78)):
    ck("sigma(vz) active %s (%s)" % (sl, a), a, float(V["svz act %s" % sl]), tol=0.006)
    ck("sigma(vz) inactive %s (%s)" % (sl, b), b, float(V["svz inact %s" % sl]), tol=0.006)
for sl, a, b in (("0-25", 8.65, 16.48), ("25-50", 9.92, 16.36), ("50-75", 10.21, 16.70)):
    ck("MAD active %s (%s)" % (sl, a), a, float(V["mad act %s" % sl]), tol=0.006)
    ck("MAD inactive %s (%s)" % (sl, b), b, float(V["mad inact %s" % sl]), tol=0.006)

# 两张归档表互校：Δσ 由弥散列重算 vs diff 表
for sl, v in (("0-25", -6.08), ("25-50", -4.03), ("50-75", -4.39)):
    ck("delta std %s from dispersions (%s)" % (sl, v), v, float(V["d recomputed %s" % sl]), tol=0.006)
for sl, v in (("0-25", -7.82), ("25-50", -6.43), ("50-75", -6.50)):
    ck("delta MAD %s from dispersions (%s)" % (sl, v), v, float(V["d mad recomputed %s" % sl]), tol=0.006)

# 误差模型重算：SE(std) 与 SE(MAD)=sqrt((1.858 sigma)^2/n) 是否复现归档值
for sl in ("0-25", "25-50", "50-75"):
    ck("SE(std) recomputed %s" % sl, float(V["se archived %s" % sl]), float(V["se recomputed %s" % sl]), tol=0.002)
    ck("SE(MAD) recomputed %s" % sl, float(V["se mad archived %s" % sl]), float(V["se mad recomputed %s" % sl]), tol=0.002)

# 40-80 pc 共同窗口（正文"0-75 pc 内两估计量一致"的独立支撑）
ck("clean window 40-80 std (-3.91)", -3.91, float(V["d recomputed 40-80"]), tol=0.006)
ck("clean window 40-80 z (-16.8)", -16.8, float(V["w40-80 z std"]), tol=0.06)
ck("clean window 40-80 MAD (-6.02)", -6.02, float(V["d mad recomputed 40-80"]), tol=0.006)

# 独立重算（从原始数据重建 6D 子样本）
ck_int("recomputed 6D subsample (82,518)", 82518, V["re 6D n"])
ck_int("recomputed n active 0-25 (2,589)", 2589, V["re n act 0-25"])
for lab, a, b in (("0-25", 13.40, 19.49), ("25-50", 15.65, 19.68), ("50-75", 15.39, 19.78)):
    ck("recomputed sigma(vz) active %s (%s)" % (lab, a), a, float(V["re svz act %s" % lab]), tol=0.006)
    ck("recomputed sigma(vz) inactive %s (%s)" % (lab, b), b, float(V["re svz inact %s" % lab]), tol=0.006)
for lab, a, b in (("0-25", 8.65, 16.48), ("25-50", 9.92, 16.36), ("50-75", 10.21, 16.70)):
    ck("recomputed MAD active %s (%s)" % (lab, a), a, float(V["re mad act %s" % lab]), tol=0.006)
    ck("recomputed MAD inactive %s (%s)" % (lab, b), b, float(V["re mad inact %s" % lab]), tol=0.006)
for lab, v in (("0-25", -6.08), ("25-50", -4.03), ("50-75", -4.39)):
    ck("recomputed delta std %s (%s)" % (lab, v), v, float(V["re d std %s" % lab]), tol=0.006)
for lab, v in (("0-25", -7.82), ("25-50", -6.43), ("50-75", -6.49)):
    ck("recomputed delta MAD %s (%s)" % (lab, v), v, float(V["re d mad %s" % lab]), tol=0.006)
for lab, a, b in (("0-25", 12.1, 12.4), ("25-50", 36.6, 37.7), ("50-75", 61.7, 62.5)):
    ck("recomputed median |z| active %s (%s)" % (lab, a), a, float(V["re zmed act %s" % lab]), tol=0.06)
    ck("recomputed median |z| inactive %s (%s)" % (lab, b), b, float(V["re zmed inact %s" % lab]), tol=0.06)
for lab, v in (("0-25", -26.0), ("25-50", -14.2), ("50-75", -15.3)):
    ck("recomputed z (std) %s (%s)" % (lab, v), v, float(V["re z std %s" % lab]), tol=0.06)
for lab, v in (("0-25", -12.7), ("25-50", -8.6), ("50-75", -8.6)):
    ck("recomputed z (MAD) %s (%s)" % (lab, v), v, float(V["re z mad %s" % lab]), tol=0.06)
ck("recomputed 40-80 delta std (-3.91)", -3.91, float(V["re d std 40-80"]), tol=0.006)
ck("recomputed 40-80 delta MAD (-6.02)", -6.02, float(V["re d mad 40-80"]), tol=0.006)
ck("recomputed 40-80 z std (-16.8)", -16.8, float(V["re z std 40-80"]), tol=0.06)
ck("recomputed far z std 75-100 (+1.9)", 1.9, float(V["re z std 75-100"]), tol=0.06)
ck("recomputed far z std 100-150 (+0.7)", 0.7, float(V["re z std 100-150"]), tol=0.06)
ck("recomputed far z std 150-200 (+5.9)", 5.9, float(V["re z std 150-200"]), tol=0.06)
ck("recomputed far z MAD 75-100 (-4.7)", -4.7, float(V["re z mad 75-100"]), tol=0.06)
ck("recomputed far z MAD 100-150 (-4.4)", -4.4, float(V["re z mad 100-150"]), tol=0.06)


# ---- §4.12 距离精度（由 step1b 从原始星表重算） ----
try:
    DB = pd.read_csv(os.path.join(OUT, "step1b_distance.csv"), encoding="utf-8-sig").set_index("item")["value"]
    ck("median relative distance error (2.3 per cent)", 2.3, float(DB["median relative distance error (per cent)"]), tol=0.05)
    ck("worst case at the plx/err>5 cut (20 per cent)", 20.0, float(DB["worst case at the plx/err>5 cut (per cent)"]), tol=0.01)
    ck("95th percentile (8.6 per cent)", 8.6, float(DB["95th percentile (per cent)"]), tol=0.1)
    ck("fraction above 10 per cent (2.9 per cent)", 2.9, 100.0 * float(DB["fraction above 10 per cent"]), tol=0.1)
except Exception as e:
    print("  (step1b_distance.csv missing: %s)" % type(e).__name__)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
