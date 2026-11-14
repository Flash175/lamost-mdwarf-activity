
# -*- coding: utf-8 -*-
"""step2_check_43.py — §4.3 的逐项比对（正文 §4.3 与表 4）。"""
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
def ck_le(item, bound, got):
    ok = abs(float(got)) <= float(bound) + 1e-9
    rows.append(dict(item=item, body=bound, computed=round(float(got), 4), diff=np.nan,
                     status="PASS" if ok else "STALE"))

LAB = ("100-200", "200-400")

# 样本与降幅
ck_int("sample (470,548)", 470548, V["n sample"])
ck("dM0-dM6 share (99.9 %)", 99.9, float(V["dM0-dM6 share"]), tol=0.05)
ck_le("Shapley sum equals full-model difference (1e-9)", 1e-9, float(V["shapley identity err"]))
ck("observed decline 100-200 (10.92)", 10.92, float(V["obs decline 100-200"]), tol=0.006)
ck("observed decline 200-400 (16.09)", 16.09, float(V["obs decline 200-400"]), tol=0.006)
ck("full-sample decline 100-200 (10.90)", 10.90, float(V["full decline 100-200"]), tol=0.006)
ck("full-sample decline 200-400 (16.03)", 16.03, float(V["full decline 200-400"]), tol=0.006)

# 表 4：六个量
vals = {
    "F1": (1.78, 2.63, 16.3, 16.3),
    "F2": (0.86, 1.99, 7.9, 12.4),
    "F3": (0.13, 0.32, 1.2, 2.0),
    "F4": (0.05, 0.16, 0.4, 1.0),
}
for k, (a, b, sa, sb) in vals.items():
    ck("%s 100-200 (%s)" % (k, a), a, float(V["contrib 100-200 %s" % k]), tol=0.006)
    ck("%s 200-400 (%s)" % (k, b), b, float(V["contrib 200-400 %s" % k]), tol=0.006)
    ck("%s share 100-200 (%s)" % (k, sa), sa, 100.0*float(V["contrib 100-200 %s" % k])/float(V["obs decline 100-200"]), tol=0.06)
    ck("%s share 200-400 (%s)" % (k, sb), sb, 100.0*float(V["contrib 200-400 %s" % k])/float(V["obs decline 200-400"]), tol=0.06)
ck_le("log g and [Fe/H] at or below 0.32 pp", 0.32,
      max(max(float(V["contrib %s F3" % l]), float(V["contrib %s F4" % l])) for l in LAB))

# F5 / F6 单列与合并的测量渠道
ck("F5 100-200 (-0.35)", -0.35, float(V["contrib 100-200 F5"]), tol=0.006)
ck("F5 200-400 (-0.84)", -0.84, float(V["contrib 200-400 F5"]), tol=0.006)
ck("F6 100-200 (+0.93)", 0.93, float(V["contrib 100-200 F6"]), tol=0.006)
ck("F6 200-400 (+1.90)", 1.90, float(V["contrib 200-400 F6"]), tol=0.006)
ck("measurement channel 100-200 (0.58)", 0.58, float(V["meas 100-200"]), tol=0.006)
ck("measurement channel 200-400 (1.06)", 1.06, float(V["meas 200-400"]), tol=0.006)

# 合计、份额、残余
ck("explained 100-200 (3.41)", 3.41, float(V["explained 100-200"]), tol=0.006)
ck("explained 200-400 (6.16)", 6.16, float(V["explained 200-400"]), tol=0.006)
ck("explained share 100-200 (31.2)", 31.2, float(V["explained share 100-200"]), tol=0.06)
ck("explained share 200-400 (38.3)", 38.3, float(V["explained share 200-400"]), tol=0.06)
ck("unexplained 100-200 (7.52)", 7.52, float(V["unexplained 100-200"]), tol=0.006)
ck("unexplained 200-400 (9.93)", 9.93, float(V["unexplained 200-400"]), tol=0.006)
ck("unexplained share 100-200 (68.8)", 68.8, float(V["unexplained share 100-200"]), tol=0.06)
ck("unexplained share 200-400 (61.7)", 61.7, float(V["unexplained share 200-400"]), tol=0.06)

# 与归档主表的一致性与 bootstrap 半宽
ck_le("matches archived table (<=0.001 pp)", 0.001, float(V["max |computed - archive| pp"]))
ck_le("bootstrap individual entry <= 0.06 pp", 0.06, float(V["boot half F max"]))
ck("bootstrap explained total (0.08)", 0.08, float(V["boot half explained max"]), tol=0.006)
ck("bootstrap unexplained (0.08)", 0.08, float(V["boot half unexplained max"]), tol=0.006)

# 稳健性
ck("R3 logit 100-200 (3.27)", 3.27, float(V["R3 memo explained 100-200"]), tol=0.006)
ck("R3 logit 200-400 (6.11)", 6.11, float(V["R3 memo explained 200-400"]), tol=0.006)
ck_le("R3 own recomputation within 0.2 pp of memo (100-200)", 0.2, float(V["R3 own vs memo diff 100-200"]))
ck_le("R3 own recomputation within 0.2 pp of memo (200-400)", 0.2, float(V["R3 own vs memo diff 200-400"]))
ck("R5 edges 125/250 explained 100-200 (4.08)", 4.08, float(V["R5 explained 100-200"]), tol=0.006)
ck("R5 edges 125/250 explained 200-400 (6.17)", 6.17, float(V["R5 explained 200-400"]), tol=0.006)
ck("R2 no distance term 100-200 (3.48)", 3.48, float(V["R2 explained 100-200"]), tol=0.006)
ck("R2 no distance term 200-400 (6.30)", 6.30, float(V["R2 explained 200-400"]), tol=0.006)
ck("explained-total change lower (0.05)", 0.05, float(V["variant change logit min pp"]), tol=0.006)
ck("explained-total change upper (0.68)", 0.68, float(V["variant change max pp"]), tol=0.006)
ck("explained-total change lower incl. R5 (0.02)", 0.02, float(V["variant change all min pp"]), tol=0.005)
ck_le("spectral subtype term moves at most 0.4 pp", 0.4,
      max(float(V["R5 F1 change max"]), float(V["logit max factor change (memo)"])))

# 高 S/N 子样本（R8）
ck_int("high-S/N subsample (42,493)", 42493, V["R8 n"])
ck("R8 F5 100-200 (+0.06)", 0.06, float(V["R8 F5 100-200"]), tol=0.006)
ck("R8 F5 200-400 (+0.23)", 0.23, float(V["R8 F5 200-400"]), tol=0.006)
ck("R8 F6 100-200 (+0.03)", 0.03, float(V["R8 F6 100-200"]), tol=0.006)
ck("R8 F6 200-400 (+0.02)", 0.02, float(V["R8 F6 200-400"]), tol=0.006)
ck("R8 observed 100-200 (6.00)", 6.00, float(V["R8 obs 100-200"]), tol=0.006)
ck("R8 observed 200-400 (8.51)", 8.51, float(V["R8 obs 200-400"]), tol=0.006)

# 与两渠道分解对齐（正文 §5 引用的 1.62 / 3.53）
ck("within-subtype absorbed 100-200 (1.62)", 1.62, float(V["within explained 100-200"]), tol=0.006)
ck("within-subtype absorbed 200-400 (3.53)", 3.53, float(V["within explained 200-400"]), tol=0.006)

# 表 4 逐格核对（印刷值 vs 复现值，按印刷精度判定：pp 两位小数、份额一位小数）
obs = {l: float(V["obs decline %s" % l]) for l in LAB}
printed = [
    ("100-200", "Spectral subtype (F1)", "contrib 100-200 F1", 1.78, 16.3),
    ("100-200", "Teff (F2)", "contrib 100-200 F2", 0.86, 7.9),
    ("100-200", "log g (F3)", "contrib 100-200 F3", 0.13, 1.2),
    ("100-200", "[Fe/H] (F4)", "contrib 100-200 F4", 0.05, 0.4),
    ("100-200", "Measurement channel (F5+F6)", "meas 100-200", 0.58, 5.3),
    ("100-200", "Explained (sum)", "explained 100-200", 3.41, 31.2),
    ("100-200", "Observed decline", "obs decline 100-200", 10.92, 100.0),
    ("100-200", "Unexplained", "unexplained 100-200", 7.52, 68.8),
    ("200-400", "Spectral subtype (F1)", "contrib 200-400 F1", 2.63, 16.3),
    ("200-400", "Teff (F2)", "contrib 200-400 F2", 1.99, 12.4),
    ("200-400", "log g (F3)", "contrib 200-400 F3", 0.32, 2.0),
    ("200-400", "[Fe/H] (F4)", "contrib 200-400 F4", 0.16, 1.0),
    ("200-400", "Measurement channel (F5+F6)", "meas 200-400", 1.06, 6.6),
    ("200-400", "Explained (sum)", "explained 200-400", 6.16, 38.3),
    ("200-400", "Observed decline", "obs decline 200-400", 16.09, 100.0),
    ("200-400", "Unexplained", "unexplained 200-400", 9.93, 61.7),
]
t4rows = []
for lab, name, key, pp, sh in printed:
    c = float(V[key]); cs = 100.0*c/obs[lab]
    ok = (abs(c - pp) <= 0.005) and (abs(cs - sh) <= 0.05)
    t4rows.append(dict(bin=lab, quantity=name, printed_pp=pp, computed_pp=round(c, 4),
                       printed_share_pct=sh, computed_share_pct=round(cs, 2),
                       status="PASS" if ok else "CHECK"))
T4 = pd.DataFrame(t4rows)
T4.to_csv(os.path.join(OUT, "表4逐格核对.csv"), index=False, encoding="utf-8-sig")
print("表 4 逐格：%d 格，%s" % (len(T4), dict(T4.status.value_counts())))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
