
# -*- coding: utf-8 -*-
"""step2_check_46.py — §4.6 的逐项比对。"""
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
def ck_bool(item, body, got):
    ok = (bool(body) == bool(round(float(got))))
    rows.append(dict(item=item, body=body, computed=bool(round(float(got))), diff=np.nan,
                     status="PASS" if ok else "STALE"))

# 冷端三档的活动率与区间（正文 72.9 (71.4-74.3)、74.7 (70.7-78.4)、65.9 (55.3-75.1)）
for sp, f, lo, hi in (("dM6", 72.9, 71.4, 74.3), ("dM7", 74.7, 70.7, 78.4), ("dM8", 65.9, 55.3, 75.1)):
    ck("%s active fraction (%s)" % (sp, f), f, float(V["frac %s" % sp]), tol=0.06)
    ck("%s Wilson lo (%s)" % (sp, lo), lo, float(V["lo %s" % sp]), tol=0.06)
    ck("%s Wilson hi (%s)" % (sp, hi), hi, float(V["hi %s" % sp]), tol=0.06)

# 区间重叠："no turnover between dM6 and dM8 can be established"
ck_bool("Wilson intervals dM6 vs dM7 overlap", True, float(V["overlap 6-7"]))
ck_bool("Wilson intervals dM7 vs dM8 overlap", True, float(V["overlap 7-8"]))

# 可探测比例：63.9 -> 30.6
ck("detectable fraction dM6 (63.9)", 63.9, float(V["detect dM6"]), tol=0.06)
ck("detectable fraction dM8 (30.6)", 30.6, float(V["detect dM8"]), tol=0.06)

# S/N 匹配后的 dM6/dM7 对
for cut, a, b in (("3-6", 73.6, 76.0), ("6-12", 71.9, 68.9), ("6-20", 73.4, 75.3)):
    ck("matched S/N %s dM6 (%s)" % (cut, a), a, float(V["matched %s dM6" % cut]), tol=0.06)
    ck("matched S/N %s dM7 (%s)" % (cut, b), b, float(V["matched %s dM7" % cut]), tol=0.06)
ck_int("dM8 survivors S/N 3-6 (15)", 15, V["matched n 3-6 dM8"])
ck_int("dM8 survivors S/N 6-12 (6)", 6, V["matched n 6-12 dM8"])
ck_int("dM8 survivors S/N 6-20 (7)", 7, V["matched n 6-20 dM8"])
ck("dM8 survivors lower bound (6)", 6.0, float(V["matched n dM8 min"]))
ck("dM8 survivors upper bound (15)", 15.0, float(V["matched n dM8 max"]))

# 第三判据：2,498 条存档谱，活动谱里 19.4% 不通过（437/2,255）
ck_int("criterion-3 archived spectra (2,498)", 2498, V["c3 n sample"])
ck_int("criterion-3 active spectra (2,255)", 2255, V["c3 n active"])
ck_int("criterion-3 active failures (437)", 437, V["c3 n fail"])
ck("criterion-3 failure rate (19.4 %)", 19.4, float(V["c3 fail pct"]), tol=0.06)
ck_int("criterion-3 inactive spectra (243)", 243, V["c3 n inactive"])
ck("inactive spectra passing (18.9 %)", 18.9, float(V["c3 inactive pass pct"]), tol=0.06)

# 第三判据：5,999 条抽样谱，|Z| 分档 13.1 / 22.2 / 36.4 / 36.8
ck_int("criterion-3 draw (5,999)", 5999, V["c3 zb n"])
for zb, v in (("<100", 13.1), ("100-200", 22.2), ("200-400", 36.4), (">400", 36.8)):
    ck("criterion-3 failure %s pc (%s)" % (zb, v), v, float(V["c3 zb fail %s" % zb]), tol=0.06)

# 第三判据：逐亚型 19.7 / 26.4 / 17.8 / 26.0 / 19.9 / 11.8 / 13.2
for sp, v in (("dM0", 19.7), ("dM1", 26.4), ("dM2", 17.8), ("dM3", 26.0),
              ("dM4", 19.9), ("dM5", 11.8), ("dM6", 13.2)):
    ck("criterion-3 failure %s (%s)" % (sp, v), v, float(V["c3 sub fail %s" % sp]), tol=0.06)

# 独立重算（不读归档表，直接从 637,973 冻结样本重建）
ck_int("recomputed 2D sample (637,973)", 637973, V["re n 2D"])
ck_int("recomputed n dM5 (22,932)", 22932, V["re n dM5"])
ck_int("recomputed n dM6 (3,671)", 3671, V["re n dM6"])
ck_int("recomputed n dM7 (483)", 483, V["re n dM7"])
ck_int("recomputed n dM8 (85)", 85, V["re n dM8"])
for sp, f in (("dM6", 72.9), ("dM7", 74.7), ("dM8", 65.9)):
    ck("recomputed fraction %s (%s)" % (sp, f), f, float(V["re frac %s" % sp]), tol=0.06)
for sp, lo, hi in (("dM6", 71.4, 74.3), ("dM7", 70.7, 78.4), ("dM8", 55.3, 75.1)):
    ck("recomputed Wilson lo %s (%s)" % (sp, lo), lo, float(V["re lo %s" % sp]), tol=0.06)
    ck("recomputed Wilson hi %s (%s)" % (sp, hi), hi, float(V["re hi %s" % sp]), tol=0.06)
for sp, v in (("dM5", 72.9), ("dM6", 63.9), ("dM7", 40.2), ("dM8", 30.6)):
    ck("recomputed detectable fraction %s (%s)" % (sp, v), v, float(V["re detect %s" % sp]), tol=0.06)

# 第三判据的定义自检：pass_c3 == (peak > 3 rms_norm)，并用该规则重算全部失败率
ck_int("pass_c3 equals peak>3*rms (mismatches, archive)", 0, V["c3 rule mismatch A"])
ck_int("pass_c3 equals peak>3*rms (mismatches, draw)", 0, V["c3 rule mismatch B"])
ck("failure rate by rule (19.4 %)", 19.4, float(V["c3 fail pct by rule"]), tol=0.06)
for zb, v in (("<100", 13.1), ("100-200", 22.2), ("200-400", 36.4), (">400", 36.8)):
    ck("|Z| failure by rule %s (%s)" % (zb, v), v, float(V["c3 zb fail rule %s" % zb]), tol=0.06)
for sp, v in (("dM0", 19.7), ("dM1", 26.4), ("dM2", 17.8), ("dM3", 26.0),
              ("dM4", 19.9), ("dM5", 11.8), ("dM6", 13.2)):
    ck("subtype failure by rule %s (%s)" % (sp, v), v, float(V["c3 sub fail rule %s" % sp]), tol=0.06)
ck("inactive passing by rule (18.9 %)", 18.9, float(V["c3 inactive pass rule pct"]), tol=0.06)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
