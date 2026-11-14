# -*- coding: utf-8 -*-
"""step5_check_34.py — 把 §3.4 的数字逐项与计算值比对。

输入：输出/step2_coolend.csv、输出/step3_matched.csv
输出：输出/step5_check.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")
ce = pd.read_csv(os.path.join(OUT, "step2_coolend.csv"), **RD).set_index("subtype")
ov = pd.read_csv(os.path.join(OUT, "step4_overlap.csv"), **RD).set_index("selection")
ri = pd.read_csv(os.path.join(OUT, "step4_rise.csv"), **RD).iloc[0]
mt = pd.read_csv(os.path.join(OUT, "step3_matched.csv"), **RD)

rows = []
def ck(item, body, got, tol=0.06):
    d = abs(float(got) - float(body))
    rows.append(dict(item=item, body=body, computed=round(float(got), 3),
                     diff=round(float(got)-float(body), 3),
                     status="PASS" if d <= tol else "FAIL"))
def ck_int(item, body, got):
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got)-int(body),
                     status="PASS" if int(got) == int(body) else "FAIL"))

# ---- 星数、活动率、Wilson 区间 ----
N = {"dM5": 22932, "dM6": 3671, "dM7": 483, "dM8": 85}
F = {"dM5": 57.0, "dM6": 72.9, "dM7": 74.7, "dM8": 65.9}
W = {"dM5": (56.3, 57.6), "dM6": (71.4, 74.3), "dM7": (70.7, 78.4), "dM8": (55.3, 75.1)}
for sp in N:
    ck_int("%s n" % sp, N[sp], ce.loc[sp, "n"])
    ck("%s active pct" % sp, F[sp], float(ce.loc[sp, "frac"]))
    ck("%s Wilson lo" % sp, W[sp][0], float(ce.loc[sp, "lo"]))
    ck("%s Wilson hi" % sp, W[sp][1], float(ce.loc[sp, "hi"]))

# ---- 可探测比例（33.0 -> 30.6 已改；dM5/dM6 用于"上升与探测能力相反"）----
ck("detect dM8", 30.6, float(ce.loc["dM8", "detect_pct"]))
ck("detect dM5", 72.9, float(ce.loc["dM5", "detect_pct"]))
ck("detect dM6", 63.9, float(ce.loc["dM6", "detect_pct"]))
ck("activity rises dM5->dM6", True, float(ce.loc["dM6","frac"]) > float(ce.loc["dM5","frac"]), tol=0)
ck("detection falls dM5->dM6", True,
   float(ce.loc["dM6","detect_pct"]) < float(ce.loc["dM5","detect_pct"]), tol=0)

# ---- S/N 匹配 ----
PAIRS = {"3-6": (73.6, 76.0), "6-12": (71.9, 68.9), "6-20": (73.4, 75.3)}
for cut, (v6, v7) in PAIRS.items():
    ck("matched %s dM6" % cut, v6, float(mt[(mt.cut==cut) & (mt.subtype=="dM6")].frac.iloc[0]))
    ck("matched %s dM7" % cut, v7, float(mt[(mt.cut==cut) & (mt.subtype=="dM7")].frac.iloc[0]))
n8 = mt[mt.subtype == "dM8"].set_index("cut").n
ck_int("dM8 matched 3-6", 15, n8.loc["3-6"])
ck_int("dM8 matched 6-12", 6, n8.loc["6-12"])
ck_int("dM8 matched 6-20", 7, n8.loc["6-20"])
ck_int("dM8 matched min (6)", 6, int(n8.min()))
ck_int("dM8 matched max (15)", 15, int(n8.max()))

# ---- 区间重叠与 dM5->dM6 的显著性 ----
def ck_bool(item, got):
    rows.append(dict(item=item, body=True, computed=bool(got), diff="", status="PASS" if bool(got) else "FAIL"))
ck_bool("overlap dM6 vs dM7 (all)", bool(ov.loc["all", "overlap"].iloc[0]) if hasattr(ov.loc["all","overlap"], "__len__") else bool(ov.loc["all", "overlap"]))
for sel in ["S/N 3-6", "S/N 6-12", "S/N 6-20"]:
    v = ov.loc[sel, "overlap"]
    ck_bool("overlap %s" % sel, bool(v.iloc[0]) if hasattr(v, "__len__") else bool(v))
rows.append(dict(item="rise z (dM5->dM6) >= 10", body=10.0, computed=round(float(ri.z_rise), 2),
                 diff=round(float(ri.z_rise)-10.0, 2),
                 status="PASS" if float(ri.z_rise) >= 10 else "FAIL"))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step5_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step5_check.csv")