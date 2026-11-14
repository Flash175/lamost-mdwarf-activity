# -*- coding: utf-8 -*-
"""step2_check_41.py — §4.1 的逐项比对（约 40 项）。

样本归属（与正文一致）：逐亚型比例、合并降幅、速率效应、占比、Teff 中位 → 637,973；
S/N 中位、600–1000 pc 星数、dM2/dM5 端点 → 639,711。

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

# 逐亚型：>400 pc 与 200–400 pc
PAIRS = {"dM0": (5.27, 5.05), "dM1": (7.01, 5.60), "dM2": (13.03, 9.25), "dM3": (15.68, 11.62),
         "dM4": (20.66, 14.53), "dM5": (33.85, 36.38), "dM6": (33.33, 52.63)}
for sp, (f4, f2) in PAIRS.items():
    ck("%s pct >400" % sp, f4, float(V["%s pct m400" % sp]))
    ck("%s pct 200-400" % sp, f2, float(V["%s pct m200" % sp]))
ck_int("dM0 n >400 (58,235)", 58235, V["dM0 n m400"])
ck_int("dM0 n 200-400 (65,497)", 65497, V["dM0 n m200"])
ck_int("dM6 n >400 (24)", 24, V["dM6 n m400"])

# 合并降幅、速率效应、权重与占比
ck("merged drop (0.76 pp)", 0.76, float(V["merged drop"]), tol=0.01)
ck("rate effect (+2.41 pp)", 2.41, float(V["rate effect"]), tol=0.01)
ck("dM4 inner weight (27.8%)", 27.8, float(V["dM4 inner weight"]), tol=0.06)
ck("dM4 share of >400 (2.1%)", 2.1, float(V["dM4 share of >400"]), tol=0.06)
ck("dM0 share 200-400 (34.3%)", 34.3, float(V["dM0 share 200-400"]), tol=0.05)
ck("dM0 share 400-600 (54.0%)", 54.0, float(V["dM0 share 400-600"]), tol=0.05)
ck("dM5 share 200-400 (1.37%)", 1.37, float(V["dM5 share 200-400"]), tol=0.01)
ck("dM5 share 400-600 (0.56%)", 0.56, float(V["dM5 share 400-600"]), tol=0.01)

# 近平面参考温度
ck("teff near dM1 (3,670)", 3670, float(V["teff near dM1"]), tol=1.0, hard=5.0)
ck("teff near dM2 (3,517)", 3517, float(V["teff near dM2"]), tol=1.0, hard=5.0)
ck("teff near dM4 (3,351)", 3351, float(V["teff near dM4"]), tol=1.0, hard=5.0)

# dM4 的 Teff 序列（637,973）与 S/N 序列（639,711）
for lab, t in (("in100",3351), ("m200",3421), ("m400",3506), ("m600",3546)):
    ck("teff dM4 %s (%d)" % (lab, t), t, float(V["teffA dM4 %s" % lab]), tol=1.0, hard=5.0)
for lab, s_ in (("in100",16.8), ("m200",8.7), ("m400",5.0), ("m600",2.8)):
    ck("snr dM4 %s (%.1f)" % (lab, s_), s_, float(V["snr dM4 %s" % lab]), tol=0.06)

# dM2 / dM5 端点（639,711）
ck("teff dM2 m600 (3,622)", 3622, float(V["teffB dM2 m600"]), tol=2.0, hard=5.0)
ck("snr dM2 m600 (7.0)", 7.0, float(V["snr dM2 m600"]), tol=0.06)
ck("teff dM5 m400 (3,471)", 3471, float(V["teffB dM5 m400"]), tol=1.0, hard=5.0)
ck("snr dM5 m400 (2.6)", 2.6, float(V["snr dM5 m400"]), tol=0.06)
ck("teff dM5 m600 (3,354)", 3354, float(V["teffB dM5 m600"]), tol=1.0, hard=5.0)
ck("snr dM5 m600 (2.5)", 2.5, float(V["snr dM5 m600"]), tol=0.06)
ck_int("nB dM4 m600 (324)", 324, V["nB dM4 m600"])
ck_int("nB dM5 m600 (143)", 143, V["nB dM5 m600"])

# 逐型贡献（内层权重固定）
for sp, v in (("dM0", 0.043), ("dM1", 0.268), ("dM2", 0.663), ("dM3", 0.261), ("dM4", 1.702), ("dM5", -0.192), ("dM6", -0.321)):
    ck("cell %s" % sp, v, float(V["cell %s" % sp]), tol=0.02)

# 组成归一到内层后的合并值
for lab, v in (("c200", 12.21), ("c400", 13.98), ("c600", 18.05), ("cFar", 14.62), ("cIn", 23.81)):
    ck("compshape %s" % lab, v, float(V["compshape %s" % lab]), tol=0.02)

# 200-400 -> 400-600 的拆分（附录 Table 9）
for lab, v in (("ss rate", 1.333), ("ss comp", -1.613), ("ss inter", -0.850), ("ss total", -1.130)):
    ck(lab, v, float(V[lab]), tol=0.003)
res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step2_check.csv")