
# -*- coding: utf-8 -*-
"""step3_table5.py — 逐格核对 Table 5（tab:metallicity，正文 _sec_metallicity.tex 第 77-82 行）的 54 格。

printed 矩阵逐字抄自 tex；recomputed 一列由 数据/金属度_清洁样本_637973.csv 逐星重算，
16-84 区间用 300 次自助重采样（种子 42）独立重算，再与归档的 金属度_垂直梯度.csv 互校。
输出：输出/表5逐格核对.csv、输出/step3_summary.csv
"""
import os, io
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

ZB = ["<100", "100-200", "200-400", "400-600", "600-1000", ">1000"]
EDG = {"<100": (0, 100), "100-200": (100, 200), "200-400": (200, 400),
       "400-600": (400, 600), "600-1000": (600, 1000), ">1000": (1000, 3000)}
COLS = ["N", "N_FeH", "cov", "med", "lo", "hi", "sig", "f_all", "f_FeH"]

# 与 _sec_metallicity.tex 第 77-82 行逐字一致
PRINTED = {
    "<100":      [177128, 163856, 92.5, -0.085, -0.086, -0.084, 0.11, 23.81, 23.29],
    "100-200":   [170974, 155789, 91.1, -0.114, -0.114, -0.113, 0.11, 12.90, 12.57],
    "200-400":   [191179, 174332, 91.2, -0.154, -0.155, -0.153, 0.12,  7.78,  7.51],
    "400-600":   [ 71688,  65743, 91.7, -0.201, -0.201, -0.200, 0.13,  6.65,  6.42],
    "600-1000":  [ 25773,  23239, 90.2, -0.203, -0.205, -0.201, 0.13,  7.76,  7.39],
    ">1000":     [  1231,   1070, 86.9, -0.212, -0.219, -0.206, 0.14, 13.24, 12.34],
}

d = pd.read_csv(os.path.join(D, "金属度_清洁样本_637973.csv"), low_memory=False,
                encoding="utf-8-sig", usecols=["absz", "active", "m_h", "m_h_err"])
mh = pd.to_numeric(d.m_h, errors="coerce"); er = pd.to_numeric(d.m_h_err, errors="coerce")
use = (er.abs() <= 0.25) & mh.between(-2.0, 0.6)
z = pd.to_numeric(d.absz, errors="coerce")

rng = np.random.default_rng(42)
rec = {}
for zb in ZB:
    lo, hi = EDG[zb]
    inm = (z > lo) & (z <= hi) if zb == "<100" else (z > lo) & (z <= hi)
    s = inm & use
    v = mh[s].values.astype(float)
    boot = np.array([np.median(rng.choice(v, v.size, replace=True)) for _ in range(300)])
    n_all = int(s.sum()); n_bin = int(inm.sum())
    rec[zb] = [n_bin, n_all, 100.0 * n_all / n_bin, float(np.median(v)),
               float(np.percentile(boot, 16)), float(np.percentile(boot, 84)),
               float(er[s].median()),
               100.0 * float(d.active[inm].mean()), 100.0 * float(d.active[s].mean())]

rows = []
for zb in ZB:
    for j, col in enumerate(COLS):
        p, r = PRINTED[zb][j], rec[zb][j]
        dif = float(r) - float(p)
        if col in ("N", "N_FeH"):
            st = "PASS" if int(r) == int(p) else "STALE"
            tol = 0
        elif col in ("lo", "hi"):
            tol = 0.0015; st = "PASS" if abs(dif) <= tol else ("ROUND" if abs(dif) <= 0.01 else "STALE")
        elif col == "med":
            tol = 0.0006; st = "PASS" if abs(dif) <= tol else ("ROUND" if abs(dif) <= 0.01 else "STALE")
        elif col == "sig":
            tol = 0.005; st = "PASS" if abs(dif) <= tol else "STALE"
        elif col == "cov":
            tol = 0.05; st = "PASS" if abs(dif) <= tol else "STALE"
        else:
            tol = 0.006; st = "PASS" if abs(dif) <= tol else ("ROUND" if abs(dif) <= 0.05 else "STALE")
        rows.append(dict(zbin=zb, column=col, printed=p, recomputed=round(float(r), 4),
                         diff=round(dif, 5), tol=tol, status=st))

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "表5逐格核对.csv"), index=False, encoding="utf-8-sig")
R.to_csv(os.path.join(D, "表4.10_表5逐格.csv"), index=False, encoding="utf-8-sig")

# 与归档的垂直梯度表互校 16-84 区间
G = pd.read_csv(os.path.join(D, "金属度_垂直梯度.csv"), encoding="utf-8-sig").set_index("zbin")
cross = []
for zb in ZB:
    cross.append(dict(zbin=zb, col="lo", archive=float(G.loc[zb, "lo"]), recomputed=rec[zb][4],
                      diff=round(rec[zb][4] - float(G.loc[zb, "lo"]), 5)))
    cross.append(dict(zbin=zb, col="hi", archive=float(G.loc[zb, "hi"]), recomputed=rec[zb][5],
                      diff=round(rec[zb][5] - float(G.loc[zb, "hi"]), 5)))
C = pd.DataFrame(cross)
C["status"] = np.where(C["diff"].abs() <= 0.002, "PASS", "ROUND")

summary = pd.DataFrame([
    dict(item="Table 5 cells checked", value=len(R)),
    dict(item="PASS", value=int((R.status == "PASS").sum())),
    dict(item="ROUND", value=int((R.status == "ROUND").sum())),
    dict(item="STALE", value=int((R.status == "STALE").sum())),
    dict(item="16-84 vs archive PASS", value=int((C.status == "PASS").sum())),
    dict(item="16-84 vs archive total", value=len(C)),
])
summary.to_csv(os.path.join(OUT, "step3_summary.csv"), index=False, encoding="utf-8-sig")
C.to_csv(os.path.join(OUT, "表5区间_与归档互校.csv"), index=False, encoding="utf-8-sig")
C.to_csv(os.path.join(D, "表4.10_表5区间与归档互校.csv"), index=False, encoding="utf-8-sig")

print("cells", len(R), R.status.value_counts().to_dict())
print("16-84 vs archive:", C.status.value_counts().to_dict())
for _, r in R[R.status != "PASS"].iterrows():
    print("  %-6s %-9s %-6s printed=%-10s recomputed=%-10s diff=%s" % (r.status, r.zbin, r.column, r.printed, r.recomputed, r["diff"]))
