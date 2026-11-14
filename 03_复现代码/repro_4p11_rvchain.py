
# -*- coding: utf-8 -*-
"""step1b_rv_chain.py — 从原始星表重跑 §4.11 的 RV 变星链（不依赖任何归档中间表）。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv、
      数据/金属度_清洁样本_637973.csv（只用来做样本内交集与降幅/活动率）
输出：输出/step1b_values.csv、数据/表4.11_RV链_原始重算.csv
口径：与 s627_rv_binary_test.py 一致（EW 闸门 -> 并 Gaia -> 视差/ruwe -> 多历元 RV 判据）。
"""
import os, io
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
RAW = os.path.join(D, "原始输入")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
SUB = ["dM0", "dM1", "dM2", "dM3", "dM4", "dM5", "dM6"]
C_KMS = 299792.458

vals = []
def put(k, v):
    vals.append(dict(item=k, value=float(v)))

# ---------- 1) 目录 -> EW 闸门 ----------
d = pd.read_csv(os.path.join(RAW, "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                encoding="utf-8-sig",
                usecols=["uid", "subclass", "ewha", "ewha_err", "z", "z_err", "snrz", "gaia_source_id"])
put("catalogue spectra", len(d))
put("catalogue uid", d.uid.nunique())
d = d[(d.ewha > -20) & (d.ewha <= 50) & (d.ewha_err > 0) & (d.ewha_err <= 5)].copy()
put("ew gate spectra", len(d))
put("ew gate uid", d.uid.nunique())

# ---------- 2) 并 Gaia（视差质量 + ruwe） ----------
g = pd.read_csv(os.path.join(RAW, "gaia", "gaia_dr3_mstars.csv"), low_memory=False, encoding="utf-8-sig",
                usecols=["source_id", "parallax", "parallax_error", "ruwe"])
d = d.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
put("after gaia spectra", len(d))
put("after gaia uid", d.uid.nunique())
d = d[(d.parallax > 0.2) & (d.parallax / d.parallax_error > 5) & (d.ruwe < 1.4)].copy()
put("after quality spectra", len(d))

# ---------- 3) 多历元 RV 判据 ----------
def rv_stats(rv, tag):
    grp = rv.groupby("uid")
    stat = grp.agg(n_ep=("rv", "size"), rv_med=("rv", "median"), rv_min=("rv", "min"),
                   rv_max=("rv", "max"), err_med=("rv_err", "median"))
    stat["ptp"] = stat.rv_max - stat.rv_min
    def chi2_red(t):
        if len(t) < 2:
            return np.nan
        s = np.sqrt(t.rv_err.values ** 2 + 1.0 ** 2)
        r = t.rv.values
        return float(np.sum(((r - r.mean()) / s) ** 2) / (len(r) - 1))
    stat["chi2_red"] = rv.groupby("uid").apply(chi2_red, include_groups=False)
    multi = stat[stat.n_ep >= 2]
    A = multi[multi.chi2_red > 3]
    B = multi[(multi.ptp > 10) & (multi.ptp > 5 * multi.err_med)]
    U = pd.Index(sorted(set(A.index) | set(B.index)))
    put("%s rv spectra" % tag, len(rv))
    put("%s multi-epoch" % tag, len(multi))
    put("%s flag chi2" % tag, len(A))
    put("%s flag ptp" % tag, len(B))
    put("%s union" % tag, len(U))
    put("%s chi2 pct of multi" % tag, 100.0 * len(A) / len(multi))
    put("%s ptp pct of multi" % tag, 100.0 * len(B) / len(multi))
    return multi, A, B, U

rvA = d[(d.z.notna()) & (d.z_err.notna()) & (d.z_err > 0) & (d.snrz >= 5)].copy()
rvA["rv"] = C_KMS * rvA.z; rvA["rv_err"] = C_KMS * rvA.z_err
multiA, flagA_A, flagB_A, UA = rv_stats(rvA, "snrz5")

rvB = rvA[rvA.rv_err < 20].copy()
multiB, flagA_B, flagB_B, UB = rv_stats(rvB, "rvlt20")

# ---------- 4) 与 637,973 星样本求交，重算降幅 ----------
cl = pd.read_csv(os.path.join(D, "金属度_清洁样本_637973.csv"), low_memory=False, encoding="utf-8-sig",
                 usecols=["uid", "subclass", "absz", "active"])
S = set(cl.uid)
put("union snrz5 in sample", len(set(UA) & S))
put("union rvlt20 in sample", len(set(UB) & S))
put("multi snrz5 in sample", len(set(multiA.index) & S))
put("multi rvlt20 in sample", len(set(multiB.index) & S))

def dec(s):
    a = s[(s.absz > 0) & (s.absz <= 100)]; b = s[(s.absz > 200) & (s.absz <= 400)]
    return 100.0 * float(a.active.mean() - b.active.mean())

put("decline baseline", dec(cl))
for tag, U in (("snrz5", UA), ("rvlt20", UB)):
    kept = cl[~cl.uid.isin(set(U))]
    put("decline removed %s" % tag, dec(kept))
    put("delta decline %s" % tag, dec(kept) - dec(cl))
    put("removed pct %s" % tag, 100.0 * len(set(U) & S) / len(cl))

cl["flag"] = cl.uid.isin(set(UA))
na, nf = int(cl.active.sum()), int(cl.flag.sum()); fa = int((cl.active & cl.flag).sum())
put("rate among active", 100.0 * fa / na)
put("rate among inactive", 100.0 * (nf - fa) / (len(cl) - na))
far = (cl.absz > 200) & (cl.absz <= 400)
sf = cl[far].copy(); naf = int(sf.active.sum()); nff = int(sf.flag.sum()); faf = int((sf.active & sf.flag).sum())
put("far rate among active", 100.0 * faf / naf)
put("far rate among inactive", 100.0 * (nff - faf) / (len(sf) - naf))

# 逐亚型：剔除前后降幅（用原始重算的并集）
for sp in SUB:
    t = cl[cl.subclass == sp]
    put("sub decl base %s" % sp, dec(t))
    put("sub decl removed %s" % sp, dec(t[~t.uid.isin(set(UA))]))
    put("sub ddecl %s" % sp, dec(t[~t.uid.isin(set(UA))]) - dec(t))
put("max abs sub ddecl", max(abs(dec(cl[(cl.subclass == sp) & (~cl.uid.isin(set(UA)))]) - dec(cl[cl.subclass == sp]))
                             for sp in SUB))
put("median rv err of multi-epoch", float(rvA[rvA.uid.isin(set(multiA.index))].rv_err.median()))

pd.DataFrame({"uid": sorted(UA)}).to_csv(os.path.join(D, "表4.11_RV变星并集_原始重算.csv"),
                                         index=False, encoding="utf-8-sig")

V = pd.DataFrame(vals)
V.to_csv(os.path.join(OUT, "step1b_values.csv"), index=False, encoding="utf-8-sig")

rows = []
def row(item, val, note):
    rows.append(dict(item=item, value=val, note=note))
row("catalogue spectra (dr13_mstellar_dM.csv 全表)", 1022559, "目录行数")
for k, v in V.set_index("item")["value"].items():
    row(k, v, "")
pd.DataFrame(rows).drop_duplicates("item").to_csv(os.path.join(D, "表4.11_RV链_原始重算.csv"),
                                                  index=False, encoding="utf-8-sig")
print(V.to_string(max_rows=60))
