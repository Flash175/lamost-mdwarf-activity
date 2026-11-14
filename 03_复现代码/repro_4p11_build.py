
# -*- coding: utf-8 -*-
"""step1_build.py — 4.11 节的量：四条通道的可检验状态（RUWE、RV 变星、Ro x [Fe/H]、残余选择）。

输入：数据/金属度_清洁样本_637973.csv（逐星：subclass/absz/active/ruwe/m_h/m_h_err）、
      数据/RV变星_SB1候选_多历元.csv（RV 变星并集清单）、数据/RV变星_剔除检验.csv（多历元计数）、
      数据/双星剔除_亚型降幅.csv（三档 RUWE 的亚型降幅）、数据/终样本_双星检验_637973.csv（四档总降幅）、
      数据/金属度_Ro试点.csv（Ro x [Fe/H] 二分斜率）、数据/图5_逐星样本_严格3340行_⑤.csv（周期样本）
输出：输出/step1_values.csv、数据/表4.11_通道预算.csv
"""
import os, io
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

SUB = ["dM0", "dM1", "dM2", "dM3", "dM4", "dM5", "dM6"]

def load(name, **kw):
    p = os.path.join(D, name)
    lines = [l for l in open(p, encoding="utf-8-sig", errors="replace").read().splitlines()
             if not l.startswith("#")]
    return pd.read_csv(io.StringIO("\n".join(lines)), **kw)

vals = []
def put(k, v):
    vals.append(dict(item=k, value=float(v)))

# ---------- 0) 逐星样本 ----------
cl = load("金属度_清洁样本_637973.csv", low_memory=False,
          usecols=["uid", "subclass", "absz", "active", "ruwe", "m_h", "m_h_err"])
put("n all", len(cl))
c6 = cl[cl.subclass.isin(SUB)].copy()
put("n dm0-dm6", len(c6))

def dec(s):
    """内盘降幅（pp）：0 < |Z| <= 100 pc 减 200 < |Z| <= 400 pc 的活动率。"""
    a = s[(s.absz > 0) & (s.absz <= 100)]
    b = s[(s.absz > 200) & (s.absz <= 400)]
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    return 100.0 * float(a.active.mean() - b.active.mean())

put("decline base all", dec(cl))
put("decline base dm0-dm6", dec(c6))

# ---------- 1) RUWE 通道 ----------
for lim, tag in ((1.2, "12"), (1.1, "11")):
    s = cl[cl.ruwe < lim]
    put("n kept ruwe%s" % tag, len(s))
    put("removed pct ruwe%s" % tag, 100.0 * (1.0 - len(s) / len(cl)))
    put("decline ruwe%s" % tag, dec(s))
    put("delta decline ruwe%s" % tag, dec(s) - dec(cl))

T = load("终样本_双星检验_637973.csv").set_index("subset")
for k, lab in (("baseline", "base"), ("RUWE<1.2", "ruwe12"), ("RUWE<1.1", "ruwe11"), ("RV removed", "rvrem")):
    put("table n %s" % lab, float(T.loc[k, "n"]))
    put("table decline %s" % lab, float(T.loc[k, "decline_pp"]))
put("rv removed count", float(T.loc["baseline", "n"] - T.loc["RV removed", "n"]))
put("rv removed pct", 100.0 * float(T.loc["baseline", "n"] - T.loc["RV removed", "n"]) / float(T.loc["baseline", "n"]))

# 逐亚型（逐星重算 + 归档表互校）
B = load("双星剔除_亚型降幅.csv").set_index("tier")
tiers = {"all(ruwe<1.4)": "base", "ruwe<1.2": "12", "ruwe<1.1": "11"}
for sp in SUB:
    t = c6[c6.subclass == sp]
    put("sub decl base %s" % sp, dec(t))
    for tier, tag in (("ruwe<1.2", "12"), ("ruwe<1.1", "11")):
        put("sub decl ruwe%s %s" % (tag, sp), dec(t[t.ruwe < float(tier[5:])]))
        put("sub ddecl ruwe%s %s" % (tag, sp), dec(t[t.ruwe < float(tier[5:])]) - dec(t))
    for tier, tag in (("all(ruwe<1.4)", "base"), ("ruwe<1.2", "12"), ("ruwe<1.1", "11")):
        put("table sub decl %s %s" % (tag, sp), float(B.loc[tier, "d_%s" % sp]))
for tier, tag in (("all(ruwe<1.4)", "base"), ("ruwe<1.2", "12"), ("ruwe<1.1", "11")):
    put("table decline %s dm0-dm6" % tag, float(B.loc[tier, "D_200-400_pp"]))

# 亚型位移的极值（用于正文的界）
for tag in ("12", "11"):
    d = {sp: abs(float(B.loc["ruwe<%.1f" % (1.2 if tag == "12" else 1.1), "d_%s" % sp]) -
                 float(B.loc["all(ruwe<1.4)", "d_%s" % sp])) for sp in SUB}
    put("max abs ddecl ruwe%s" % tag, max(d.values()))
    put("max abs ddecl ruwe%s dm0-dm4" % tag, max(d[sp] for sp in SUB[:5]))
    put("max abs ddecl ruwe%s excl dm6" % tag, max(d[sp] for sp in SUB[:-1]))

# dM6 的外档星数（正文以「少数星」为由把 dM6 单列为例外）
put("dm6 far n", len(c6[(c6.subclass == "dM6") & (c6.absz > 200) & (c6.absz <= 400)]))
put("dm6 below400 n", len(c6[(c6.subclass == "dM6") & (c6.absz <= 400)]))
for lim, tag in ((1.2, "12"), (1.1, "11")):
    put("dm6 far n ruwe%s" % tag, len(c6[(c6.subclass == "dM6") & (c6.absz > 200) & (c6.absz <= 400) & (c6.ruwe < lim)]))

# ---------- 2) RV 变星通道 ----------
rv = load("RV变星_SB1候选_多历元.csv")
put("rv union rows", len(rv))
multi = int((rv.n_ep >= 2).sum())
put("rv file multi rows", multi)
A = rv[rv.chi2_red > 3]
Bb = rv[(rv.ptp > 10) & (rv.ptp > 5 * rv.err_med)]
put("flag chi2 rows", len(A))
put("flag ptp rows", len(Bb))
RC = load("RV变星_剔除检验.csv").set_index("tag")
put("rv multi-epoch archived", float(RC.loc["multi-epoch only", "n"]))
put("rv baseline archived", float(RC.loc["all (baseline)", "n"]))
put("flag chi2 pct", 100.0 * len(A) / float(RC.loc["multi-epoch only", "n"]))
put("flag ptp pct", 100.0 * len(Bb) / float(RC.loc["multi-epoch only", "n"]))
put("union pct of sample", 100.0 * len(rv) / len(cl))

cl["flag"] = cl.uid.isin(set(rv.uid))
kept = cl[~cl.flag]
put("decline after rv removal", dec(kept))
put("delta decline rv removal", dec(kept) - dec(cl))
for sp in SUB:
    t = c6[c6.subclass == sp]
    put("sub decl rvrem %s" % sp, dec(t[~t.uid.isin(set(rv.uid))]))
    put("sub ddecl rvrem %s" % sp, dec(t[~t.uid.isin(set(rv.uid))]) - dec(t))
put("max abs ddecl rvrem", max(abs(float(dec(c6[(c6.subclass == sp) & (~c6.uid.isin(set(rv.uid)))]) - dec(c6[c6.subclass == sp])))
                               for sp in SUB if len(c6[(c6.subclass == sp) & (~c6.uid.isin(set(rv.uid))) & (c6.absz > 200) & (c6.absz <= 400)]) > 0))
for lab, mask in (("all", cl.absz.notna()), ("far", (cl.absz > 200) & (cl.absz <= 400))):
    s = cl[mask]
    na, nf = int(s.active.sum()), int(s.flag.sum())
    fa = int((s.active & s.flag).sum())
    put("flag pct among active %s" % lab, 100.0 * fa / na)
    put("flag pct among inactive %s" % lab, 100.0 * (nf - fa) / (len(s) - na))
    put("flag pct of all %s" % lab, 100.0 * nf / len(s))

# ---------- 3) Ro x [Fe/H] 通道 ----------
P = load("图5_逐星样本_严格3340行_⑤.csv", usecols=["uid", "spt", "P", "Ro", "RHa", "teff"])
put("period sample n", len(P))
# Ro 段：未饱和段（Ro > 0.132）逐亚型计数、[Fe/H] 两半的对比度与加权均值一致性
Pj = P.merge(cl[["uid", "m_h", "m_h_err", "subclass"]], on="uid", how="left")
Pj = Pj[Pj.m_h.notna() & (Pj.m_h_err.abs() <= 0.25) & Pj.m_h.between(-2.0, 0.6)].copy()
Pj["Ro"] = pd.to_numeric(Pj["Ro"], errors="coerce")
for sp in ("dM3", "dM5", "dM6"):
    put("unsaturated %s" % sp, int((Pj[(Pj.subclass == sp) & (Pj.Ro > 0.132)]).shape[0]))
Rt = load("金属度_Ro试点.csv").set_index("subtype")
for sp in Rt.index:
    put("mh contrast %s" % sp, float(Rt.loc[sp, "mh_high"] - Rt.loc[sp, "mh_low"]))
    cutv = float(Rt.loc[sp, "mh_split"])
    tt = Pj[Pj.subclass == sp]
    put("fit half lo %s" % sp, int(((tt.m_h <= cutv) & (tt.Ro > 0.132)).sum()))
    put("fit half hi %s" % sp, int(((tt.m_h > cutv) & (tt.Ro > 0.132)).sum()))
hcnt = [float(v["value"]) for v in vals if v["item"].startswith("fit half")]
put("fit half min", min(hcnt)); put("fit half max", max(hcnt))
wgt = 1.0 / Rt.dbeta_se ** 2
mn = float((Rt.dbeta * wgt).sum() / wgt.sum())
put("ro chi2 about mean", float((wgt * (Rt.dbeta - mn) ** 2).sum()))
mh = cl[["uid", "m_h", "m_h_err"]].drop_duplicates("uid")
j = P[["uid"]].merge(mh, on="uid", how="left")
use = j.m_h.notna() & (j.m_h_err.abs() <= 0.25) & j.m_h.between(-2.0, 0.6)
put("period sample with mh", int(use.sum()))
put("period sample with mh pct", 100.0 * float(use.mean()))
put("period sample median mh err", float(j.m_h_err[use].median()))

R = load("金属度_Ro试点.csv").set_index("subtype")
w = 1.0 / R.dbeta_se ** 2
put("ro weighted mean dbeta", float((R.dbeta * w).sum() / w.sum()))
put("ro weighted mean dbeta se", float(np.sqrt(1.0 / w.sum())))
put("ro weighted mean dbeta 2sigma", abs(float((R.dbeta * w).sum() / w.sum())) + 2 * float(np.sqrt(1.0 / w.sum())))
put("ro halves min n", float(min(R.n_low.min(), R.n_high.min())))
put("ro halves max n", float(max(R.n_low.max(), R.n_high.max())))
for sp in R.index:
    put("ro n %s" % sp, float(R.loc[sp, "n"]))
    put("ro dbeta %s" % sp, float(R.loc[sp, "dbeta"]))
    put("ro dbeta se %s" % sp, float(R.loc[sp, "dbeta_se"]))

# ---------- 写出 ----------
V = pd.DataFrame(vals)
V.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")

rows = []
def row(ch, item, val, unit, note):
    rows.append(dict(channel=ch, item=item, value=val, unit=unit, note=note))
row("sample", "N (all dM0-dM9)", len(cl), "stars", "637,973 星样本")
row("sample", "N (dM0-dM6)", len(c6), "stars", "亚型分解只用 dM0-dM6")
row("RUWE", "decline baseline", round(dec(cl), 3), "pp", "0<|Z|<=100 vs 200<|Z|<=400")
row("RUWE", "RUWE<1.2 removed", round(100.0 * (1 - len(cl[cl.ruwe < 1.2]) / len(cl)), 2), "per cent", "终样本_双星检验 n = 610,537")
row("RUWE", "decline RUWE<1.2", round(dec(cl[cl.ruwe < 1.2]), 3), "pp", "终样本_双星检验 15.794")
row("RUWE", "RUWE<1.1 removed", round(100.0 * (1 - len(cl[cl.ruwe < 1.1]) / len(cl)), 2), "per cent", "去除 13.6%")
row("RUWE", "decline RUWE<1.1", round(dec(cl[cl.ruwe < 1.1]), 3), "pp", "终样本_双星检验 15.571")
for sp in SUB:
    row("RUWE", "subtype decline %s" % sp, round(float(B.loc["all(ruwe<1.4)", "d_%s" % sp]), 3), "pp", "dM0-dM6 各亚型基线")
row("RV", "multi-epoch with sigma_rv<20", int(RC.loc["multi-epoch only", "n"]), "stars", "20.6 per cent")
row("RV", "flagged chi2_red>3", len(A), "stars", "%.1f per cent of multi-epoch" % (100.0 * len(A) / float(RC.loc["multi-epoch only", "n"])))
row("RV", "flagged ptp>10 & >5 err", len(Bb), "stars", "%.1f per cent of multi-epoch" % (100.0 * len(Bb) / float(RC.loc["multi-epoch only", "n"])))
row("RV", "union flagged", len(rv), "stars", "%.2f per cent of 637,973" % (100.0 * len(rv) / len(cl)))
row("RV", "decline after removal", float(T.loc["RV removed", "decline_pp"]), "pp", "baseline %.3f" % float(T.loc["baseline", "decline_pp"]))
row("Ro x [Fe/H]", "periodic stars with mh", int(use.sum()), "stars", "of %d" % len(P))
row("Ro x [Fe/H]", "median mh error", round(float(j.m_h_err[use].median()), 3), "dex", "")
row("Ro x [Fe/H]", "precision-weighted mean dbeta", round(float((R.dbeta * w).sum() / w.sum()), 4),
    "dex^-1", "+- %.3f" % float(np.sqrt(1.0 / w.sum())))
for sp in R.index:
    row("Ro x [Fe/H]", "dbeta %s" % sp, round(float(R.loc[sp, "dbeta"]), 3), "dex^-1",
        "+- %.3f, %d stars" % (float(R.loc[sp, "dbeta_se"]), int(R.loc[sp, "n"])))
pd.DataFrame(rows).to_csv(os.path.join(D, "表4.11_通道预算.csv"), index=False, encoding="utf-8-sig")
print("quantities:", len(V), "| table rows:", len(rows))
