
# -*- coding: utf-8 -*-
"""step1_build.py — 4.7 节的量：|Z| 两半的零结果、断点扫描、四估计量对照、tau_c 约定与适用域。

输入：数据/表3.5_Z分层拟合.csv、表3.5_面板拟合.csv、表3.5_Tau约定对照.csv、
      数据/fig5_rss_profile_all.csv（1,101 点断点剖面）、数据/fig5_estimator_2x2.csv、
      数据/图5_逐星样本_严格3340行_⑤.csv、数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv
输出：输出/step1_values.csv、数据/表4.7_零结果与断点.csv
"""
import os, io
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

def load(name):
    p = os.path.join(D, name)
    lines = [l for l in open(p, encoding="utf-8-sig", errors="replace").read().splitlines() if not l.startswith("#")]
    return pd.read_csv(io.StringIO("\n".join(lines)))

vals = []
def put(k, v):
    vals.append(dict(item=k, value=float(v)))

# 1) |Z| 两半
Z = load("表3.5_Z分层拟合.csv").set_index("half")
bl, el = float(Z.loc["low", "beta"]), float(Z.loc["low", "err"])
bh, eh = float(Z.loc["high", "beta"]), float(Z.loc["high", "err"])
db, ed = bh - bl, float(np.hypot(el, eh))
put("zmed low", float(Z.loc["low", "zmed"]))
put("zmed high", float(Z.loc["high", "zmed"]))
put("beta low", bl); put("err low", el)
put("beta high", bh); put("err high", eh)
put("dbeta halves", db); put("err dbeta", ed); put("sigma halves", abs(db)/ed)
put("two sigma bound", abs(db) + 2.0*ed)
put("upper95 bound", abs(db) + 1.96*ed)

# 2) 断点剖面
P = load("fig5_rss_profile_all.csv")
ro, rss = P.Ro_sat.values, P.RSS.values
i0 = int(np.argmin(rss))
put("rss min", float(rss[i0])); put("ro argmin", float(ro[i0]))
for lab, thr in (("0.1pct", 0.001), ("0.5pct", 0.005), ("1pct", 0.010), ("1sig", 1.0/float(rss[i0])), ("4sig", 4.0/float(rss[i0]))):
    m = (rss - rss[i0]) <= thr*float(rss[i0]) if lab.endswith("pct") else (rss - rss[i0]) <= thr*float(rss[i0])
    if m.sum():
        put("range %s lo" % lab, float(ro[m].min())); put("range %s hi" % lab, float(ro[m].max()))
gap = (ro > 0.155) & (ro < 0.160)
put("gap max drss", float((rss[gap] - rss[i0]).max()))
put("gap argmax ro", float(ro[gap][np.argmax(rss[gap] - rss[i0])]))
put("rss at 0.1543", float(rss[np.argmin(np.abs(ro - 0.1543))]))

m1 = (rss - rss[i0]) <= 1.0
ro1 = ro[m1]
d = np.diff(ro1)
brk = np.where(d > 0.002)[0]
put("1sig contiguous hi", float(ro1[brk[0]]) if len(brk) else float(ro1.max()))
put("1sig isolated hi", float(ro1.max()) if len(brk) else float("nan"))

put("rss at 0.1598", float(rss[np.argmin(np.abs(ro - 0.1598))]))

# 3) 面板与估计量 2x2
Pp = load("表3.5_面板拟合.csv").set_index("panel")
put("beta M0M3", float(Pp.loc["dM0-dM3", "beta"])); put("err M0M3", float(Pp.loc["dM0-dM3", "err"]))
put("beta M4M6", float(Pp.loc["dM4-dM6", "beta"])); put("err M4M6", float(Pp.loc["dM4-dM6", "err"]))
E = load("fig5_estimator_2x2.csv")
EA = E[E["口径"].astype(str).str.startswith("A")]
for est, br, key in (("自由截距", "固定 0.132", "fi_fix"), ("自由截距", "自由断点", "fi_free"),
                     ("锚定连续", "固定 0.132", "ac_fix"), ("锚定连续", "自由断点", "ac_free")):
    r = EA[(EA["估计量"] == est) & (EA["断点"] == br)].iloc[0]
    put("2x2 dbeta %s" % key, float(r.dbeta)); put("2x2 sigma %s" % key, float(r.sigma))
    put("2x2 sig %s" % key, float(r.sig))
    if est == "锚定连续" and br == "自由断点":
        put("2x2 bp M0M3 ac_free", float(r.bp_M0M3)); put("2x2 bp M4M6 ac_free", float(r.bp_M4M6))

# 4) tau_c 约定
T = load("表3.5_Tau约定对照.csv").set_index("scheme")
put("tau ratio", float(T.loc["Wright (+0.31)", "tau_med"])/float(T.loc["Li (-0.30)", "tau_med"]))
put("tau ro sat wright", float(T.loc["Wright (+0.31)", "ro_sat_argmin"]))
put("tau dbeta wright", float(T.loc["Wright (+0.31)", "dbeta"]))
put("tau sigma wright", float(T.loc["Wright (+0.31)", "sigma"]))
put("tau panel wright", float(T.loc["Wright (+0.31)", "dbeta_panel"]))
put("tau dbeta li", float(T.loc["Li (-0.30)", "dbeta"]))

# 5) 29 颗缺口星
S = load("图5_逐星样本_严格3340行_⑤.csv")
S = S[np.isfinite(pd.to_numeric(S.M, errors="coerce")) & np.isfinite(pd.to_numeric(S.P, errors="coerce"))]
ro_s = pd.to_numeric(S.P, errors="coerce").values/(10.0**(2.33 - 1.50*pd.to_numeric(S.M, errors="coerce").values - 0.30*pd.to_numeric(S.M, errors="coerce").values**2))
g = S[(ro_s >= 0.155) & (ro_s <= 0.160)]
put("gap n stars", len(g))
for sp in ("dM0", "dM1", "dM2", "dM4"):
    put("gap n %s" % sp, int((g.spt == sp).sum()))

# 6) 适用域：周期子样本 vs 母样本
cat = load(os.path.join("原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")) if False else pd.read_csv(
    os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
    usecols=["uid", "ewha", "ewha_err", "gaia_g_mean_mag"])
core = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)]
u = core.sort_values("ewha", ascending=False).drop_duplicates("uid")
Gp = pd.to_numeric(u.gaia_g_mean_mag, errors="coerce")
sub = u[u.uid.isin(set(S.uid))]
Gs = pd.to_numeric(sub.gaia_g_mean_mag, errors="coerce")
put("parent n uid", len(u))
put("parent ew median", float(u.ewha.median()))
put("sample ew median", float(S.ewha.median()))
put("ew offset", float(S.ewha.median() - u.ewha.median()))
put("parent g median", float(Gp.median()))
put("sample g median", float(Gs.median()))
put("g offset", float(Gs.median() - Gp.median()))
put("frac below 0.75", 100.0*float((S.ewha < 0.75).mean()))

V = pd.DataFrame(vals).round(6)
V.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")
hdr = ["# 4.7 sum: rotation-activity null result (|Z| halves), break scan and its tolerances, four-estimator 2x2,",
       "#        convective-turnover-time convention, the 29-star gap, and the applicability domain of the periodic subsample",
       "# inputs: 表3.5_Z分层拟合.csv, 表3.5_面板拟合.csv, 表3.5_Tau约定对照.csv, fig5_rss_profile_all.csv,",
       "#         fig5_estimator_2x2.csv, 图5_逐星样本_严格3340行_⑤.csv, 原始输入/lamost/dr13/dr13_mstellar_dM.csv",
       "# generated by 复现流程/复现流程4.7/step1_build.py (代码/repro_4p7_build.py)"]
with open(os.path.join(D, "表4.7_零结果与断点.csv"), "w", encoding="utf-8-sig", newline="") as f:
    f.write("\n".join(hdr) + "\n")
    V.to_csv(f, index=False)
print(V.to_string(index=False)[:2200])
print("  -> 输出/step1_values.csv（%d 项）" % len(V))
