
# -*- coding: utf-8 -*-
"""step1_build.py — 4.6 节的量：冷端（dM6-dM8）的"零结果"与第三判据的噪声底。

输入：数据/表3.4_冷端.csv、表3.4_SN匹配.csv、表3.4_区间检验.csv（§3.4 的三张表）
      数据/criterion3_sample.csv（2,498 条存档谱）、criterion3_by_zb_sample.csv（5,999 条抽样谱）
输出：输出/step1_values.csv、数据/表4.6_冷端与第三判据.csv
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

# 1) 冷端三档的活动率与 Wilson 区间（§3.4 的表）
C = load("表3.4_冷端.csv").set_index("subtype")
for sp in ("dM6", "dM7", "dM8"):
    put("frac %s" % sp, float(C.loc[sp, "frac"]))
    put("lo %s" % sp, float(C.loc[sp, "lo"]))
    put("hi %s" % sp, float(C.loc[sp, "hi"]))
    put("n %s" % sp, float(C.loc[sp, "n"]))
put("detect dM6", float(C.loc["dM6", "detect_pct"]))
put("detect dM8", float(C.loc["dM8", "detect_pct"]))

# 2) S/N 匹配后的 dM6/dM7 与 dM8 的存活谱数
M = load("表3.4_SN匹配.csv")
for cut in ("3-6", "6-12", "6-20"):
    for sp in ("dM6", "dM7"):
        put("matched %s %s" % (cut, sp), float(M[(M.cut == cut) & (M.subtype == sp)].frac.iloc[0]))
    put("matched n %s dM8" % cut, float(M[(M.cut == cut) & (M.subtype == "dM8")].n.iloc[0]))
put("matched n dM8 min", float(M[M.subtype == "dM8"].n.min()))
put("matched n dM8 max", float(M[M.subtype == "dM8"].n.max()))

# 3) 区间重叠检验
O = load("表3.4_区间检验.csv")
for pair, key in (("dM6 vs dM7", "overlap 6-7"), ("dM7 vs dM8", "overlap 7-8")):
    v = O[(O.selection == "all") & (O.pair == pair)].overlap.iloc[0]
    put(key, 1.0 if str(v).strip().lower() in ("true", "1", "1.0") else 0.0)

# 4) 第三判据：2,498 条存档谱（"活动"用目录判据 active_cat）
A = load("criterion3_sample.csv")
A["fail"] = ~A.pass_c3.astype(bool)
act = A[A.active_cat.astype(bool)]
inact = A[~A.active_cat.astype(bool)]
put("c3 n sample", len(A))
put("c3 n active", len(act))
put("c3 n fail", int(act.fail.sum()))
put("c3 fail pct", 100.0*act.fail.mean())
put("c3 n inactive", len(inact))
put("c3 inactive pass pct", 100.0*inact.pass_c3.astype(bool).mean())

# 5) 第三判据：5,999 条抽样谱（|Z| 分档与逐亚型都在活动谱内统计）
B = load("criterion3_by_zb_sample.csv")
B["fail"] = ~B.pass_c3.astype(bool)
put("c3 zb n", len(B))
actb = B[B.active_cat.astype(bool)]
zb = pd.cut(actb.absz, (0, 100, 200, 400, 1e9), labels=["<100", "100-200", "200-400", ">400"], right=False)
g = actb.assign(zb=zb).groupby("zb", observed=True).fail.mean().mul(100.0)
for k in ("<100", "100-200", "200-400", ">400"):
    put("c3 zb fail %s" % k, float(g.loc[k]))
# 逐亚型与"非活动谱也通过"的比例都取自 2,498 条存档谱（A），不是 5,999 抽样（B）
for sp in ("dM0", "dM1", "dM2", "dM3", "dM4", "dM5", "dM6"):
    sub = act[act.subclass == sp]
    put("c3 sub fail %s" % sp, 100.0*sub.fail.mean())


# 6) 独立重算：从 637,973 冻结样本重做冷端的 n、活动率、Wilson 区间与可探测比例
cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                  usecols=["uid", "subclass", "ewha", "ewha_err", "ra", "dec", "gaia_source_id"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
gg = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                 usecols=["source_id", "parallax", "parallax_error", "ruwe", "phot_g_mean_mag"])
mm = cat.merge(gg, left_on="gaia_source_id", right_on="source_id", how="inner")
mm = mm[(mm.parallax > 0.2) & (mm.parallax/mm.parallax_error > 5) & (mm.ruwe < 1.4)]
uu = mm.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
dd = 1000.0/uu.parallax.values
from astropy.coordinates import SkyCoord
import astropy.units as un
sc2 = SkyCoord(ra=uu.ra.values*un.deg, dec=uu.dec.values*un.deg, distance=dd*un.pc, frame="icrs")
uu["absz"] = np.abs(dd*np.sin(np.deg2rad(sc2.galactic.b.deg)))
uu["MG"] = uu.phot_g_mean_mag.values - 5*np.log10(np.clip(dd, 1, None)) + 5
uu = uu[(uu.MG >= 6) & (uu.absz <= 3000)].copy()
uu["detect"] = 3.0*uu.ewha_err <= 0.75
put("re n 2D", len(uu))

def wilson(k, n, z=1.96):
    p = k/n
    d = 1.0 + z*z/n
    c = p + z*z/(2*n)
    h = z*np.sqrt(p*(1 - p)/n + z*z/(4*n*n))
    return 100.0*(c - h)/d, 100.0*(c + h)/d

for sp in ("dM6", "dM7", "dM8"):
    sub = uu[uu.subclass == sp]
    k = int(sub.active.sum()); n = len(sub)
    lo2, hi2 = wilson(k, n)
    put("re n %s" % sp, n)
    put("re frac %s" % sp, 100.0*k/n)
    put("re lo %s" % sp, lo2)
    put("re hi %s" % sp, hi2)
    put("re detect %s" % sp, 100.0*sub.detect.mean())
put("re n dM5", int((uu.subclass == "dM5").sum()))
put("re detect dM5", 100.0*uu[uu.subclass == "dM5"].detect.mean())

# 7) 第三判据的定义自检：pass_c3 应当等于 peak > 3*rms_norm
A["rule"] = A.peak > 3.0*A.rms_norm
B["rule"] = B.peak > 3.0*B.rms_norm
put("c3 rule mismatch A", int((A.rule.astype(bool) != A.pass_c3.astype(bool)).sum()))
put("c3 rule mismatch B", int((B.rule.astype(bool) != B.pass_c3.astype(bool)).sum()))
A["fail_rule"] = ~A.rule.astype(bool)
B["fail_rule"] = ~B.rule.astype(bool)
actr = A[A.active_cat.astype(bool)]
put("c3 fail pct by rule", 100.0*actr.fail_rule.mean())
actr_b = B[B.active_cat.astype(bool)]
zbr = pd.cut(actr_b.absz, (0, 100, 200, 400, 1e9), labels=["<100", "100-200", "200-400", ">400"], right=False)
gr = actr_b.assign(zb=zbr).groupby("zb", observed=True).fail_rule.mean().mul(100.0)
for kk in ("<100", "100-200", "200-400", ">400"):
    put("c3 zb fail rule %s" % kk, float(gr.loc[kk]))
for sp in ("dM0", "dM1", "dM2", "dM3", "dM4", "dM5", "dM6"):
    put("c3 sub fail rule %s" % sp, 100.0*actr[actr.subclass == sp].fail_rule.mean())
put("c3 inactive pass rule pct", 100.0*A[~A.active_cat.astype(bool)].rule.astype(bool).mean())

V = pd.DataFrame(vals).round(6)
V.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")
hdr = ["# 4.6 sum: cool-end null result (dM6-dM8 fractions, Wilson intervals, S/N-matched pairs, detection fraction)",
       "#        and the third criterion (line height above 3x off-line noise) measured on archived spectra",
       "# inputs: 表3.4_冷端.csv, 表3.4_SN匹配.csv, 表3.4_区间检验.csv, criterion3_sample.csv, criterion3_by_zb_sample.csv",
       "# active = catalogue criterion (active_cat); the |Z|-binned and per-subtype rates are computed among ACTIVE spectra",
       "# generated by 复现流程/复现流程4.6/step1_build.py (代码/repro_4p6_build.py)"]
with open(os.path.join(D, "表4.6_冷端与第三判据.csv"), "w", encoding="utf-8-sig", newline="") as f:
    f.write("\n".join(hdr) + "\n")
    V.to_csv(f, index=False)
print(V.to_string(index=False)[:1500])
print("  -> 输出/step1_values.csv（%d 项）" % len(V))
