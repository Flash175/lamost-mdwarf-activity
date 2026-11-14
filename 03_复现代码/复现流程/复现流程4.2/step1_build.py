# -*- coding: utf-8 -*-
"""step1_build.py — §4.2 的全部量：占比、逐型倍数与 Wilson、单历元、相关性、C1/C2、L 判据。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv、
      数据/typing_tests.csv、数据/表2.6_Lbol_逐亚型.csv、数据/chi_lbol_by_zbin.csv
输出：输出/step1_values.csv（长表）
"""
import os
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                  usecols=["uid","subclass","ewha","ewha_err","ra","dec","gaia_source_id","gaia_g_mean_mag"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                usecols=["source_id","parallax","parallax_error","ruwe"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)]
ep = m.groupby("uid").size().rename("n_ep")
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
dist = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
u["b"] = sc.galactic.b.deg
u["absz"] = np.abs(dist*np.sin(np.deg2rad(u.b.values)))
u["zgc"] = np.abs(sc.transform_to("galactocentric").z.to(un.pc).value)
MG = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(1000.0/u.parallax.values, 1, None)) + 5
u = u[(MG >= 6) & (u.absz <= 3000)].copy()
print("  2D 样本 %d" % len(u))

vals = []
def put(k, v): vals.append(dict(item=k, value=float(v)))
def wilson(k, n, z=1.96):
    p = k/n; dd = 1 + z*z/n; c = (p + z*z/(2*n))/dd
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/dd
    return 100*(c-h), 100*(c+h)

# 1) dM0 占比
for lo, hi, lab in ((0,100,"in100"), (200,400,"m200"), (400,600,"m400"), (600,1000,"m600")):
    s = u[(u.absz >= lo) & (u.absz < hi)]
    put("dM0 share %s" % lab, 100*(s.subclass == "dM0").mean())

# 2) 逐型倍数与 Wilson 是否不相交
disjoint_all = 1.0
for sp in ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]:
    a = u[(u.absz < 100) & (u.subclass == sp)]; b = u[(u.absz >= 200) & (u.absz < 400) & (u.subclass == sp)]
    fa, fb = a.active.mean(), b.active.mean()
    put("factor %s" % sp, fa/fb)
    la, ha = wilson(int(a.active.sum()), len(a)); lb, hb = wilson(int(b.active.sum()), len(b))
    disjoint_all *= float(hb < la)
put("wilson disjoint all", disjoint_all)

# 3) 单历元
u2 = u.join(ep, on="uid")
one = u2[u2.n_ep == 1]
put("single epoch n", len(one)); put("single epoch pct", 100*len(one)/len(u2))
for sp in ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]:
    a = one[(one.subclass == sp) & (one.absz < 100)]; b = one[(one.subclass == sp) & (one.absz >= 200) & (one.absz < 400)]
    if len(a) >= 30 and len(b) >= 30:
        put("single-epoch factor %s" % sp, a.active.mean()/b.active.mean())
        put("single-epoch nfar %s" % sp, len(b))
facs = [v["value"] for v in vals if v["item"].startswith("single-epoch factor")]
put("single-epoch factor min", min(facs)); put("single-epoch factor max", max(facs))

# 4) 相关性
dist = 1000.0/u.parallax.values
dsinb = dist*np.sin(np.deg2rad(u.b.values))
put("corr zgc dsinb", np.corrcoef(u.zgc.values, np.abs(dsinb))[0,1])

# 5) C1 / C2（typing_tests.csv）
T = pd.read_csv(os.path.join(D, "typing_tests.csv"), **RD)
c1 = T[T.test.str.startswith("C1")].set_index("subtype"); c2 = T[T.test.str.startswith("C2")].set_index("subtype")
for sp in c1.index:
    put("C1 ratio %s" % sp, float(c1.loc[sp, "ratio_test"]))
    put("C2 ratio %s" % sp, float(c2.loc[sp, "ratio_test"]))
put("C1 dM5 nfar", float(c1.loc["dM5", "n_far_test"]))

# 6) L 判据（表2.6 与 chi_lbol_by_zbin）
L = pd.read_csv(os.path.join(D, "表2.6_Lbol_逐亚型.csv"), **RD).set_index("spt")
for sp in ["dM0","dM4","dM5","dM6","dM7","dM8"]:
    put("L pct %s" % sp, float(L.loc[sp, "pct_lB"]))
put("spread EW dM0-dM6", float(L.loc["dM6","pct_EW"] - L.loc["dM0","pct_EW"]))
put("spread L dM0-dM6", float(L.loc["dM6","pct_lB"] - L.loc["dM0","pct_lB"]))
put("spread EW incl dM7", float(L.loc["dM7","pct_EW"] - L.loc["dM0","pct_EW"]))
put("spread L incl dM7", float(L.loc["dM6","pct_lB"] - L.loc["dM7","pct_lB"]))
Z = pd.read_csv(os.path.join(D, "chi_lbol_by_zbin.csv"), **RD).set_index("zbin")
put("L factor in100 over m200", float(Z.loc["<100","pct_lB"])/float(Z.loc["200-400","pct_lB"]))

# 7) 现象数：分箱活动率与两个倍数（新图1_数据_2D分箱.csv）
B = pd.read_csv(os.path.join(D, '新图1_数据_2D分箱.csv'), **RD).set_index('zbin')
for zb in ('<100','100-200','200-400','400-600','600-1000'):
    put('bin f %s' % zb, float(B.loc[zb, 'frac']))
put('factor 3.06', float(B.loc['<100','frac'])/float(B.loc['200-400','frac']))
put('factor 3.58', float(B.loc['<100','frac'])/float(B.loc['400-600','frac']))

# 8) 0.75 A 在各亚型 EW 分布中的分位（2D 逐星）
for sp in ('dM0','dM6'):
    x = u[u.subclass == sp]
    put('pct075 %s' % sp, 100*float(np.mean(x.ewha <= 0.75)))

# 9) 分位判据下的倍数范围（分位判据_下降倍数.csv）
Q = pd.read_csv(os.path.join(D, '分位判据_下降倍数.csv'), **RD).set_index('subtype')
r = Q['ratio_000_to_200_400'].astype(float)
put('quantile ratio min', float(r.min())); put('quantile ratio max', float(r.max()))

# 10) 16 点分解（多因素_渠道对齐_snrr版_降幅口径.csv）
M = pd.read_csv(os.path.join(D, '多因素_渠道对齐_snrr版_降幅口径.csv'), **RD).set_index('bin')
put('drop 200-400', float(M.loc['200-400','drop_pp']))
put('within 200-400', float(M.loc['200-400','within_twochannel']))
put('comp 200-400', float(M.loc['200-400','comp_twochannel']))
put('sum components', float(M.loc['200-400','within_twochannel']) + float(M.loc['200-400','comp_twochannel']))
V = pd.DataFrame(vals).round(4)
V.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")
print(V.to_string(index=False)[:1800])
print("  -> 输出/step1_values.csv（%d 项）" % len(V))
