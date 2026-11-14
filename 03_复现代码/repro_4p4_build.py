# -*- coding: utf-8 -*-
"""step1_build.py — §4.4 的量：内盘指数拟合与外推、S/N>=15、分类检验、星等限制。

输入：数据/新图1_数据_2D分箱.csv、lnf_binfit.csv、highsn_far_bins.csv、typing_tests.csv、
      数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv
输出：输出/step1_values.csv
"""
import os, io
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

def load(name):
    p = os.path.join(D, name)
    lines = [l for l in open(p, encoding="utf-8-sig", errors="replace").read().splitlines() if not l.startswith("#")]
    return pd.read_csv(io.StringIO("\n".join(lines)))

vals = []
def put(k, v): vals.append(dict(item=k, value=float(v)))

# 1) 内盘指数拟合（对数线性，用 |Z| 中位）
B = load("新图1_数据_2D分箱.csv")
b3 = B[B.zbin.isin(["<100","100-200","200-400"])]
x = b3.absz_med.values.astype(float); y = np.log(b3.frac.values.astype(float)/100)
bb, aa = np.polyfit(x, y, 1)
put("exp a", aa); put("exp b", bb); put("exp H", -1/bb)
put("exp A", np.exp(aa))
# 外推（在 400-600 与 600-1000 的 |Z| 中位）
z400 = float(B[B.zbin=="400-600"].absz_med.iloc[0]); z600 = float(B[B.zbin=="600-1000"].absz_med.iloc[0])
put("z med 400-600", z400); put("z med 600-1000", z600)
put("extrap 400-600", 100*np.exp(aa+bb*z400)); put("extrap 600-1000", 100*np.exp(aa+bb*z600))
put("ratio 400-600", float(B[B.zbin=="400-600"].frac.iloc[0])/(100*np.exp(aa+bb*z400)))
put("ratio 600-1000", float(B[B.zbin=="600-1000"].frac.iloc[0])/(100*np.exp(aa+bb*z600)))

# 2) 二项加权拟合（直接解析 fits 段，该文件有两段结构）
lines = [l for l in open(os.path.join(D, 'lnf_binfit.csv'), encoding='utf-8-sig', errors='replace').read().splitlines() if l.startswith('fits,')]
fits = pd.DataFrame([l.split(',') for l in lines],
                    columns=['section','zbin','a','a_se','b_per_kpc','b_se','H_pc','H_se','chi2'])
for c in ('a','a_se','b_per_kpc','b_se','H_pc','H_se','chi2'):
    fits[c] = fits[c].astype(float)
row = fits[fits.zbin == 'binomial@zmed'].iloc[0]
put('binom a', row.a); put('binom b', row.b_per_kpc); put('binom H', row.H_pc); put('binom H_se', row.H_se)
row2 = fits[fits.zbin == 'unweighted@zmed'].iloc[0]
put('unw H', row2.H_pc)

# 3) S/N>=15
H = load("highsn_far_bins.csv").set_index("zbin")
put("hi 400-600", float(H.loc["400-600","frac_hi"])); put("hi 600-1000", float(H.loc["600-1000","frac_hi"]))

# 4) 分类检验（§4.4 引用 C1/C2）
T = load("typing_tests.csv")
c1 = T[T.test.str.startswith("C1")].set_index("subtype"); c2 = T[T.test.str.startswith("C2")].set_index("subtype")
for sp in ["dM0","dM1","dM2","dM3","dM4","dM5"]:
    put("C1 %s" % sp, float(c1.loc[sp,"ratio_test"])); put("C2 %s" % sp, float(c2.loc[sp,"ratio_test"]))
put("C1 dM5 nfar", float(c1.loc["dM5","n_far_test"]))

# 5) 星等限制（重建样本）
cat = pd.read_csv(os.path.join(D,"原始输入","lamost","dr13","dr13_mstellar_dM.csv"), low_memory=False,
                  usecols=["uid","subclass","ewha","ewha_err","ra","dec","gaia_source_id","gaia_g_mean_mag"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(os.path.join(D,"原始输入","gaia","gaia_dr3_mstars.csv"), low_memory=False,
                usecols=["source_id","parallax","parallax_error","ruwe"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)]
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
dist = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
u["absz"] = np.abs(dist*np.sin(np.deg2rad(sc.galactic.b.deg)))
MG = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(1000.0/u.parallax.values, 1, None)) + 5
u["MG"] = MG
u = u[(MG >= 6) & (u.absz <= 3000)].copy()
for lo, hi, lab in ((400,600,"400-600"), (600,1000,"600-1000")):
    s = u[(u.absz>=lo)&(u.absz<hi)]; q = s[s.gaia_g_mean_mag <= 16.5]
    put("G16.5 n %s" % lab, len(q)); put("G16.5 pct %s" % lab, 100*q.active.mean())
s = u[(u.absz>=600)&(u.absz<1000)]
for mLo, mHi, lab in ((6,7,"MG6-7"), (7,99,"MGge7")):
    q = s[(s.MG>=mLo)&(s.MG<mHi)]
    put("%s n" % lab, len(q)); put("%s pct" % lab, 100*q.active.mean())
    put("%s share" % lab, 100*len(q)/len(s))
for lo, hi, lab in ((0,100,"in100"), (200,400,"m200"), (400,600,"m400"), (600,1000,"m600"), (1000,3000,"far")):
    s = u[(u.absz>=lo)&(u.absz<hi)]
    put("MGge7 shift %s" % lab, 100*s[s.MG>=7].active.mean() - 100*s.active.mean())

# 6) Wilson 区间（外两档，来自 新图1 分箱表）
for zb in ('400-600', '600-1000'):
    r = B[B.zbin == zb].iloc[0]
    put('wilson lo %s' % zb, float(r.lo)); put('wilson hi %s' % zb, float(r.hi))

# 7) 逐亚型外推倍数范围（holdout_prediction.csv）
SP = load('holdout_prediction.csv'); sp = SP[SP.kind == 'subtype']
for zb in ('400-600', '600-1000'):
    g = sp[sp.zbin == zb]
    put('subtype ratio min %s' % zb, float(g.ratio.min()))
    put('subtype ratio max %s' % zb, float(g.ratio.max()))
    put('subtype ratio n %s' % zb, float(len(g)))

# 8) 拟合形式诊断（lnf_fit_diagnostics.csv）
DG = load('lnf_fit_diagnostics.csv').set_index('model')
def pick(df, key):
    return df.loc[df.index[df.index.astype(str).str.startswith(key)][0]]
r1 = pick(DG, 'L1'); r2 = pick(DG, 'L2')
put('L1 rms pp', float(r1.zbin_rms_pp)); put('L1 max pp', float(r1.zbin_max_pp))
put('L2 rms pp', float(r2.zbin_rms_pp)); put('L2 max pp', float(r2.zbin_max_pp))
put('L2 n_par', float(r2.n_par))
put('AUC min', float(DG.auc.min())); put('AUC max', float(DG.auc.max()))
put('tjur max', float(DG.tjur.max()))

# 9) 内盘重拟合的尺度高（lnf_fit_coefs.csv，region=inner）
C = load('lnf_fit_coefs.csv')
inner = C[C.region.astype(str).str.contains('inner')]
def ih(model):
    q = inner[(inner.model.astype(str).str.startswith(model)) & (inner.param == 'b_kpc')].iloc[0]
    return float(q.H_pc), float(q.H_se), float(q.coef)
h1, e1, b1 = ih('M1'); h3, e3, b3 = ih('M3')
a1 = float(inner[(inner.model.astype(str).str.startswith('M1')) & (inner.param == 'a')].coef.iloc[0])
put('inner H exp', h1); put('inner H exp se', e1)
put('inner H subtypeFE', h3); put('inner H subtypeFE se', e3)
put('H rel diff pct', 100*((-1.0/bb) - h1)/h1)

# 10) 逐星平均 vs 档中位代入（lnf_holdout_starlevel.csv + 内盘 M1 系数）
S = load('lnf_holdout_starlevel.csv')
m1 = S[S.model.astype(str).str.startswith('M1')]
for zb in ('400-600', '600-1000'):
    put('starlevel pred %s' % zb, float(m1[m1.zbin == zb].pred.iloc[0]))
    z = float(B[B.zbin == zb].absz_med.iloc[0])
    put('median pred %s' % zb, 100*np.exp(a1 + b1*z/1000.0))

# 11) 独立子样本留出（子样本 A = 全样本）
I = load('lnf_holdout_independent.csv')
ia = I[I.subsample.astype(str).str.startswith('A')]
for zb in ('400-600', '600-1000'):
    put('holdoutA pred %s' % zb, float(ia[ia.zbin == zb].pred.iloc[0]))

# 12) 星等限通道：maglimit_tests.csv 的 M1_G<=16.5 与 M2 的 AIC
MT = load('maglimit_tests.csv')
g165m = MT[MT.test == 'M1_G<=16.5'].set_index('zbin')
for zb in ('400-600', '600-1000'):
    put('maglim n %s' % zb, float(g165m.loc[zb, 'n']))
    put('maglim pct %s' % zb, float(g165m.loc[zb, 'frac']))
m2t = MT[MT.test.astype(str).str.startswith('M2')].set_index('test')
put('AIC dG', float(m2t.loc['M2_+G', 'delta_aic']))
put('AIC dMG', float(m2t.loc['M2_+M_G', 'delta_aic']))
put('AIC dlogd', float(m2t.loc['M2_+ln d', 'delta_aic']))
put('AIC pair dMG', float(m2t.loc['M2_+ln d +M_G', 'delta_aic']))
put('AIC pair dG', float(m2t.loc['M2_+ln d +G', 'delta_aic']))
put('AIC pair spread', float(m2t.loc[['M2_+G +M_G','M2_+ln d +M_G','M2_+ln d +G'], 'delta_aic'].max() - m2t.loc[['M2_+G +M_G','M2_+ln d +M_G','M2_+ln d +G'], 'delta_aic'].min()))

# 13) M_G>=7 后 400-600 档的活动率（地板）
s400 = u[(u.absz >= 400) & (u.absz < 600)]
put('MGge7 pct m400', 100*float(s400[s400.MG >= 7].active.mean()))

# 14) Hβ 3 万谱分层样本的三种口径
HB = load('hbeta30k_per_spectrum.csv')
put('hb30k total', float(len(HB)))
ZB = ['<100', '100-200', '200-400', '400-600', '600-1000', '>1000']
ratio = []
for zb in ZB:
    g = HB[HB.zbin == zb]
    put('hb30k n %s' % zb, float(len(g)))
    put('hb30k cat %s' % zb, 100*float(g.act_cat.mean()))
    put('hb30k locA %s' % zb, 100*float(g.act_ha.mean()))
    put('hb30k locB %s' % zb, 100*float(g.act_hb.mean()))
    ratio.append((zb, float(g.act_hb.mean())/float(g.act_ha.mean())))
inner5 = [r for zb, r in ratio if zb != '>1000']
put('hb30k ratio min', min(inner5)); put('hb30k ratio max', max(inner5))
put('hb30k ratio far', dict(ratio)['>1000'])

# 15) S/N>=20 子样本 <100 pc 的本地 Hβ（正文 29.12）
HA = load('hbeta_activity_per_spectrum.csv')
put('hb20k locB in100', 100*float(HA[HA.zbin == '<100'].act_hb.mean()))

# 16) 外两档的活动率本身（本节的 6.65 / 7.76）
for zb in ("400-600", "600-1000"):
    put("bin f %s" % zb, float(B[B.zbin == zb].frac.iloc[0]))

# 17) 打字漂移能吃掉晚型"超出 1 的部分"多少（C1 -> C2，以及 原始 -> C2）
share_c1c2 = {}
share_raw = {}
for sp in ("dM2", "dM3", "dM4", "dM5"):
    c1v = float(c1.loc[sp, "ratio_test"]); c2v = float(c2.loc[sp, "ratio_test"])
    c0v = float(T[(T.subtype == sp) & (T.test.str.startswith("C1"))].ratio.iloc[0])
    share_c1c2[sp] = 100.0*(c1v - c2v)/(c1v - 1.0)
    share_raw[sp] = 100.0*(c0v - c2v)/(c0v - 1.0)
    put("drift share C1C2 %s" % sp, share_c1c2[sp])
    put("drift share raw %s" % sp, share_raw[sp])
put("drift share C1C2 min", min(share_c1c2.values()))
put("drift share C1C2 max", max(share_c1c2.values()))
put("drift share raw min", min(share_raw.values()))
put("drift share raw max", max(share_raw.values()))


# 18) 阈值不完备带来的加噪位移（正文 "at most 0.16 percentage points"）
TC = load("阈值模型_目录_逐档.csv"); TL = load("阈值模型_本地_逐档.csv")
nc = float(TC.rel_pp.abs().max()); nl = float(TL.rel_pp.abs().max())
put("noise shift catalogue max pp", nc)
put("noise shift local max pp", nl)
put("noise shift both max pp", max(nc, nl))
f0 = float(B[B.zbin == "<100"].frac.iloc[0]); f2 = float(B[B.zbin == "200-400"].frac.iloc[0])
f4 = float(B[B.zbin == "400-600"].frac.iloc[0]); f6 = float(B[B.zbin == "600-1000"].frac.iloc[0])
put("noise vs decline ratio", (f0 - f2)/nc)
put("noise vs outer gap ratio", (f6 - f4)/nc)

# 19) Hbeta/Halpha 逐档比值（正文 0.60--0.69）
for zb in ("<100", "100-200", "200-400", "400-600", "600-1000"):
    g = HB[HB.zbin == zb]
    put("hb ratio %s" % zb, float(g.act_hb.mean())/float(g.act_ha.mean()))

# 20) 指数拟合随外档加入而崩掉（正文 "holds only for |Z| < 400 pc"）
order = ["<100", "100-200", "200-400", "400-600", "600-1000", ">1000"]
Bo = B.set_index("zbin").loc[order].reset_index()
for k in (3, 4, 6):
    s2 = Bo.iloc[:k]
    x2 = s2.absz_med.values.astype(float); y2 = np.log(s2.frac.values.astype(float)/100)
    b2, a2 = np.polyfit(x2, y2, 1)
    pred = 100*np.exp(a2 + b2*x2)
    res2 = pred - s2.frac.values.astype(float)
    put("exp H %dbins" % k, -1.0/b2)
    put("exp rms %dbins" % k, float(np.sqrt((res2**2).mean())))

V = pd.DataFrame(vals).round(6)
V.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")
print(V.to_string(index=False)[:1600])
print("  -> 输出/step1_values.csv（%d 项）" % len(V))
