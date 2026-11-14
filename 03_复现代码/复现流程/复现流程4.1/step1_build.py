# -*- coding: utf-8 -*-
"""step1_build.py — §4.1 的两个样本与全部数量。

样本 A（637,973，去巨星 + |Z|<=3 kpc）：逐亚型比例、合并比例、速率效应、占比、
                                          Teff 中位、近平面参考温度
样本 B（639,711，去巨星前）：S/N 中位、600-1000 pc 的星数、dM2/dM5 端点

输出：输出/step1_values.csv（长表：item, value）
"""
import os
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                  usecols=["uid","subclass","ewha","ewha_err","snrr","teff","ra","dec",
                           "gaia_source_id","gaia_g_mean_mag"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                usecols=["source_id","parallax","parallax_error","ruwe"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)]
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
dist = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
u["absz"] = np.abs(dist*np.sin(np.deg2rad(sc.galactic.b.deg)))
MG = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(1000.0/u.parallax.values, 1, None)) + 5
B = u.copy()                                  # 639,711
A = u[(MG >= 6) & (u.absz <= 3000)].copy()    # 637,973
print("  A(637,973) = %d ; B(639,711) = %d" % (len(A), len(B)))

vals = []
def put(k, v):
    vals.append(dict(item=k, value=float(v)))

# 逐亚型比例与星数（A）
for sp in ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]:
    for lab, lo, hi in (("in100",0,100), ("m200",200,400), ("m400",400,3000)):
        s = A[(A.absz >= lo) & (A.absz < hi) & (A.subclass == sp)]
        put("%s pct %s" % (sp, lab), 100*s.active.mean() if len(s) else np.nan)
        put("%s n %s" % (sp, lab), len(s))

# 合并比例与降幅（A）
f200 = A[(A.absz >= 200) & (A.absz < 400)].active.mean()
f400 = A[A.absz >= 400].active.mean()
put("merged 200-400", 100*f200); put("merged >400", 100*f400); put("merged drop", 100*(f200-f400))

# 速率效应（A，内层权重固定）
inner = A[A.absz < 100]; w = inner.subclass.value_counts(normalize=True)
fa = A[A.absz >= 400].groupby("subclass").active.mean()
fb = A[(A.absz >= 200) & (A.absz < 400)].groupby("subclass").active.mean()
put("rate effect", 100*sum(w[k]*(fa[k]-fb[k]) for k in w.index if k in fa.index and k in fb.index))
put("dM4 inner weight", 100*w.get("dM4", np.nan))
put("dM4 share of >400", 100*(A[A.absz >= 400].subclass == "dM4").mean())
# 逐型贡献与组成归一值
for k in w.index:
    if k in fa.index and k in fb.index:
        put("cell %s" % k, 100*w[k]*(fa[k]-fb[k]))
for lo, hi, lab in ((200,400,"c200"), (400,600,"c400"), (600,1000,"c600"), (400,3000,"cFar"), (0,100,"cIn")):
    s2 = A[(A.absz >= lo) & (A.absz < hi)]
    fr = s2.groupby("subclass").active.mean()
    cm = [k for k in w.index if k in fr.index]
    ww = w[cm]/w[cm].sum()
    put("compshape %s" % lab, 100*float((ww*fr[cm]).sum()))

# 200-400 -> 400-600 的 rate/comp/inter 拆分
AA = A[(A.absz >= 200) & (A.absz < 400)]; BB = A[(A.absz >= 400) & (A.absz < 600)]
wa = AA.subclass.value_counts(normalize=True); wb = BB.subclass.value_counts(normalize=True)
qa = AA.groupby("subclass").active.mean(); qb = BB.groupby("subclass").active.mean()
ks = [k for k in wa.index if k in wb.index and k in qa.index and k in qb.index]
r2 = sum(wa[k]*(qb[k]-qa[k]) for k in ks)
c2 = sum((wb[k]-wa[k])*qa[k] for k in ks)
i2 = sum((wb[k]-wa[k])*(qb[k]-qa[k]) for k in ks)
put("ss rate", 100*r2); put("ss comp", 100*c2); put("ss inter", 100*i2)
put("ss total", 100*(r2+c2+i2))

# 亚型占比（A）
for lab, lo, hi in (("200-400",200,400), ("400-600",400,600)):
    s = A[(A.absz >= lo) & (A.absz < hi)]
    sh = s.subclass.value_counts(normalize=True)
    put("dM0 share %s" % lab, 100*sh.get("dM0", np.nan))
    put("dM5 share %s" % lab, 100*sh.get("dM5", np.nan))

# 近平面参考温度（A）
for sp in ["dM1","dM2","dM4"]:
    s = A[(A.absz < 100) & (A.subclass == sp)]
    put("teff near %s" % sp, s.teff.median())

# Teff 中位（A）与 S/N 中位、星数（B）
for sp in ["dM2","dM4","dM5"]:
    for lab, lo, hi in (("in100",0,100), ("m200",200,400), ("m400",400,600), ("m600",600,1000)):
        a = A[(A.absz >= lo) & (A.absz < hi) & (A.subclass == sp)]
        b = B[(B.absz >= lo) & (B.absz < hi) & (B.subclass == sp)]
        put("teffA %s %s" % (sp, lab), a.teff.median())
        put("teffB %s %s" % (sp, lab), b.teff.median())
        put("snr %s %s" % (sp, lab), b.snrr.median())
        put("nB %s %s" % (sp, lab), len(b))

T = pd.DataFrame(vals).round(4)
T.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False)[:1500])
print("  -> 输出/step1_values.csv（%d 项）" % len(T))