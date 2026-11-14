
# -*- coding: utf-8 -*-
"""step1_build.py — 4.5 节的量：活动 vs 非活动的 σ(vz) 对比（固定 |Z|）、亚型分层、估计量分歧区间、亮端偏置。

输入：数据/kinematic_diff_errors.csv、kinematic_diff_by_subtype.csv、zmatch_kinematics.csv、
      数据/subtype_kinematics_inner.csv、新图1_数据_6D分箱.csv、
      数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv
输出：输出/step1_values.csv、数据/表4.5_运动学对比.csv、数据/RV覆盖_by_G.csv
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
def put(k, v):
    vals.append(dict(item=k, value=float(v)))

# 1) 三个内盘切片：两种估计量的差与显著性
E = load("kinematic_diff_errors.csv").set_index("slice")
for sl, lab in (("0-25", "0-25"), ("25-50", "25-50"), ("50-75", "50-75")):
    put("std diff %s" % lab, -float(E.loc[sl, "d_std"]))
    put("std se %s" % lab, float(E.loc[sl, "se_std"]))
    put("std z %s" % lab, float(E.loc[sl, "z_std"]))
    put("mad diff %s" % lab, -float(E.loc[sl, "d_mad"]))
    put("mad se %s" % lab, float(E.loc[sl, "se_mad"]))
    put("mad z %s" % lab, float(E.loc[sl, "z_mad"]))
for sl in ("75-100", "100-150", "150-200"):
    put("std z far %s" % sl, float(E.loc[sl, "z_std"]))
    put("mad z far %s" % sl, float(E.loc[sl, "z_mad"]))

# 2) 档内 |Z| 是否对齐（zmatch）
Z = load("zmatch_kinematics.csv").set_index("slice")
for sl in ("0-25", "25-50", "50-75"):
    put("zmed act %s" % sl, float(Z.loc[sl, "zmed_act"]))
    put("zmed inact %s" % sl, float(Z.loc[sl, "zmed_inact"]))
    put("zmed diff %s" % sl, abs(float(Z.loc[sl, "zmed_inact"]) - float(Z.loc[sl, "zmed_act"])))

# 3) 亚型分层
S = load("kinematic_diff_by_subtype.csv")
for grp in ("dM0-dM2", "dM3-dM5"):
    zz = S[S.group == grp].z_std.astype(float).abs()
    put("strata %s zmin" % grp, float(zz.min()))
    put("strata %s zmax" % grp, float(zz.max()))

# 4) 内盘非活动的 σ(vz) 随亚型（检验"活动寿命序列"）
K = load("subtype_kinematics_inner.csv").set_index("subtype")
put("inactive svz min", float(K.svz_inact.min()))
put("inactive svz max", float(K.svz_inact.max()))
put("inactive svz range", float(K.svz_inact.max()) - float(K.svz_inact.min()))
put("active frac dM0", float(K.loc["dM0", "frac_active"]))
put("active frac dM5", float(K.loc["dM5", "frac_active"]))
put("active frac factor", float(K.loc["dM5", "frac_active"])/float(K.loc["dM0", "frac_active"]))

# 5) 6D 子样本规模与 |Z| 分布（新图1 的 6D 分箱表）
S6 = load("新图1_数据_6D分箱.csv")
n6 = float(S6.n.sum())
put("6D n", n6)
put("6D frac inside 400", 100.0*float(S6.loc[S6.zbin != ">400", "n"].sum())/n6)
put("6D n beyond 400", float(S6.loc[S6.zbin == ">400", "n"].iloc[0]))

# 6) 亮端偏置：按 G 统计 Gaia RV 覆盖率（在 637,973 冻结样本上重算）
cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                  usecols=["uid", "subclass", "ewha", "ewha_err", "ra", "dec", "gaia_source_id"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                usecols=["source_id", "parallax", "parallax_error", "ruwe", "phot_g_mean_mag", "radial_velocity", "pmra", "pmdec"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)]
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
dist = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
u["absz"] = np.abs(dist*np.sin(np.deg2rad(sc.galactic.b.deg)))
u["MG"] = u.phot_g_mean_mag.values - 5*np.log10(np.clip(dist, 1, None)) + 5
u = u[(u.MG >= 6) & (u.absz <= 3000)].copy()
G = pd.to_numeric(u.phot_g_mean_mag, errors="coerce")
rv = pd.to_numeric(u.radial_velocity, errors="coerce")
put("n 2D", len(u))
rows = []
for lo, hi in ((12, 14), (14, 16), (16, 18)):
    s = (G >= lo) & (G < hi)
    f = 100.0*rv[s].notna().mean()
    put("rv cover %d-%d" % (lo, hi), f)
    put("n G %d-%d" % (lo, hi), int(s.sum()))
    rows.append(dict(g_lo=lo, g_hi=hi, n=int(s.sum()), n_rv=int(rv[s].notna().sum()), cover_pct=round(f, 4)))
put("rv cover all", 100.0*rv.notna().mean())
RV = pd.DataFrame(rows)
hdr = ["# Gaia DR3 RV 覆盖率（637,973 冻结 2D 样本，按 G 分档）",
       "# 生成：复现流程/复现流程4.5/step1_build.py（代码/repro_4p5_build.py）"]
with open(os.path.join(D, "RV覆盖_by_G.csv"), "w", encoding="utf-8-sig", newline="") as f:
    f.write("\n".join(hdr) + "\n")
    RV.to_csv(f, index=False)


# 7) 两类各自的绝对弥散、样本量，以及误差模型的自洽重算
for sl in ("0-25", "25-50", "50-75", "40-80"):
    r = Z.loc[sl]
    na, ni = float(r["n_act"]), float(r["n_inact"])
    sa, si = float(r["svz_act"]), float(r["svz_inact"])
    ma, mi = float(r["svz_act_mad"]), float(r["svz_inact_mad"])
    put("n act %s" % sl, na)
    put("n inact %s" % sl, ni)
    put("svz act %s" % sl, sa)
    put("svz inact %s" % sl, si)
    put("mad act %s" % sl, ma)
    put("mad inact %s" % sl, mi)
    put("d recomputed %s" % sl, sa - si)
    put("d mad recomputed %s" % sl, ma - mi)
    se_std = np.sqrt(sa*sa/(2.0*(na - 1.0)) + si*si/(2.0*(ni - 1.0)))
    se_mad = np.sqrt((1.858*sa)**2/na + (1.858*si)**2/ni)
    put("se recomputed %s" % sl, se_std)
    put("se mad recomputed %s" % sl, se_mad)
    if sl != "40-80":
        put("se archived %s" % sl, float(E.loc[sl, "se_std"]))
        put("se mad archived %s" % sl, float(E.loc[sl, "se_mad"]))
put("w40-80 z std", float(E.loc["40-80", "z_std"]))
put("w40-80 z mad", float(E.loc["40-80", "z_mad"]))


# 8) 独立重算：从原始数据重建 6D 子样本，自己算 σ(vz)、MAD、档中位与显著性
from astropy.coordinates import Galactocentric
sel = rv.notna() & (rv.abs() < 500)
s6 = u[sel].copy()
rv6 = rv[sel].values
dist6 = 1000.0/s6.parallax.values
c6 = SkyCoord(ra=s6.ra.values*un.deg, dec=s6.dec.values*un.deg, distance=dist6*un.pc,
              pm_ra_cosdec=s6.pmra.values*un.mas/un.yr, pm_dec=s6.pmdec.values*un.mas/un.yr,
              radial_velocity=rv6*un.km/un.s, frame="icrs")
gc = c6.transform_to(Galactocentric())
s6["zgc"] = np.abs(gc.z.to_value(un.pc))
s6["vz"] = gc.v_z.to_value(un.km/un.s)
s6 = s6[s6.zgc <= 3000].copy()
put("re 6D n", len(s6))

def mad(x):
    return 1.4826*np.median(np.abs(x - np.median(x)))

SL = ((0, 25, "0-25"), (25, 50, "25-50"), (50, 75, "50-75"), (40, 80, "40-80"),
      (75, 100, "75-100"), (100, 150, "100-150"), (150, 200, "150-200"))
for lo, hi, lab in SL:
    sl = s6[(s6.zgc >= lo) & (s6.zgc < hi)]
    a = sl[sl.active]
    i2 = sl[~sl.active]
    na, ni = len(a), len(i2)
    sa, si = float(a.vz.std()), float(i2.vz.std())
    ma, mi = mad(a.vz), mad(i2.vz)
    put("re n act %s" % lab, na); put("re n inact %s" % lab, ni)
    put("re svz act %s" % lab, sa); put("re svz inact %s" % lab, si)
    put("re mad act %s" % lab, ma); put("re mad inact %s" % lab, mi)
    put("re d std %s" % lab, sa - si); put("re d mad %s" % lab, ma - mi)
    put("re zmed act %s" % lab, float(a.zgc.median()))
    put("re zmed inact %s" % lab, float(i2.zgc.median()))
    se = float(np.sqrt(sa*sa/(2.0*(na - 1.0)) + si*si/(2.0*(ni - 1.0))))
    sem = float(np.sqrt((1.858*sa)**2/na + (1.858*si)**2/ni))
    put("re se std %s" % lab, se); put("re z std %s" % lab, (sa - si)/se)
    put("re se mad %s" % lab, sem); put("re z mad %s" % lab, (ma - mi)/sem)

V = pd.DataFrame(vals).round(6)
V.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")
hdr2 = ["# 4.5 sum: active-vs-inactive sigma(vz) difference (std and MAD), |Z| matching, subtype strata,",
       "# estimator disagreement beyond 75 pc, inactive dispersion per subtype, 6D subsample and RV coverage",
       "# inputs: kinematic_diff_errors.csv, kinematic_diff_by_subtype.csv, zmatch_kinematics.csv,",
       "#         subtype_kinematics_inner.csv, 新图1_数据_6D分箱.csv + 原始输入（重算 RV 覆盖率）",
       "# generated by 复现流程/复现流程4.5/step1_build.py (代码/repro_4p5_build.py)"]
with open(os.path.join(D, "表4.5_运动学对比.csv"), "w", encoding="utf-8-sig", newline="") as f:
    f.write("\n".join(hdr2) + "\n")
    V.to_csv(f, index=False)
print(V.to_string(index=False)[:1500])
print("  -> 输出/step1_values.csv（%d 项）" % len(V))
