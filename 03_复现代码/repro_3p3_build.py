# -*- coding: utf-8 -*-
"""step1_build.py — §3.3 的两个样本：6D（Gaia RV）与 LAMOST-RV。

6D：2D 样本（637,973）+ |RV|<500 km/s + |z_gc|<3 kpc  -> 82,518 颗
LAMOST-RV：2D 样本 + z!=0, z_err>0, sigma_rv<20, |RV|<500 -> 611,094 颗
速度：astropy 标准 Galactocentric 变换，取 v_z 分量（与正文一致）。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv
输出：输出/step1_6d.csv、输出/step1_lamost.csv
"""
import os
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

C_KMS = 299792.458
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
CAT = os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")
GAI = os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv")

cat = pd.read_csv(CAT, low_memory=False,
                  usecols=["uid","subclass","ewha","ewha_err","z","z_err","ra","dec",
                           "gaia_source_id","gaia_g_mean_mag"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50)]
cat = cat[(cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(GAI, low_memory=False,
                usecols=["source_id","parallax","parallax_error","pmra","pmdec",
                         "radial_velocity","ruwe"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)].copy()
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
dist = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=dist*un.pc, frame="icrs")
u["absz"] = np.abs(dist*np.sin(np.deg2rad(sc.galactic.b.deg)))
u = u[u.absz <= 3000].copy()
MG = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(1000.0/u.parallax.values, 1, None)) + 5
u = u[MG >= 6].copy()
print("  2D sample %d (expect 637,973)" % len(u))

def transform(df, rv):
    c = SkyCoord(ra=df.ra.values*un.deg, dec=df.dec.values*un.deg,
                 distance=(1000.0/df.parallax.values)*un.pc,
                 pm_ra_cosdec=df.pmra.values*un.mas/un.yr, pm_dec=df.pmdec.values*un.mas/un.yr,
                 radial_velocity=rv*un.km/un.s, frame="icrs")
    gc = c.transform_to("galactocentric")
    return gc.velocity.d_z.to(un.km/un.s).value, np.abs(gc.z.to(un.kpc).value)*1000

rv = u.radial_velocity.values
six = u[np.isfinite(rv) & (np.abs(rv) < 500)].copy()
vz, zgc = transform(six, six.radial_velocity.values)
six["vz"], six["zgc"] = vz, zgc
six = six[six.zgc < 3000].copy()
print("  6D sample  %d (expect 82,518); active %d" % (len(six), int(six.active.sum())))
six[["uid","subclass","active","vz","zgc"]].to_csv(
    os.path.join(OUT, "step1_6d.csv"), index=False, encoding="utf-8-sig")

lam = u[(u.z != 0) & (u.z_err > 0) & (C_KMS*u.z_err.abs() < 20) & (np.abs(C_KMS*u.z) < 500)].copy()
lam["rv_lam"] = C_KMS*lam.z.values
lam["sig_rv"] = C_KMS*lam.z_err.abs().values
vzl, _ = transform(lam, lam.rv_lam.values)
lam["vz"] = vzl
vzg, _ = transform(lam, np.nan_to_num(lam.radial_velocity.values, nan=0.0))
vzg[~np.isfinite(lam.radial_velocity.values)] = np.nan
lam["vz_gaia"] = vzg
print("  LAMOST-RV %d (expect 611,094); median sigma_rv %.2f km/s" % (len(lam), float(np.median(lam.sig_rv))))
lam[["uid","subclass","active","absz","sig_rv","rv_lam","radial_velocity","vz","vz_gaia"]].to_csv(
    os.path.join(OUT, "step1_lamost.csv"), index=False, encoding="utf-8-sig")
print("  -> 输出/step1_6d.csv、输出/step1_lamost.csv")
