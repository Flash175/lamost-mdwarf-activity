# -*- coding: utf-8 -*-
"""step1_sample_chain.py — 复现 Table 1 的九步样本链（§2.1）。"""
import os
import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as un

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
CAT = os.path.join(ROOT, "数据", "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")
GAIA = os.path.join(ROOT, "数据", "原始输入", "gaia", "gaia_dr3_mstars.csv")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

rows = []
def add(step, remaining, unit, removed, note):
    rows.append(dict(step=step, remaining=int(remaining), unit=unit,
                     removed=None if removed is None else int(removed), note=note))
    print("  %-46s %10d %-7s %s" % (note, remaining, unit, "" if removed is None else "removed %d" % removed))

cat = pd.read_csv(CAT, low_memory=False,
                  usecols=["uid","ewha","ewha_err","ra","dec","gaia_source_id","gaia_g_mean_mag"])
n1 = len(cat); add(1, n1, "spectra", None, "DR13 dM catalogue")
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50)]
add(2, len(cat), "spectra", n1-len(cat), "ewha in (-20,50]")
prev = len(cat); cat = cat[(cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
add(3, len(cat), "spectra", prev-len(cat), "ewha_err in (0,5]")
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat["gaia_source_id"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
g = pd.read_csv(GAIA, low_memory=False,
                usecols=["source_id","parallax","parallax_error","ruwe","radial_velocity"])
m = cat.merge(g, left_on="gaia_source_id", right_on="source_id", how="inner")
add(4, len(m), "spectra", len(cat)-len(m), "Gaia DR3 cross-match")
prev = len(m); m = m[(m.parallax > 0.2) & (m.parallax/m.parallax_error > 5) & (m.ruwe < 1.4)].copy()
add(5, len(m), "spectra", prev-len(m), "parallax>0.2 mas, >5sigma, ruwe<1.4")
u = m.sort_values("ewha", ascending=False).drop_duplicates("uid").copy()
add(6, len(u), "stars", None, "one maximum-EW spectrum per star")
d = 1000.0/u.parallax.values
sc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg, distance=d*un.pc, frame="icrs")
u["absz"] = np.abs(d*np.sin(np.deg2rad(sc.galactic.b.deg)))
prev = len(u); u = u[u.absz <= 3000].copy()
add(7, len(u), "stars", prev-len(u), "|Z| <= 3 kpc (heliocentric)")
MG = u.gaia_g_mean_mag.values - 5*np.log10(np.clip(1000.0/u.parallax.values, 1, None)) + 5
prev = len(u); u = u[MG >= 6].copy()
add(8, len(u), "stars", prev-len(u), "bright objects removed, M_G<6  (2D sample)")
gc = SkyCoord(ra=u.ra.values*un.deg, dec=u.dec.values*un.deg,
              distance=(1000.0/u.parallax.values)*un.pc, frame="icrs")
zg = np.abs(gc.transform_to("galactocentric").z.to(un.pc).value)
rv = u.radial_velocity.values
sel = np.isfinite(rv) & (np.abs(rv) < 500) & (zg < 3000)
add(9, int(sel.sum()), "stars", None, "6D: Gaia RV, |RV|<500 km/s, |z_gc|<3 kpc")
pd.DataFrame(rows).to_csv(os.path.join(OUT, "step1_chain.csv"), index=False, encoding="utf-8-sig")
print("  -> 输出/step1_chain.csv")
