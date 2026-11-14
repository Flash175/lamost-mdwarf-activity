# -*- coding: utf-8 -*-
"""step4_break_scan.py — A.6 / Table 10 断点扫描（3,340 颗周期样本）。"""
import os, json
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
df = pd.read_csv(os.path.join(D, "原始输入", "processed", "dr13", "fig5_sample_strict_v4.csv"),
                 low_memory=False, encoding="utf-8-sig", comment="#")
ro = pd.to_numeric(df["Ro"], errors="coerce").values; rh = pd.to_numeric(df["RHa"], errors="coerce").values
m = np.isfinite(ro) & np.isfinite(rh) & (ro > 0) & (rh > 0)
x, y = np.log10(ro[m]), np.log10(rh[m]); n = x.size
o = np.argsort(x); xs, ys = x[o], y[o]
Px = np.concatenate(([0], np.cumsum(xs))); Pxx = np.concatenate(([0], np.cumsum(xs**2)))
Py = np.concatenate(([0], np.cumsum(ys))); Pyy = np.concatenate(([0], np.cumsum(ys**2)))
Pxy = np.concatenate(([0], np.cumsum(xs*ys)))
grid = np.linspace(np.log10(0.004), np.log10(0.5), 1101); ks = np.searchsorted(xs, grid, side="right")
rec = []
for rs, k in zip(10**grid, ks):
    if k < 10 or n-k < 10: continue
    s1 = Pyy[k] - Py[k]**2/k; n2 = n-k
    sxx = Pxx[n]-Pxx[k]-(Px[n]-Px[k])**2/n2
    sxy = Pxy[n]-Pxy[k]-(Px[n]-Px[k])*(Py[n]-Py[k])/n2
    syy = Pyy[n]-Pyy[k]-(Py[n]-Py[k])**2/n2
    rec.append((rs, s1+syy-sxy**2/sxx, sxy/sxx))
prof = pd.DataFrame(rec, columns=["Ro_sat","RSS","beta_unsat"])
prof.to_csv(os.path.join(OUT, "step4_rss_profile.csv"), index=False, encoding="utf-8-sig")
mn = prof["RSS"].min(); best = float(prof.loc[prof["RSS"].idxmin(), "Ro_sat"])
def rng(cond):
    s = prof.loc[cond, "Ro_sat"]; return [float(s.min()), float(s.max())]
res = dict(n=int(n), rss_min=float(mn), best=best, pct_1=rng(prof["RSS"] <= mn*1.01),
           pct_05=rng(prof["RSS"] <= mn*1.005), pct_01=rng(prof["RSS"] <= mn*1.001),
           d1=rng(prof["RSS"] <= mn+1), d4=rng(prof["RSS"] <= mn+4))
json.dump(res, open(os.path.join(OUT, "step4_break.json"), "w", encoding="utf-8"), indent=2)
print("  n=%d  RSS_min=%.2f @ %.4f  1%%=%s" % (n, mn, best, [round(v,3) for v in res["pct_1"]]))
