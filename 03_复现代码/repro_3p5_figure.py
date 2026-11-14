# -*- coding: utf-8 -*-
"""step3_figure.py — 用冻结样本重画 Figure 2（R'Hα vs Ro，两面板，断点固定 0.132）。

输入：输出/step1_sample.csv
输出：输出/fig_rossby_repro.png、输出/step3_figure.csv
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_sample.csv"), low_memory=False, encoding="utf-8-sig")
RMIN = 0.132

def fitline(x, y):
    A = np.vstack([x, np.ones_like(x)]).T
    b, a = np.linalg.lstsq(A, y, rcond=None)[0]
    return a, b

rows = []
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.0), sharey=True)
for ax, pan in zip(axes, ["dM0-dM3", "dM4-dM6"]):
    s = df[df.panel == pan]
    ro = pd.to_numeric(s.Ro, errors="coerce"); rh = pd.to_numeric(s.RHa, errors="coerce")
    m = np.isfinite(ro) & np.isfinite(rh) & (ro > 0) & (rh > 0)
    ax.scatter(ro[m], rh[m], s=4, alpha=0.45, lw=0, color="0.25")
    unsat = m & (ro > RMIN)
    a, b = fitline(np.log10(ro[unsat]), np.log10(rh[unsat]))
    xs = np.logspace(np.log10(RMIN), np.log10(ro[m].max()), 50)
    ax.plot(xs, 10**(a + b*np.log10(xs)), color="crimson", lw=1.6)
    ax.axvline(RMIN, color="0.6", ls=":", lw=1.0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$Ro = P_{\rm rot}/\tau_c$", fontsize=9)
    ax.set_title("%s   N=%d" % (pan, len(s)), fontsize=9)
    ax.tick_params(labelsize=8)
    rows.append(dict(panel=pan, n_total=len(s), n_used=int(m.sum()), n_unsat=int(unsat.sum()),
                     beta=b, intercept=a, ro_min=float(ro[m].min()), ro_max=float(ro[m].max())))
axes[0].set_ylabel(r"$R'_{\rm H\alpha}$", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_rossby_repro.png"), dpi=200); plt.close(fig)
T = pd.DataFrame(rows).round(4)
T.to_csv(os.path.join(OUT, "step3_figure.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))
print("  -> 输出/fig_rossby_repro.png、输出/step3_figure.csv")
