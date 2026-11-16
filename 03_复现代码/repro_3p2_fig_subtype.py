# -*- coding: utf-8 -*-
"""make_fig_subtype.py — 重画正文 Figure 4（近/远 |Z| x 全样本/S/N 10-30 匹配）。

数据源：数据/sn_matched_by_bin_with_n.csv（与 §3.2 的多项检查同源）
输出：04_图件/fig1_subtype_sn.png（覆盖；旧图存 04_图件/_superseded/）
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))     # 03_复现代码
ROOT = HERE
D = os.path.join(ROOT, "数据")
FIG = os.path.join(os.path.dirname(ROOT), "04_图件")  # 仓库根/04_图件

S = pd.read_csv(os.path.join(D, "sn_matched_by_bin_with_n.csv"), comment="#", encoding="utf-8-sig")
def get(zb, sp):
    return S[(S.zbin == zb) & (S.subtype == sp)].iloc[0]

subs = ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]
w = 0.2
fig, ax = plt.subplots(figsize=(8.6, 4.6))
for i, sp in enumerate(subs):
    n = get("<100", sp); f = get("200-400", sp)
    ax.bar(i-1.5*w, n.frac_all, w, color="0.55", label="near (<100 pc), full" if i == 0 else None)
    ax.bar(i-0.5*w, n.frac_match, w, color="#c0392b", label="near, S/N 10-30" if i == 0 else None)
    ax.bar(i+0.5*w, f.frac_all, w, color="0.8", label="far (200-400 pc), full" if i == 0 else None)
    ax.bar(i+1.5*w, f.frac_match, w, color="#2c7fb8", label="far, S/N 10-30" if i == 0 else None)
    ax.text(i, max(n.frac_all, n.frac_match) + 2, "n$_{far}$=%d" % int(f.n_match), ha="center", fontsize=6, color="0.3")
    print("  %-4s near %.2f->%.2f  far %.2f->%.2f  n_far=%d" % (sp, n.frac_all, n.frac_match, f.frac_all, f.frac_match, int(f.n_match)))
ax.set_xticks(range(len(subs))); ax.set_xticklabels(subs)
ax.set_xlabel("Spectral subtype"); ax.set_ylabel("Active fraction (%)")
ax.set_title("Active fraction by subtype: near vs far |Z|, before/after S/N matching")
ax.legend(fontsize=7, ncol=2, loc="upper left"); ax.set_ylim(0, 85)
fig.tight_layout()
os.makedirs(os.path.join(FIG, "_superseded"), exist_ok=True)
old = os.path.join(FIG, "fig1_subtype_sn.png")
if os.path.exists(old):
    import shutil
    shutil.copyfile(old, os.path.join(FIG, "_superseded", "fig1_subtype_sn_before_2026.png"))
fig.savefig(os.path.join(FIG, "fig1_subtype_sn.png"), dpi=200); plt.close(fig)
print("  -> 04_图件/fig1_subtype_sn.png 已更新（旧图在 04_图件/_superseded/）")
