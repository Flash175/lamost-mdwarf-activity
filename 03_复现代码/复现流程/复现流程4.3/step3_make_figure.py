
# -*- coding: utf-8 -*-
"""step3_make_figure.py — 用本流程算出的量重画 §4.3 的瀑布图（Figure 8）。

输入：输出/step1_values.csv
输出：../figures/fig_sixfactor.png|pdf（正文 Figure 8）、输出/图8_瀑布图_数据.csv（画图用的逐条数值）
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "输出")
FIG = os.path.join(os.path.dirname(os.path.dirname(HERE)), "figures")

V = pd.read_csv(os.path.join(OUT, "step1_values.csv"), encoding="utf-8-sig").set_index("item")["value"]
LAB = ("100-200", "200-400")
KEYS = ("F1", "F2", "F3", "F4", "F5", "F6")
NAMES = ("F1 subtype", "F2 $T_{\\rm eff}$", "F3 log g", "F4 [Fe/H]", "F5 log S/N", "F6 log d")

plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "axes.labelsize": 10,
                     "xtick.labelsize": 9.5, "ytick.labelsize": 9.5})
fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.6), gridspec_kw={"wspace": 0.14})
rows = []
for ax, lab in zip(axes, LAB):
    vals = [float(V["contrib %s %s" % (lab, k)]) for k in KEYS]
    obs = float(V["obs decline %s" % lab])
    expl = float(V["explained %s" % lab])
    unex = float(V["unexplained %s" % lab])
    vals.append(unex)
    labs = list(NAMES) + ["unexplained"]
    cum = 0.0
    for i, v in enumerate(vals):
        colour = "#c0392b" if i == 6 else "#2c3e50"
        ax.bar(i, v, bottom=cum, color=colour, width=0.68, zorder=3)
        if abs(v) >= 0.6:
            ytxt = cum + (0.70*v if v > 0 else 0.30*v)
        else:
            ytxt = cum + v + 0.35
        ax.text(i, ytxt, "%+.2f" % v, ha="center", va="center", fontsize=7.9, color="black",
                zorder=5, bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.55))
        rows.append(dict(bin=lab, quantity=labs[i], value_pp=round(v, 4),
                         bar_bottom_pp=round(cum, 4), bar_top_pp=round(cum + v, 4)))
        cum += v
    ax.axhline(obs, color="k", ls="--", lw=1.0, zorder=2)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(labs, rotation=34, ha="right")
    if lab == LAB[0]:
        ax.set_ylabel("contribution to the decline (pp)")
    ax.set_title("%s pc: observed %.2f = %.3f + %.3f pp" % (lab, obs, expl, unex))
    ax.set_ylim(0, obs*1.08)
    ax.grid(axis="y", ls=":", lw=0.4, alpha=0.5, zorder=0)
plt.tight_layout(pad=0.5)
os.makedirs(FIG, exist_ok=True)
fig.savefig(os.path.join(FIG, "fig_sixfactor.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(FIG, "fig_sixfactor.pdf"), bbox_inches="tight")
plt.close(fig)
T = pd.DataFrame(rows)
T.to_csv(os.path.join(OUT, "图8_瀑布图_数据.csv"), index=False, encoding="utf-8-sig")
print(T.head(8).to_string(index=False))
print("  -> figures/fig_sixfactor.png")
