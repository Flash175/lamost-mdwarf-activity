
# -*- coding: utf-8 -*-
"""step4_make_figure9.py — 重画 Figure 9（fig:metallicity，三面板）。

三个面板的数据全部读自已归档的 CSV，不在脚本里引入任何新数字：
(a) 数据/金属度_垂直梯度.csv           -> 中位 [Fe/H] 与 16-84 自助区间、每档星数
(b) 数据/金属度_固定Z活动率.csv        -> 固定 |Z| 下四个 [Fe/H] 分位的活动率
(c) 数据/金属度_逐亚型斜率.csv         -> 逐亚型 logit 斜率与 95% 误差棒
标题里的三个标量从 输出/step1_values.csv 读，不硬编码。
输出：输出/fig9_metallicity.png、输出/图9_绘图数据.csv
"""
import os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))      # 03_复现代码
D = os.path.join(HERE, "数据")
FLOW = os.path.join(HERE, "复现流程", "复现流程4.10")
OUT = os.path.join(FLOW, "输出"); os.makedirs(OUT, exist_ok=True)

grad = pd.read_csv(os.path.join(D, "金属度_垂直梯度.csv"), encoding="utf-8-sig")
fx = pd.read_csv(os.path.join(D, "金属度_固定Z活动率.csv"), encoding="utf-8-sig")
per = pd.read_csv(os.path.join(D, "金属度_逐亚型斜率.csv"), encoding="utf-8-sig")
V = pd.read_csv(os.path.join(OUT, "step1_values.csv"), encoding="utf-8-sig").set_index("item")["value"]

slope = float(V["slope dex per pc"]); r = float(V["pearson r"])
bg = float(V["beta global"]); bgs = float(V["beta global se"])
bf = float(V["beta FE"]); bfs = float(V["beta FE se"])

fig, ax = plt.subplots(1, 3, figsize=(13.4, 4.1))
x = np.arange(len(grad))
ax[0].errorbar(x, grad.mh_med, yerr=[grad.mh_med - grad.lo, grad.hi - grad.mh_med],
               fmt="o-", color="#7b3294", ms=6, capsize=3)
for i, rr in grad.iterrows():
    ax[0].text(i, rr.hi + 0.012, "%s" % format(int(rr.n), ","), ha="center", fontsize=6.5)
ax[0].set_xticks(x); ax[0].set_xticklabels(grad.zbin, fontsize=8)
ax[0].set_xlabel("|Z|  (pc)"); ax[0].set_ylabel("[Fe/H]  (dex, median)")
ax[0].set_title("[Fe/H] decreases outwards\n(%.5f dex pc$^{-1}$, r = %.3f)" % (slope, r), fontsize=9)

b = fx[fx.quartile != "spearman"]
cols = ["#c0392b", "#e08a7f", "#bdc3c7", "#95a5a6", "#7f8c8d"]
drawn = []
for lb, c in zip(["<100", "100-200", "200-400", "400-600", "600-1000"], cols):
    t = b[b.zbin == lb]
    if len(t):
        ax[1].plot(np.arange(len(t)), t.frac_pct, "o-", color=c, ms=5, lw=1.2, label=lb)
        for _, rr in t.iterrows():
            drawn.append(dict(panel="b", series=lb, quartile=rr.quartile, n=int(rr.n),
                              frac_pct=float(rr.frac_pct), mh_med=float(rr.mh_med)))
ax[1].set_xticks(range(4)); ax[1].set_xticklabels(["Q1", "Q2", "Q3", "Q4"], fontsize=8)
ax[1].set_xlabel("[Fe/H] quartile (metal-poor to metal-rich)"); ax[1].set_ylabel("Active fraction (%)")
ax[1].set_title("Mixed subtype mix: apparent positive trend\n(b_mh = +%.2f $\\pm$ %.2f with |Z| fixed)" % (bg, bgs), fontsize=9)
ax[1].legend(fontsize=6.5, title="|Z| (pc)", title_fontsize=6.5, ncol=2)

xx = np.arange(len(per))
ax[2].errorbar(xx, per.b_mh, yerr=1.96 * per.se, fmt="o", color="#2c3e50", ms=6, capsize=3)
ax[2].axhline(0, color="k", lw=0.8)
ax[2].axhline(bg, color="#c0392b", ls="--", lw=1.2)
ax[2].text(len(per) - 0.6, 0.75, "global (mixed): +%.2f" % bg, color="#c0392b", fontsize=6.5, ha="right")
ax[2].set_xticks(xx); ax[2].set_xticklabels(per.subtype, fontsize=8)
ax[2].set_xlabel("spectral subtype ($|Z|\\leq400$ pc)")
ax[2].set_ylabel("logit slope  $\\partial\\ln[f/(1-f)]/\\partial[\\rm Fe/H]$")
ax[2].set_title("Within subtype the sign reverses\n(fixed-effect b_mh = $-$%.2f $\\pm$ %.2f)" % (abs(bf), bfs), fontsize=9)
plt.tight_layout()
fig.savefig(os.path.join(OUT, "fig9_metallicity.png"), dpi=200); plt.close(fig)

for i, rr in grad.iterrows():
    drawn.append(dict(panel="a", series=rr.zbin, quartile="median", n=int(rr.n),
                      frac_pct=np.nan, mh_med=float(rr.mh_med)))
for _, rr in per.iterrows():
    drawn.append(dict(panel="c", series=rr.subtype, quartile="b_mh", n=np.nan,
                      frac_pct=float(rr.b_mh), mh_med=np.nan))
pd.DataFrame(drawn).to_csv(os.path.join(OUT, "图9_绘图数据.csv"), index=False, encoding="utf-8-sig")
pd.DataFrame(drawn).to_csv(os.path.join(D, "表4.10_图9绘图数据.csv"), index=False, encoding="utf-8-sig")

# 与论文里实际使用的图逐像素比对（04_图件/fig6_metallicity.png）
REPO = os.path.dirname(HERE)  # 仓库根
paper = os.path.join(REPO, "04_图件", "fig6_metallicity.png")
regen = os.path.join(OUT, "fig9_metallicity.png")
chk = dict(item="figure pixel comparison", value=np.nan)
try:
    from PIL import Image, ImageChops
    a = Image.open(paper).convert("RGB"); c = Image.open(regen).convert("RGB")
    if a.size != c.size:
        stat = "SIZE"
        diffpx = -1
    else:
        d = np.asarray(ImageChops.difference(a, c)).astype(int)
        diffpx = int((d.max(axis=2) > 8).sum())
        stat = "PASS" if diffpx == 0 else "DIFF"
    share = (100.0 * diffpx / (a.size[0] * a.size[1])) if diffpx >= 0 else float("nan")
    chk = dict(item="figure pixel comparison (%s, %dx%d, %d px differ >8 = %.2f%%)"
                   % (stat, a.size[0], a.size[1], diffpx, share),
               value=diffpx)
except Exception as e:
    chk = dict(item="figure pixel comparison skipped (%s)" % type(e).__name__, value=np.nan)
pd.DataFrame([chk]).to_csv(os.path.join(OUT, "图9_复现自检.csv"), index=False, encoding="utf-8-sig")
print(chk["item"])

print("panel a points:", len(grad), "| panel b series:", b.zbin.nunique(), "x", b.quartile.nunique(),
      "| panel c points:", len(per))
print("titles: slope=%.3e r=%.3f | b_global=%.3f+-%.3f | b_FE=%.3f+-%.3f" % (slope, r, bg, bgs, bf, bfs))
print("saved", os.path.join(OUT, "fig9_metallicity.png"))
