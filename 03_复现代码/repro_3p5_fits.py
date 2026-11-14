# -*- coding: utf-8 -*-
"""step2_fits.py — §3.5 与 Figure 2 图注的斜率：面板拟合、|z| 分层拟合、Δβ。

口径（与正文/图注一致）：
  · 断点固定 Ro_sat = 0.132，只用 Ro > 0.132 的未饱和段；
  · 斜率 = log R'Hα 对 log Ro 的 OLS 斜率；
  · 误差 = bootstrap 16–84 百分位半宽（1,000 次，固定随机种子 20261112）；
  · 面板拟合用全部 3,340 颗；|z| 分层用 3,142（可用天体测量）子样本，按银心 |z| 中位分半。

输入：输出/step1_sample.csv
输出：输出/step2_panels.csv、输出/step2_halves.csv、输出/step2_summary.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_sample.csv"), low_memory=False, encoding="utf-8-sig")
RMIN, SEED, NBOOT = 0.132, 20261112, 1000

def fit(df, rmin=RMIN):
    ro = pd.to_numeric(df.Ro, errors="coerce").values
    rh = pd.to_numeric(df.RHa, errors="coerce").values
    s = np.isfinite(ro) & np.isfinite(rh) & (ro > rmin)
    x, y = np.log10(ro[s]), np.log10(rh[s])
    def ols(x, y):
        A = np.vstack([x, np.ones_like(x)]).T
        return np.linalg.lstsq(A, y, rcond=None)[0][0]
    b = ols(x, y)
    rng = np.random.default_rng(SEED); bs = np.empty(NBOOT)
    for i in range(NBOOT):
        k = rng.integers(0, len(x), len(x)); bs[i] = ols(x[k], y[k])
    p16, p84 = np.percentile(bs, [16, 84])
    return dict(beta=b, err=(p84-p16)/2, n=len(x))

rows = []
for p in ["dM0-dM3", "dM4-dM6"]:
    r = fit(df[df.panel == p]); r["panel"] = p; rows.append(r)
r = fit(df); r["panel"] = "all"; rows.append(r)
P = pd.DataFrame(rows)[["panel","n","beta","err"]].round(3)
P.to_csv(os.path.join(OUT, "step2_panels.csv"), index=False, encoding="utf-8-sig")
print("面板拟合（全部 3,340 颗，断点 0.132）："); print(P.to_string(index=False))

U = df[df.usable].copy()
med = float(U.z_gc.median())
lo, hi = U[U.z_gc <= med], U[U.z_gc > med]
rl, rh = fit(lo), fit(hi)
H = pd.DataFrame([
    dict(half="low", zmed=float(lo.z_gc.median()), n_stars=len(lo), **rl),
    dict(half="high", zmed=float(hi.z_gc.median()), n_stars=len(hi), **rh)]).round(3)
H.to_csv(os.path.join(OUT, "step2_halves.csv"), index=False, encoding="utf-8-sig")
print("\n|z| 分层（可用天体测量 %d 颗，中位 %.2f pc）：" % (len(U), med))
print(H.to_string(index=False))

db = rl["beta"] - rh["beta"]
ed = float(np.hypot(rl["err"], rh["err"]))
S = pd.DataFrame([dict(n_usable=len(U), z_median=med, beta_low=rl["beta"], err_low=rl["err"],
                       beta_high=rh["beta"], err_high=rh["err"], dbeta=db, err_dbeta=ed,
                       sigma=abs(db)/ed)]).round(3)
S.to_csv(os.path.join(OUT, "step2_summary.csv"), index=False, encoding="utf-8-sig")
print("\n  Δβ = β_low - β_high = %+.3f +- %.3f  (%.2f sigma)" % (db, ed, abs(db)/ed))
print("  -> 输出/step2_panels.csv、输出/step2_halves.csv、输出/step2_summary.csv")
