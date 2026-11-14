# -*- coding: utf-8 -*-
"""step4_threshold.py — §3.2 第二段：3sigma 阈值造成的探测不完备有多大？

两条误差模型各自把 sigma(EW)=K/(S/N) 的噪声加回目录 EW，再按同一判据重判活动：
  (a) K_cat = median(ewha_err * snrr)，用整份 EW 门内目录（s576 流程）；
  (b) K_loc = 2.177 A，用本地 5,999 条谱的 EW 分布（s575 流程）。
报告"同一批星若改在 19.7 的 S/N 下观测，测得活动率会变多少"。

输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/criterion3_by_zb_sample.csv、
      输出/step2_bins.csv、数据/snr_bias_catalog_pipeline.csv
输出：输出/step4_catalog.csv、输出/step4_local.csv、输出/step4_local_bias.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据"); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")
CAT = os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv")

bins = pd.read_csv(os.path.join(OUT, "step2_bins.csv"), **RD)
bins = bins[bins.zb != ">1000"].copy()
BINS = [(r.zb, float(r.snr_med), round(float(r.frac), 2)) for r in bins.itertuples()]
print("  各档输入（|Z| 档，S/N 中位，实测活动率）:", BINS)

# ---------------- (a) 目录误差模型 ----------------
cat = pd.read_csv(CAT, low_memory=False, usecols=["ewha","ewha_err","snrr","subclass"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) &
          (cat.ewha_err <= 5) & (cat.snrr > 0)]
K_cat = float(np.median(cat.ewha_err.values*cat.snrr.values))
act = cat[(cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)].ewha.values
ina = cat[cat.ewha < 0.75].ewha.values
print("\n  (a) EW 门内谱数 %d；K_cat = %.3f A；活动谱 %d，非活动谱 %d"
      % (len(cat), K_cat, len(act), len(ina)))

RNG = np.random.default_rng(20261029); NSIM = 300
def pass_prob_cat(pool, snr, npick=4000):
    e0 = RNG.choice(pool, size=min(npick, len(pool)), replace=False)
    sig = K_cat/snr
    e = e0[:, None] + RNG.normal(0.0, sig, size=(len(e0), NSIM))
    ok = (e >= 0.75) & (e > 3*sig)
    return ok.mean()

rows = []; ref = None
for lb, s, fobs in BINS:
    pa = pass_prob_cat(act, s); pi = pass_prob_cat(ina, s)
    ft = (fobs/100 - pi)/(pa - pi)
    rel = np.nan if ref is None else 100*((ft*pa + (1-ft)*pi) - (ft*ref[0] + (1-ft)*ref[1]))
    if ref is None: ref = (pa, pi)
    rows.append(dict(zbin=lb, snr_med=s, frac_obs=fobs, pass_act=100*pa,
                     pass_inact=100*pi, f_true_solved=100*ft, rel_pp=rel))
A = pd.DataFrame(rows)
A.round(4).to_csv(os.path.join(OUT, "step4_catalog.csv"), index=False, encoding="utf-8-sig")
print(A.round(3).to_string(index=False))

frz = pd.read_csv(os.path.join(D, "snr_bias_catalog_pipeline.csv"), **RD)
j = A.merge(frz, on="zbin", suffixes=("","_frz"))
print("  与冻结表 snr_bias_catalog_pipeline.csv 的最大差：pass_act %.3f pp，pass_inact %.4f pp，K_cat %.3f A"
      % (float(np.max(np.abs(j.pass_act - j.pass_act_frz))),
         float(np.max(np.abs(j.pass_inact - j.pass_inact_frz))),
         abs(K_cat - 0.749)))

# ---------------- (b) 本地误差模型 ----------------
T = pd.read_csv(os.path.join(D, "criterion3_by_zb_sample.csv"), comment="#")
T = T[np.isfinite(T.ew) & np.isfinite(T.err) & (T.err > 0) &
      np.isfinite(T.snrr) & (T.snrr > 0)].copy()
K_loc = float(np.median(T.err.values*T.snrr.values))
print("\n  (b) 本地谱数 %d；K_loc = %.3f A；活动谱 %d，非活动谱 %d"
      % (len(T), K_loc, int(T.active_cat.sum()), int((~T.active_cat).sum())))
act_ew = T[T.active_cat].ew.values; ina_ew = T[~T.active_cat].ew.values

RNG2 = np.random.default_rng(20261029); NSIM2 = 400
SNRS = np.array([40, 30, 20, 15, 12, 10, 8, 6, 5, 4, 3, 2.5, 2.0])
def pass_prob(ew_true, snr):
    sig = K_loc/snr
    e = ew_true[:, None] + RNG2.normal(0.0, sig, size=(len(ew_true), NSIM2))
    ok = (e >= 0.75) & (e > 3*sig)
    return ok.mean()

grid = np.array([0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0])
resp = []
for gg in grid:
    for S in SNRS:
        resp.append(dict(ew_true=gg, snr=float(S), sigma=K_loc/S,
                         pass_pct=100*pass_prob(np.full(200, gg), S)))
R = pd.DataFrame(resp)
R.round(3).to_csv(os.path.join(OUT, "step4_local_response.csv"), index=False, encoding="utf-8-sig")

bias = []
for S in SNRS:
    pa = pass_prob(act_ew, S); pi = pass_prob(ina_ew, S)
    for f_true in (0.05, 0.078, 0.20, 0.50):
        f_meas = f_true*pa + (1-f_true)*pi
        bias.append(dict(snr=float(S), f_true=f_true, pass_act=100*pa,
                         pass_inact=100*pi, f_meas=100*f_meas,
                         bias_pp=100*(f_meas-f_true), ratio=f_meas/f_true))
BB = pd.DataFrame(bias)
BB.round(4).to_csv(os.path.join(OUT, "step4_local_bias.csv"), index=False, encoding="utf-8-sig")

rows = []; ref = None
for lb, s, fobs in BINS:
    pa = pass_prob(act_ew, s); pi = pass_prob(ina_ew, s)
    ft = (fobs/100 - pi)/(pa - pi)
    rel = np.nan if ref is None else 100*((ft*pa + (1-ft)*pi) - (ft*ref[0] + (1-ft)*ref[1]))
    if ref is None: ref = (pa, pi)
    rows.append(dict(zbin=lb, snr_med=s, frac_obs=fobs, pass_act=100*pa,
                     pass_inact=100*pi, f_true_solved=100*ft, rel_pp=rel))
C = pd.DataFrame(rows)
C.round(4).to_csv(os.path.join(OUT, "step4_local.csv"), index=False, encoding="utf-8-sig")
print(C.round(3).to_string(index=False))
print("\n  目录模型 rel 范围 %+.3f ~ %+.3f pp ; 本地模型 rel 范围 %+.3f ~ %+.3f pp"
      % (A.rel_pp.min(), A.rel_pp.max(), C.rel_pp.min(), C.rel_pp.max()))
print("  -> 输出/step4_catalog.csv、输出/step4_local.csv、输出/step4_local_bias.csv、输出/step4_local_response.csv")
