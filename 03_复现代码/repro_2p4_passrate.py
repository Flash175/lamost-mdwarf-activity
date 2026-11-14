# -*- coding: utf-8 -*-
"""step2_passrate.py — §2.4 的“压线星过线率随 S/N 下降”。

模型（正文 §2.4）：真实 EW 固定在 0.75 Å，测量噪声 sigma_EW = 2.18 / (S/N)；
判据要求测量值 >= 0.75 Å 且 > 3 sigma，因此过线率 = 1 - Phi( max(0.75, 3 sigma) ; 0.75, sigma )。
输出：输出/step2_passrate.csv（解析值 + 蒙特卡洛 + 正文值）
"""
import os
import math
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
os.makedirs(OUT, exist_ok=True)

K = 2.18                     # sigma_EW = K / (S/N)，与正文 §2 的 EW 误差关系一致
TRUE_EW = 0.75
SNRS = [12, 8, 6, 5, 3]
PAPER = {12: 49.9, 8: 40.2, 6: 17.5, 5: 10.3, 3: 2.5}

def phi(z):                  # 标准正态 CDF
    return 0.5 * (1.0 + math.erf(z / np.sqrt(2.0)))

rng = np.random.default_rng(20261112)
rows = []
for sn in SNRS:
    sig = K / sn
    thr = max(TRUE_EW, 3 * sig)
    analytic = 100 * (1 - phi((thr - TRUE_EW) / sig))
    m = TRUE_EW + rng.normal(0, sig, 400_000)
    mc = 100 * np.mean((m >= TRUE_EW) & (m > 3 * sig))
    rows.append(dict(snr=sn, sigma_ew=sig, threshold_used=thr,
                     pass_analytic_pct=analytic, pass_mc_pct=mc, paper_pct=PAPER[sn]))
df = pd.DataFrame(rows)
df.to_csv(os.path.join(OUT, "step2_passrate.csv"), index=False, encoding="utf-8-sig")
print(df.round(2).to_string(index=False))
print("  -> 输出/step2_passrate.csv")
