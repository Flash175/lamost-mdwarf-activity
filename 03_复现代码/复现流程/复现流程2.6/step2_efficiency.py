# -*- coding: utf-8 -*-
"""step2_efficiency.py — §2.6 的探测效率与两线一致性。

输入：../../数据/hbeta_feasibility.csv（960 条谱，S/N 四档各 240）
输出：输出/step2_efficiency.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
os.makedirs(OUT, exist_ok=True)
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
f = pd.read_csv(os.path.join(D, "hbeta_feasibility.csv"), low_memory=False,
                encoding="utf-8-sig", comment="#")

rows = []
for b in ("5-10", "10-20", "20-40", ">40"):
    g = f[f.snr_bin == b]; ha = g[g.ha_det]
    rows.append(dict(snr_bin=b, n=len(g), n_ha_det=len(ha),
                     p_hb_given_ha_pct=100*ha.hb_det.mean() if len(ha) else np.nan))
hi = f[(f.snr_bin.isin(["20-40", ">40"])) & f.ha_det]
lo = f[(f.snr_bin == "5-10") & f.ha_det]
rows.append(dict(snr_bin=">=20", n=int((f.snr_bin.isin(["20-40",">40"])).sum()), n_ha_det=len(hi),
                 p_hb_given_ha_pct=100*hi.hb_det.mean()))
pd.DataFrame(rows).to_csv(os.path.join(OUT, "step2_efficiency.csv"), index=False, encoding="utf-8-sig")

both = f[f.ha_det & f.hb_det]
x, y = np.log10(both.ew_ha.values), np.log10(both.ew_hb.values)
r = float(np.corrcoef(both.ew_ha, both.ew_hb)[0,1])
slope = float(np.polyfit(x, y, 1)[0])
pd.DataFrame([dict(n_both=len(both), pearson_r=r, loglog_slope=slope)]).to_csv(
    os.path.join(OUT, "step2_consistency.csv"), index=False, encoding="utf-8-sig")

print(pd.DataFrame(rows).round(1).to_string(index=False))
print("  两条线都 >3sigma: %d 条; r = %.3f; log-log slope = %.3f" % (len(both), r, slope))
print("  -> 输出/step2_efficiency.csv, step2_consistency.csv")
