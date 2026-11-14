# -*- coding: utf-8 -*-
"""step1_teffaxis.py — §3.7 的清理后 Teff 轴：K 侧 5,047 谱 + M 侧 642,596 星。

注意 M 侧口径（与 s539 一致）：EW 门 → 按 uid 去重 → 左连 Gaia ruwe → ruwe<1.4。
（若先去巨星或先并 Gaia 再去重，会得到 644,632，与正文不符。）
输入：数据/原始输入/lamost/dr13/dr13_mstellar_dM.csv、数据/原始输入/gaia/gaia_dr3_mstars.csv、
      数据/原始输入/processed/dr13/ktype_pilot_clean.csv
输出：输出/step1_teffaxis.csv、输出/step1_mclean.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出")

cat = pd.read_csv(os.path.join(D, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"), low_memory=False,
                  usecols=["uid","subclass","teff","ewha","ewha_err","gaia_source_id","logg"])
cat = cat[(cat.ewha > -20) & (cat.ewha <= 50) & (cat.ewha_err > 0) & (cat.ewha_err <= 5)].copy()
cat["active"] = (cat.ewha >= 0.75) & (cat.ewha > 3*cat.ewha_err)
cat = cat.sort_values("ewha", ascending=False).drop_duplicates("uid")
print("  EW 门 + 去重: %d 颗，活动 %.2f%%（正文 14.35%%）" % (len(cat), 100*cat.active.mean()))
gm = pd.read_csv(os.path.join(D, "原始输入", "gaia", "gaia_dr3_mstars.csv"), low_memory=False,
                 usecols=["source_id","ruwe"])
gm["sid"] = gm.source_id.astype("int64")
cat["sid"] = pd.to_numeric(cat.gaia_source_id, errors="coerce")
cat = cat.merge(gm[["sid","ruwe"]].drop_duplicates("sid"), on="sid", how="left")
catq = cat[cat.ruwe < 1.4].copy()
print("  + ruwe<1.4: %d 颗，活动 %.2f%%（正文 642,596 / 13.53%%）" % (len(catq), 100*catq.active.mean()))

kk = pd.read_csv(os.path.join(D, "原始输入", "processed", "dr13", "ktype_pilot_clean.csv"),
                 low_memory=False, encoding="utf-8-sig", comment="#")
kk["active"] = kk.active.astype(bool)
print("  K 侧清理后 %d 条谱，活动 %.2f%%（正文 5,047 / 1.74%%）" % (len(kk), 100*kk.active.mean()))

edges = np.arange(3100, 5700, 100)
rows = []
for lo in edges[:-1]:
    a = kk[(kk.teff >= lo) & (kk.teff < lo+100)]
    b = catq[(catq.teff >= lo) & (catq.teff < lo+100)]
    if len(a) < 20 and len(b) < 20: continue
    rows.append(dict(teff_ctr=lo+50, n_K=len(a), pct_K=100*a.active.mean() if len(a) else np.nan,
                     n_M=len(b), pct_M=100*b.active.mean() if len(b) else np.nan))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step1_teffaxis.csv"), index=False, encoding="utf-8-sig")
print(); print(T.to_string(index=False))
S = pd.DataFrame([dict(n_dup=len(cat), pct_dup=100*cat.active.mean(),
                       n_clean=len(catq), pct_clean=100*catq.active.mean(),
                       n_M_in_bins=int(T.n_M.sum()), n_M_below3100=int((catq.teff < 3100).sum()),
                       n_M_above4150=int((catq.teff > 4150).sum()),
                       n_bins_K=int((T.n_K > 0).sum()), n_bins_M=int((T.n_M > 0).sum()),
                       k_total=int(T.n_K.sum()), k_pct=100*kk.active.mean())]).round(3)
S.to_csv(os.path.join(OUT, "step1_mclean.csv"), index=False, encoding="utf-8-sig")
print(); print(S.to_string(index=False))
print("  -> 输出/step1_teffaxis.csv、输出/step1_mclean.csv")
