# -*- coding: utf-8 -*-
"""step1_kpilot.py — §3.6 的 K 侧：未清理 6,111 谱的逐型统计 + 清理后 5,047 谱的 Teff 轴与 S/N 匹配。

输入：数据/原始输入/processed/dr13/ktype_pilot_ew.csv（6,111，未清理）
      数据/原始输入/processed/dr13/ktype_pilot_clean.csv（5,047，清理后）
输出：输出/step1_ktype.csv、输出/step1_teff.csv、输出/step1_kstats.csv
"""
import os
import numpy as np, pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
P = os.path.join(D, "原始输入", "processed", "dr13")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

raw = pd.read_csv(os.path.join(P, "ktype_pilot_ew.csv"), low_memory=False, encoding="utf-8-sig", comment="#")
raw["active"] = raw.active.astype(bool)
cl = pd.read_csv(os.path.join(P, "ktype_pilot_clean.csv"), low_memory=False, encoding="utf-8-sig", comment="#")
cl["active"] = cl.active.astype(bool)
print("  未清理 %d 条；清理后 %d 条" % (len(raw), len(cl)))

def perbin(df, order=("K0","K1","K2","K3","K4","K5","K6","K7")):
    rows = []
    for s in order:
        x = df[df.subclass == s]
        if len(x) == 0:
            rows.append(dict(subtype=s, n=0, n_act=0, act_pct=np.nan)); continue
        rows.append(dict(subtype=s, n=len(x), n_act=int(x.active.sum()), act_pct=100*x.active.mean()))
    return pd.DataFrame(rows)

T = perbin(raw).merge(perbin(cl), on="subtype", suffixes=("_raw", "_clean"))
T.columns = ["subtype","n_raw","n_act_raw","pct_raw","n_clean","n_act_clean","pct_clean"]
T.round(3).to_csv(os.path.join(OUT, "step1_ktype.csv"), index=False, encoding="utf-8-sig")
print(T.round(2).to_string(index=False))

def chisq(df):
    rows = perbin(df); rows = rows[rows.n > 0]
    ks = rows.n_act.values.astype(float); ns = rows.n.values.astype(float)
    p = ks.sum()/ns.sum()
    chi2 = ((ks - ns*p)**2/(ns*p*(1-p))).sum()
    return chi2, stats.chi2.sf(chi2, len(rows)-1), len(rows)-1, 100*p
stats_rows = []
for lab, df in (("raw", raw), ("clean", cl)):
    chi2, pv, dof, pct = chisq(df)
    k1 = df[df.subclass=="K1"].active; k4 = df[df.subclass=="K4"].active
    d = 100*(k4.mean()-k1.mean()); se = 100*np.sqrt(k1.mean()*(1-k1.mean())/len(k1) + k4.mean()*(1-k4.mean())/len(k4))
    stats_rows.append(dict(sample=lab, n=len(df), act_pct=pct, chi2=chi2, dof=dof, p=pv,
                           k1_pct=100*k1.mean(), k4_pct=100*k4.mean(), d_k1k4=d, se_d=se))
S = pd.DataFrame(stats_rows).round(4)
S.to_csv(os.path.join(OUT, "step1_kstats.csv"), index=False, encoding="utf-8-sig")
print(); print(S.to_string(index=False))

rows = []
for lo in range(3900, 5500, 100):
    b = cl[(cl.teff >= lo) & (cl.teff < lo+100)]
    m = b[(b.snrr >= 6) & (b.snrr <= 20)] if len(b) else b
    rows.append(dict(teff_lo=lo, teff_ctr=lo+50, n=len(b), act_pct=100*b.active.mean() if len(b) else np.nan,
                     n_matched=len(m), act_pct_matched=100*m.active.mean() if len(m) else np.nan))
X = pd.DataFrame(rows).round(3)
X.to_csv(os.path.join(OUT, "step1_teff.csv"), index=False, encoding="utf-8-sig")
shift = []
for s_ in ["K0","K1","K2","K3","K4","K5","K7"]:
    x = cl[cl.subclass == s_]; m_ = x[(x.snrr >= 6) & (x.snrr <= 20)]
    if len(m_) == 0: continue
    shift.append(abs(100*m_.active.mean() - 100*x.active.mean()))
S["matched_max_shift_pp"] = max(shift)
S.to_csv(os.path.join(OUT, "step1_kstats.csv"), index=False, encoding="utf-8-sig")
print(); print("  S/N 6-20 最大位移 = %.2f pp（正文 <=1.9）" % max(shift))
print(); print(X[X.n > 0].to_string(index=False))
print("  -> 输出/step1_ktype.csv、输出/step1_kstats.csv、输出/step1_teff.csv")