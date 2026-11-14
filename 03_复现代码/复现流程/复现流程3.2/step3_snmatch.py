# -*- coding: utf-8 -*-
"""step3_snmatch.py — §3.2 第一段的 S/N 匹配：在各 |Z| 档 x 亚型内取 S/N(10-30)，
看近档/远档的活动率是否仍差。含 Wilson 95% 区间与冻结表核对。

输入：输出/step1_2d.csv、数据/sn_matched_by_bin_with_n.csv、数据/coolend_sn_cuts_with_n.csv
输出：输出/step3_snmatch.csv、输出/step3_cells.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据"); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return (100*(c-h), 100*(c+h))

df = pd.read_csv(os.path.join(OUT, "step1_2d.csv"), **RD)
subs = ["dM0","dM1","dM2","dM3","dM4","dM5","dM6","dM7","dM8"]
rows = []
for lb, (zl, zh) in (("<100",(0,100)), ("100-200",(100,200)),
                     ("200-400",(200,400)), (">400",(400,3000))):
    for sp in subs:
        s = df[(df.absz >= zl) & (df.absz < zh) & (df.subclass == sp)]
        if len(s) < 30: continue
        sm = s[(s.snrr >= 10) & (s.snrr <= 30)]
        k, n = int(sm.active.sum()), len(sm)
        lo, hi = wilson(k, n)
        rows.append(dict(zbin=lb, subtype=sp, n_all=len(s), frac_all=100*s.active.mean(),
                         n_match=n, frac_match=100*k/n if n else np.nan, lo=lo, hi=hi))
C = pd.DataFrame(rows)
C.round(4).to_csv(os.path.join(OUT, "step3_snmatch.csv"), index=False, encoding="utf-8-sig")
print(C.round(2).to_string(index=False))

frz = pd.read_csv(os.path.join(D, "sn_matched_by_bin_with_n.csv"), **RD)
j = C.merge(frz, on=["zbin","subtype"], suffixes=("","_frz"))
diff = {c: float(np.nanmax(np.abs(j[c] - j[c+"_frz"]))) for c in
        ["n_all","frac_all","n_match","frac_match","lo","hi"]}
print("\n  与冻结表 sn_matched_by_bin_with_n.csv 的最大绝对差：")
print("   " + "  ".join("%s=%.4f" % (k, v) for k, v in diff.items()))

# 正文引用的四对格子
cells = [("dM0","<100"), ("dM0","200-400"), ("dM3","<100"), ("dM3","200-400"),
         ("dM5","<100"), ("dM5","200-400"), ("dM6","<100"), ("dM6","200-400")]
t = C.set_index(["subtype","zbin"])
out = []
for sp, zb in cells:
    r = t.loc[(sp, zb)]
    out.append(dict(subtype=sp, zbin=zb, n_all=int(r.n_all), frac_all=r.frac_all,
                    n_match=int(r.n_match), frac_match=r.frac_match, lo=r.lo, hi=r.hi,
                    match_share_pct=100*r.n_match/r.n_all))
O = pd.DataFrame(out).round(3)
O.to_csv(os.path.join(OUT, "step3_cells.csv"), index=False, encoding="utf-8-sig")
print()
print(O.to_string(index=False))

ce = pd.read_csv(os.path.join(D, "coolend_sn_cuts_with_n.csv"), **RD)
d5 = ce[(ce.subtype == "dM5") & (ce.cut == "all")].iloc[0]
print("\n  coolend 表核对：dM5 detect %.2f%%, dM6 %.2f%%, dM7 %.2f%%, dM8 %.2f%%"
      % (d5.detect_pct,
         float(ce[(ce.subtype=="dM6")&(ce.cut=="all")].detect_pct.iloc[0]),
         float(ce[(ce.subtype=="dM7")&(ce.cut=="all")].detect_pct.iloc[0]),
         float(ce[(ce.subtype=="dM8")&(ce.cut=="all")].detect_pct.iloc[0])))
print("  -> 输出/step3_snmatch.csv、输出/step3_cells.csv")
