# -*- coding: utf-8 -*-
"""step4_overlap.py — §3.4 尚未被覆盖的三条断言：区间重叠与 dM5->dM6 的显著性。

① dM6 与 dM7 在总样本上区间重叠；
② 三档 S/N 匹配下 dM6 与 dM7 的 Wilson 区间仍重叠（"overlapping intervals in all three cases"）；
③ dM7 与 dM8 区间重叠；dM5 与 dM6 的上升显著（两比例 z 检验）。

输入：输出/step1_2d.csv
输出：输出/step4_overlap.csv、输出/step4_rise.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
df = pd.read_csv(os.path.join(OUT, "step1_2d.csv"), low_memory=False, encoding="utf-8-sig")

def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*np.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return (100*(c-h), 100*(c+h))

def cell(sub, lo=None, hi=None, closed=False):
    s = df[df.subclass == sub]
    if lo is not None:
        s = s[(s.snrr >= lo) & ((s.snrr <= hi) if closed else (s.snrr < hi))]
    k, n = int(s.active.sum()), len(s)
    a, b = wilson(k, n)
    return dict(subtype=sub, n=n, frac=100*k/n if n else np.nan, lo=a, hi=b)

rows = []
pairs = [("all", "dM6", "dM7", cell("dM6"), cell("dM7")),
         ("all", "dM7", "dM8", cell("dM7"), cell("dM8")),
         ("S/N 3-6", "dM6", "dM7", cell("dM6",3,6), cell("dM7",3,6)),
         ("S/N 6-12", "dM6", "dM7", cell("dM6",6,12), cell("dM7",6,12)),
         ("S/N 6-20", "dM6", "dM7", cell("dM6",6,20,True), cell("dM7",6,20,True))]
for tag, s1, s2, a, b in pairs:
    overlap = (a["lo"] <= b["hi"]) and (b["lo"] <= a["hi"])
    rows.append(dict(selection=tag, pair="%s vs %s" % (s1, s2),
                     f1=a["frac"], lo1=a["lo"], hi1=a["hi"], n1=a["n"],
                     f2=b["frac"], lo2=b["lo"], hi2=b["hi"], n2=b["n"],
                     overlap=bool(overlap)))
T = pd.DataFrame(rows).round(3)
T.to_csv(os.path.join(OUT, "step4_overlap.csv"), index=False, encoding="utf-8-sig")
print(T.to_string(index=False))

# dM5 -> dM6 的两比例 z（合并比例）
a, b = cell("dM5"), cell("dM6")
n1, n2 = a["n"], b["n"]; p1, p2 = a["frac"]/100, b["frac"]/100
p = (p1*n1 + p2*n2)/(n1+n2)
se = np.sqrt(p*(1-p)*(1/n1 + 1/n2))
z = (p2-p1)/se
R = pd.DataFrame([dict(n_dM5=n1, f_dM5=100*p1, n_dM6=n2, f_dM6=100*p2, z_rise=z)]).round(3)
R.to_csv(os.path.join(OUT, "step4_rise.csv"), index=False, encoding="utf-8-sig")
print()
print(R.to_string(index=False))
print("  -> 输出/step4_overlap.csv、输出/step4_rise.csv")
