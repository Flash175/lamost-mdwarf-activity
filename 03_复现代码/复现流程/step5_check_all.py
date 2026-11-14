# -*- coding: utf-8 -*-
"""step5_check_all.py — 汇总前四步并与正文/摘要逐项比对。"""
import os, json
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
ch = pd.read_csv(os.path.join(OUT, "step1_chain.csv"))
t2 = pd.read_csv(os.path.join(OUT, "step2_table2.csv"))
d3 = json.load(open(os.path.join(OUT, "step3_decomp.json"), encoding="utf-8"))
d4 = json.load(open(os.path.join(OUT, "step4_break.json"), encoding="utf-8"))
checks = []
def ck(item, exp, got, tol=0.0):
    ok = (abs(float(got)-float(exp)) <= tol) if tol else (int(got) == int(exp))
    checks.append(dict(item=item, expected=exp, computed=got, status="PASS" if ok else "FAIL"))
for i, e in enumerate([1022559,1018728,1016104,997156,842150,639711,639702,637973,82518], 1):
    ck("Table1 step%d" % i, e, int(ch.loc[ch.step == i, "remaining"].iloc[0]))
for _, r in t2.iterrows():
    for k, exp, got, tol in (("n", r.n_paper, r.n_src, 0), ("active%", r.act_paper, r.act_src, 0.005),
                             ("band", r.band_paper, r.band_src, 0.005), ("half", r.half_paper, r.half_src, 0.005),
                             ("wide", r.wide_paper, r.wide_src, 0.005), ("widebound", r.wb_paper, r.wb_src, 0.005)):
        ck("Table2 %s %s" % (r.zbin, k), exp, got, tol)
for b, exp in (("100-200", (10.92, 31.0, 7.52)), ("200-400", (16.09, 38.0, 9.93))):
    ck("Table4 %s observed" % b, exp[0], d3[b]["observed"], 0.01)
    ck("Table4 %s explained%%" % b, exp[1], d3[b]["explained_pct"], 0.5)
    ck("Table4 %s unexplained" % b, exp[2], d3[b]["unexplained_pp"], 0.01)
ss = d3["shift_share_200-400>400-600"]
for k, e in (("rate", 1.333), ("comp", -1.613), ("inter", -0.850), ("total", -1.130)):
    ck("Table9 %s" % k, e, ss[k], 0.002)
ck("Table10 RSS_min", 1132.2, d4["rss_min"], 0.25)
ck("Table10 best", 0.132, d4["best"], 0.001)
ck("Table10 1% low", 0.032, d4["pct_1"][0], 0.001)
ck("Table10 1% high", 0.252, d4["pct_1"][1], 0.002)
ck("Table10 1sigma low", 0.092, d4["d1"][0], 0.001)
ck("Table10 1sigma high", 0.160, d4["d1"][1], 0.001)
res = pd.DataFrame(checks)
res.to_csv(os.path.join(OUT, "step5_check.csv"), index=False, encoding="utf-8-sig")
bad = res[res.status == "FAIL"]
print(bad.to_string(index=False) if len(bad) else "  (no mismatch)")
print("  PASS %d / %d" % (int((res.status == "PASS").sum()), len(res)))
raise SystemExit(0 if len(bad) == 0 else 2)
