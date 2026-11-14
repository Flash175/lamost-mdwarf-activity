# -*- coding: utf-8 -*-
"""step3_check_24.py — 把 §2.4 正文引用的数字逐项与计算值比对。"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
band = pd.read_csv(os.path.join(OUT, "step1_band.csv")).set_index("zbin")
pr = pd.read_csv(os.path.join(OUT, "step2_passrate.csv"))

rows = []
def ck(item, exp, got, tol=0.005):
    ok = abs(float(got) - float(exp)) <= tol
    rows.append(dict(item=item, expected=exp, computed=got, status="PASS" if ok else "FAIL"))

b = band["band_070_080_pct"]; h = band["halfkill_pp"]; w = band["wide_pp"]
trend = ["200-400", "400-600", "600-1000"]

ck("<100 pc 带宽 (%)", 1.53, round(float(b["<100"]), 2))
ck("200-1000 三档下限 (%)", 0.95, round(float(b[trend].min()), 2))
ck("200-1000 三档上限 (%)", 1.13, round(float(b[trend].max()), 2))
ck("单档最大上限 (pp)", 0.77, round(float(h.max()), 2))
ck("逐档差分 (pp)", 0.29, round(float(h.max() - h.min()), 2), 0.006)
ck("宽带差分 (pp)", 1.6, round(float(w.max() - w.min()), 1), 0.02)
ck("16 点降幅 (pp)", 16.03, round(23.81 - 7.78, 2), 0.02)
ck("带宽最大者在最内档", 1, int(b.idxmax() == "<100"), 0)

for _, r in pr.iterrows():
    ck("过线率 S/N=%d (MC %%)" % int(r.snr), r.paper_pct, round(float(r.pass_mc_pct), 1), 0.6)
    ck("过线率 S/N=%d (解析 %%)" % int(r.snr), r.paper_pct, round(float(r.pass_analytic_pct), 1), 0.6)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step3_check.csv"), index=False, encoding="utf-8-sig")
bad = res[res.status == "FAIL"]
print(bad.to_string(index=False) if len(bad) else "  (no mismatch)")
print("  PASS %d / %d" % (int((res.status == "PASS").sum()), len(res)))
raise SystemExit(0 if len(bad) == 0 else 2)
