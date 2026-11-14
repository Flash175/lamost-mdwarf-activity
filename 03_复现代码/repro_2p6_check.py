# -*- coding: utf-8 -*-
"""step4_check_26.py — 把 §2.6 正文的数字逐项与计算值比对。"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
tr = pd.read_csv(os.path.join(OUT, "step1_hbeta_tracks.csv")).set_index("zbin")
ef = pd.read_csv(os.path.join(OUT, "step2_efficiency.csv")).set_index("snr_bin")
co = pd.read_csv(os.path.join(OUT, "step2_consistency.csv")).iloc[0]
lb = pd.read_csv(os.path.join(OUT, "step3_lbol_bins.csv")).set_index("zbin")
ls = pd.read_csv(os.path.join(OUT, "step3_lbol_subtype.csv")).set_index("spt")

rows = []
def ck(item, exp, got, tol=0.005):
    ok = abs(float(got) - float(exp)) <= tol
    rows.append(dict(item=item, expected=exp, computed=got, status="PASS" if ok else "FAIL"))

order = ["<100","100-200","200-400","400-600","600-1000",">1000"]
EXP = {  # (n, cat, ha, hb)
    "<100": (800, 35.38, 35.00, 29.13), "100-200": (800, 18.38, 17.88, 15.75),
    "200-400": (800, 7.38, 7.38, 5.75), "400-600": (800, 5.63, 5.00, 3.38),
    "600-1000": (800, 6.25, 6.13, 3.75), ">1000": (133, 18.05, 18.05, 4.51),
}
for b in order:
    n, c, h1, h2 = EXP[b]
    ck("%s n" % b, n, int(tr.loc[b,"n"]), 0)
    ck("%s cat Ha" % b, c, float(tr.loc[b,"f_cat"]), 0.006)
    ck("%s local Ha" % b, h1, float(tr.loc[b,"f_ha"]), 0.006)
    ck("%s local Hb" % b, h2, float(tr.loc[b,"f_hb"]), 0.006)

ck("ratio 400-600/600-1000 (cat)", 0.90, round(float(tr.loc["400-600","f_cat"]/tr.loc["600-1000","f_cat"]),2))
ck("ratio 400-600/600-1000 (Ha)", 0.82, round(float(tr.loc["400-600","f_ha"]/tr.loc["600-1000","f_ha"]),2))
ck("ratio 400-600/600-1000 (Hb)", 0.90, round(float(tr.loc["400-600","f_hb"]/tr.loc["600-1000","f_hb"]),2))
ck("Hb deficit near plane (pp)", 6.3, float(tr.loc["<100","f_cat"]-tr.loc["<100","f_hb"]), 0.07)
ck("Hb deficit outermost (pp)", 13.5, float(tr.loc[">1000","f_cat"]-tr.loc[">1000","f_hb"]), 0.07)

ck("efficiency 20-40 pc (%)", 75, round(float(ef.loc["20-40","p_hb_given_ha_pct"]),0), 0.5)
ck("efficiency >40 pc (%)", 85, round(float(ef.loc[">40","p_hb_given_ha_pct"]),0), 0.5)
ck("efficiency 5-10 pc (%)", 44, round(float(ef.loc["5-10","p_hb_given_ha_pct"]),0), 0.5)
ck("both lines >3sigma (n)", 327, int(co["n_both"]), 0)
ck("Pearson r", 0.77, round(float(co["pearson_r"]),2))
ck("log-log slope", 0.98, round(float(co["loglog_slope"]),2))

for b, exp in (("<100",18.58),("100-200",9.71),("200-400",5.78),("400-600",5.48),("600-1000",6.76),(">1000",12.27)):
    ck(f"L/Lbol {b} (%)", exp, round(float(lb.loc[b,"pct_lB"]),2))
ck("L/Lbol decline factor", 2.75, round(float(lb.loc["<100","pct_lB"]/lb.loc["600-1000","pct_lB"]),2))
ck("EW decline factor", 3.06, float(lb.loc["<100","pct_EW"]/lb.loc["600-1000","pct_EW"]), 0.01)
ck("EW at 1e-4, dM0 (A)", 0.86, round(float(ls.loc["dM0","ew_at_1e4"]),2))
ck("EW at 1e-4, dM7 (A)", 19.08, round(float(ls.loc["dM7","ew_at_1e4"]),2))
ck("criterion gap dM0 (pp)", 0.6, round(float(ls.loc["dM0","diff_EW_minus_lB"]),1), 0.05)
ck("criterion gap dM6 (pp)", 36.5, round(float(ls.loc["dM6","diff_EW_minus_lB"]),1), 0.05)
mono = all(ls.loc[a,"diff_EW_minus_lB"] <= ls.loc[b,"diff_EW_minus_lB"]
           for a, b in zip(["dM0","dM1","dM2","dM3","dM4","dM5"], ["dM1","dM2","dM3","dM4","dM5","dM6"]))
ck("criterion gap monotone", 1, int(mono), 0)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step4_check.csv"), index=False, encoding="utf-8-sig")
bad = res[res.status == "FAIL"]
print(bad.to_string(index=False) if len(bad) else "  (no mismatch)")
print("  PASS %d / %d" % (int((res.status=="PASS").sum()), len(res)))
raise SystemExit(0 if len(bad) == 0 else 2)
