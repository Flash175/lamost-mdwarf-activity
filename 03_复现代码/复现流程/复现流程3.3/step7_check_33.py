# -*- coding: utf-8 -*-
"""step7_check_33.py — 把 §3.3（含附录 app:rv）的数字逐项与计算值比对。

输入：输出/step1_6d.csv、step2_bins.csv、step2_slices.csv、step3_errors.csv、
      step3_window.csv、step3_strata.csv、step4_lamost.csv、step4_counts.csv、
      step5_crosscheck.csv、step5_summary.csv
输出：输出/step7_check.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

bins = pd.read_csv(os.path.join(OUT, "step2_bins.csv"), **RD).set_index("zbin")
sl = pd.read_csv(os.path.join(OUT, "step2_slices.csv"), **RD).set_index("zslice")
er = pd.read_csv(os.path.join(OUT, "step3_errors.csv"), **RD).set_index("slice")
win = pd.read_csv(os.path.join(OUT, "step3_window.csv"), **RD).iloc[0]
gr = pd.read_csv(os.path.join(OUT, "step3_strata.csv"), **RD)
lam = pd.read_csv(os.path.join(OUT, "step4_lamost.csv"), **RD).set_index("z_lo")
lc = pd.read_csv(os.path.join(OUT, "step4_counts.csv"), **RD).iloc[0]
cc = pd.read_csv(os.path.join(OUT, "step5_crosscheck.csv"), **RD).set_index("z_lo")
cs = pd.read_csv(os.path.join(OUT, "step5_summary.csv"), **RD).iloc[0]
n6d = len(pd.read_csv(os.path.join(OUT, "step1_6d.csv"), **RD))
t5 = pd.read_csv(os.path.join(OUT, "step6_table5.csv"), **RD).set_index("zbin")
t5s = pd.read_csv(os.path.join(OUT, "step6_summary.csv"), **RD).iloc[0]

rows = []
def ck(item, body, got, tol=0.06):
    d = abs(float(got) - float(body))
    rows.append(dict(item=item, body=body, computed=round(float(got), 3),
                     diff=round(float(got)-float(body), 3),
                     status="PASS" if d <= tol else "FAIL"))
def ck_int(item, body, got):
    rows.append(dict(item=item, body=int(body), computed=int(got),
                     diff=int(got)-int(body),
                     status="PASS" if int(got) == int(body) else "FAIL"))

# ---------- 第一段：6D 子样本 ----------
ck_int("6D subsample", 82518, n6d)
SIG = {"0-50": (14.4, 19.6), "50-100": (18.3, 20.3), "100-150": (21.5, 21.2), "150-200": (26.2, 21.8)}
for k, (va, vi) in SIG.items():
    ck("sigma_vz act %s" % k, va, float(bins.loc[k, "sd_act"]))
    ck("sigma_vz inact %s" % k, vi, float(bins.loc[k, "sd_inact"]))
ZMED = {"0-50": (22.1, 25.6), "50-100": (71.6, 74.5), "100-150": (118.7, 122.4), "150-200": (169.8, 170.8)}
for k, (va, vi) in ZMED.items():
    ck("median |Z| act %s" % k, va, float(bins.loc[k, "zmed_act"]))
    ck("median |Z| inact %s" % k, vi, float(bins.loc[k, "zmed_inact"]))

DSTD = {"0-25": -6.1, "25-50": -4.0, "50-75": -4.4, "75-100": 0.8, "100-150": 0.3, "150-200": 4.4}
DMAD = {"0-25": -7.8, "25-50": -6.4, "50-75": -6.5, "75-100": -5.4, "100-150": -4.4, "150-200": 2.5}
for k, v in DSTD.items():
    ck("d_sigma_std %s" % k, v, float(sl.loc[k, "d_std"]), tol=0.06)
for k, v in DMAD.items():
    ck("d_sigma_MAD %s" % k, v, float(sl.loc[k, "d_mad"]), tol=0.06)
w = (sl.n_act + sl.n_inact).astype(float)
ck("count-weighted mean", -1.56, float((sl.d_std*w).sum()/w.sum()), tol=0.006)

# ---------- 第一段/§4：25 pc 片的 |Z| 中位（两类各自的中位，非与档边界比）----------
for k, (za, zi, dd) in {"0-25": (12.1, 12.4, 0.3), "25-50": (36.6, 37.7, 1.1), "50-75": (61.7, 62.5, 0.9)}.items():
    ck("zmed act %s" % k, za, float(sl.loc[k, "zmed_act"]), tol=0.06)
    ck("zmed inact %s" % k, zi, float(sl.loc[k, "zmed_inact"]), tol=0.06)
    ck("zmed gap %s" % k, dd, float(sl.loc[k, "zmed_inact"] - sl.loc[k, "zmed_act"]), tol=0.06)

# ---------- 第一段：40-80 pc 窗口 ----------
ck_int("window n_act", 2957, win.n_act); ck_int("window n_inact", 16399, win.n_inact)
ck("window |Z| med act", 59.0, float(win.zmed_act), tol=0.06)
ck("window |Z| med inact", 60.3, float(win.zmed_inact), tol=0.06)
ck("window sigma act", 15.8, float(win.sd_act)); ck("window sigma inact", 19.7, float(win.sd_inact))
ck("window MAD act", 10.6, float(win.mad_act)); ck("window MAD inact", 16.6, float(win.mad_inact))

# ---------- 第一段：标准误与显著性 ----------
for k, (d, z) in {"0-25": (-6.08, -26.0), "25-50": (-4.03, -14.2), "50-75": (-4.39, -15.3)}.items():
    ck("d_std %s" % k, d, float(er.loc[k, "d_std"]), tol=0.006)
    ck("z_std %s" % k, z, float(er.loc[k, "z_std"]), tol=0.06)
for k, (d, z) in {"0-25": (-7.82, -12.7), "25-50": (-6.43, -8.6), "50-75": (-6.49, -8.6)}.items():
    ck("d_MAD %s" % k, d, float(er.loc[k, "d_mad"]), tol=0.006)
    ck("z_MAD %s" % k, z, float(er.loc[k, "z_mad"]), tol=0.06)
for k, (zs, zm) in {"75-100": (1.9, -4.7), "100-150": (0.7, -4.4), "150-200": (5.9, 1.3)}.items():
    ck("z_std %s" % k, zs, float(er.loc[k, "z_std"]), tol=0.06)
    ck("z_MAD %s" % k, zm, float(er.loc[k, "z_mad"]), tol=0.06)
ck("min |z_MAD| in 0-75 pc", 8.6, float(min(abs(er.loc[k, "z_mad"]) for k in ("0-25","25-50","50-75"))), tol=0.06)
ck("window d_std", -3.91, float(win.d_std), tol=0.006)
ck("window z_std", -16.8, float(win.z_std), tol=0.06)
ck("window d_MAD", -6.02, float(win.d_mad), tol=0.006)
ck("window z_MAD", -9.8, float(win.z_mad), tol=0.06)

# ---------- 第三段：亚型分层 ----------
EARLY = {"0-25": -4.96, "25-50": -3.77, "50-75": -3.25, "40-80": -3.32}
LATE = {"0-25": -6.84, "25-50": -3.97, "50-75": -6.46, "40-80": -5.14}
for k, v in EARLY.items():
    ck("early d_std %s" % k, v, float(gr[(gr.group=="dM0-dM2") & (gr.slice==k)].d_std.iloc[0]), tol=0.006)
for k, v in LATE.items():
    ck("late d_std %s" % k, v, float(gr[(gr.group=="dM3-dM5") & (gr.slice==k)].d_std.iloc[0]), tol=0.006)
ze = gr[gr.group=="dM0-dM2"].z_std.abs()
zl = gr[gr.group=="dM3-dM5"].z_std.abs()
ck("early z min", 7.7, float(ze.min()), tol=0.06); ck("early z max", 14.2, float(ze.max()), tol=0.06)
ck("late z min", 9.1, float(zl.min()), tol=0.06);  ck("late z max", 19.7, float(zl.max()), tol=0.06)

# ---------- 第二段：LAMOST-RV ----------
ck_int("LAMOST-RV stars", 611094, lc.n_stars)
ck("median sigma_rv", 3.1, float(lc.median_sig_rv), tol=0.05)
ck_int("|Z| 400-600 stars", 68622, lc.n_400_600)
ck_int("|Z| 600-1000 stars", 24495, lc.n_600_1000)
ck_int("both velocities", 81719, cs.n_both)
ck("dRV>50 km/s per cent", 2.2, float(cs.pct_gt50), tol=0.06)
for lo, v in zip([0,25,50,75,100], [-5.20,-5.06,-4.68,-4.02,-3.64]):
    ck("LAMOST d_std %d" % lo, v, float(lam.loc[lo, "d_std"]), tol=0.006)
for lo, v in zip([0,25,50,75,100], [-6.89,-6.69,-6.72,-7.28,-6.43]):
    ck("LAMOST d_MAD %d" % lo, v, float(lam.loc[lo, "d_mad"]), tol=0.006)
ck("LAMOST z range lo", -37.0, float(lam.loc[[0,25,50,75,100], "z_std"].min()), tol=0.6)
ck("LAMOST z range hi", -21.0, float(lam.loc[[0,25,50,75,100], "z_std"].max()), tol=0.6)
ck("LAMOST z_std 125-150", -0.5, float(lam.loc[125, "z_std"]), tol=0.06)
ck("LAMOST z_std 150-175", 4.9, float(lam.loc[150, "z_std"]), tol=0.06)
ck("LAMOST z_std 175-200", 9.3, float(lam.loc[175, "z_std"]), tol=0.06)
sub = lam[lam.index >= 200]
ck("beyond 200 d_std lo", 4.27, float(sub.d_std.min()), tol=0.006)
ck("beyond 200 d_std hi", 10.63, float(sub.d_std.max()), tol=0.006)
ck("beyond 200 z lo", 3.4, float(sub.z_std.min()), tol=0.06)
ck("beyond 200 z hi", 10.7, float(sub.z_std.max()), tol=0.06)

# ---------- 第二段：交叉核对 ----------
ck("X d_std Gaia 0-50", -5.88, float(cc.loc[0, "d_std_gaia"]), tol=0.006)
ck("X d_std LAMOST 0-50", -5.70, float(cc.loc[0, "d_std_lam"]), tol=0.006)
ck("X d_std Gaia 25-50", -4.57, float(cc.loc[25, "d_std_gaia"]), tol=0.006)
ck("X d_std LAMOST 25-50", -5.13, float(cc.loc[25, "d_std_lam"]), tol=0.006)
ck("X d_std Gaia 175-200", 13.46, float(cc.loc[175, "d_std_gaia"]), tol=0.006)
ck("X d_std LAMOST 175-200", 6.16, float(cc.loc[175, "d_std_lam"]), tol=0.006)

# ---------- Table 5（附录 tab:kin6d）----------
T5 = {
 "0-50":    (4609, 31.3, 20.0, 14.4, 37.3, 25.1, 19.6, 22.1, 25.6),
 "50-100":  (3288, 31.4, 21.0, 18.3, 36.4, 24.7, 20.3, 71.6, 74.5),
 "100-150": (1722, 33.0, 24.7, 21.5, 37.7, 26.1, 21.2, 118.7, 122.4),
 "150-200": ( 641, 38.0, 26.2, 26.2, 39.3, 26.9, 21.8, 169.8, 170.8),
 "200-400": ( 473, 41.2, 29.5, 39.0, 40.0, 27.5, 22.5, 242.8, 232.6),
 "400-3000":(   8, 38.0, 16.9, 47.2, 39.8, 28.4, 22.4, 411.8, 421.8),
}
for k, v in T5.items():
    n_act, ax, ay, az, ix, iy, iz, za, zi = v
    ck_int("Table5 n_act %s" % k, n_act, t5.loc[k, "n_act"])
    for lab, val in (("svx_act", ax), ("svy_act", ay), ("svz_act", az),
                     ("svx_inact", ix), ("svy_inact", iy), ("svz_inact", iz)):
        ck("Table5 %s %s" % (lab, k), val, float(t5.loc[k, lab]))
    ck("Table5 zmed_act %s" % k, za, float(t5.loc[k, "zmed_act"]))
    ck("Table5 zmed_inact %s" % k, zi, float(t5.loc[k, "zmed_inact"]))
ck("Table5 far fraction (per cent)", 0.06, float(t5s.pct_beyond400), tol=0.006)
ck_int("Table5 n beyond 400 pc", 47, t5s.n_beyond400)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step7_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step6_check.csv")