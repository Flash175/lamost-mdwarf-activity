# -*- coding: utf-8 -*-
"""step2_check_23.py — 把 §2.3（2D vs 6D 可用性）的数字逐项与产出表比对。"""
import os, json
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

six = json.load(open(os.path.join(OUT, "step1_6d_far.json"), encoding="utf-8"))
aud = pd.read_csv(os.path.join(D, "far400_2d_audit.csv"), **RD).set_index("zbin")
gm = pd.read_csv(os.path.join(D, "giant_cut_by_MG.csv"), **RD).set_index("zbin")
far = ["400-600", "600-1000", ">1000"]

rows = []
def ck(item, exp, got, tol=0.0, note=""):
    ok = abs(float(got)-float(exp)) <= tol if tol else int(got) == int(exp)
    rows.append(dict(item=item, expected=exp, computed=got, status="PASS" if ok else "FAIL", note=note))

# --- 6D far region（逐星表重算） ---
ck("6D far n", 316, six["n"])
ck("6D far active", 61, six["n_act"])
ck("6D far active S/N<3", 34, six["n_snr_lt3_act"])
ck("6D far active logg sentinel", 19, six["n_logg_sentinel_act"])
ck("6D far active genuine logg<3.5", 1, six["n_logg_lt35_real_act"])
ck("6D far active median G", 14.0, six["g_med_act"], 0.06)
ck("  ... of S/N<3 subset", 14.2, six["g_med_act_snr_lt3"], 0.06, "data median 14.147 -> paper 14.2 (double rounding)")
ck("  ... of sentinel subset", 13.8, six["g_med_act_sentinel"], 0.06, "data median 13.765")
ck("6D far M_G<6 (all)", 290, six["n_mg_lt6"])
ck("6D far M_G<6 (active)", 57, six["n_mg_lt6_act"])
ck("6D far inactive", 255, six["n_inact"])
ck("6D far M_G<6 (inactive)", 233, six["n_mg_lt6_inact"])
ck("6D far median M_G (bin)", 3.2, six["mg_med_bin"], 0.06)
ck("  ... of the 290 bright", 2.9, six["mg_med_bright"], 0.06)
ck("6D far median distance /kpc", 1.3, round(six["dist_med_pc"]/1000,1), 0.05)
ck("6D far after cut (n)", 26, six["n_postcut"])
ck("6D far after cut (active)", 4, six["n_postcut_act"])

# --- 2D far region（三档合并 98,692） ---
ck("2D far n (post cut)", 98692, int(aud.loc[far,"n"].sum()))
ck("2D far n_logg_valid", 96870, int(aud.loc[far,"n_logg_valid"].sum()))
ck("2D far genuine logg<3.5", 15, int(aud.loc[far,"n_logg_lt35_real"].sum()))
ck("2D far logg sentinel", 1822, int(aud.loc[far,"n_logg_missing"].sum()))
for b, exp in (("400-600",4.642),("600-1000",4.608),(">1000",4.585)):
    ck("2D logg median %s" % b, exp, round(float(aud.loc[b,"logg_med_masked"]),3), 0.001)

# --- 亮天体统计（giant_cut_by_MG 的 >400 行） ---
ck("2D far n before cut", 99684, int(gm.loc[">400","n"]))
ck("2D far bright objects", 992, int(gm.loc[">400","n_bright"]))
ck("2D far bright fraction %", 1.0, round(float(gm.loc[">400","frac_bright"]),1), 0.05)
ck("2D far active % before", 7.11, round(float(gm.loc[">400","frac_all"]),2), 0.005)
ck("2D far active % after (all 992 removed)", 7.02, round(float(gm.loc[">400","frac_clean"]),2), 0.005)
ck(">1000 pc bright objects", 369, int(992 - gm.loc["400-600","n_bright"] - gm.loc["600-1000","n_bright"]))
ck(">1000 pc n before cut", 1600, int(1231 + 369))


# --- §2.3 引用的 Gaia RV 覆盖率（与 §4.5、§4.13 同源，均由 数据/RV覆盖_by_G.csv 守护） ---
cov = pd.read_csv(os.path.join(D, "RV覆盖_by_G.csv"), **RD).set_index("g_lo")
ck("Gaia RV coverage G=12-14", 97.0, round(float(cov.loc[12, "cover_pct"]), 3), 0.5, "§2.3 两处 / §4.5 / §4.13 同值")
ck("Gaia RV coverage G=14-16", 53.0, round(float(cov.loc[14, "cover_pct"]), 3), 0.5, "同上")
ck("Gaia RV coverage G=16-18", 0.003, round(float(cov.loc[16, "cover_pct"]), 4), 0.001, "同上")

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
bad = res[res.status == "FAIL"]
print(bad.to_string(index=False) if len(bad) else "  (no mismatch)")
print("  PASS %d / %d" % (int((res.status=="PASS").sum()), len(res)))
print("  -> 输出/step2_check.csv")
