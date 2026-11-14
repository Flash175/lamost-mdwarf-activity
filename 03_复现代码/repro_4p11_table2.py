
# -*- coding: utf-8 -*-
"""step3_table2.py — §4.11 的 Table 2（tab:channels）逐行核对 + 两张交付表。

Table 2 有 5 行通道，其中 3 行带可核对的数：RUWE<1.1、11 per cent / 0.1 pp、|Delta beta| <~ 0.35。
本步把这三行逐格核对，并写出逐亚型位移表与 Ro x [Fe/H] 试探表。
输出：输出/step3_table2_check.csv、输出/step3_summary.csv、
      数据/表4.11_亚型位移.csv、数据/表4.11_Ro金属度试探.csv
"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出")

V = pd.read_csv(os.path.join(OUT, "step1_values.csv"), encoding="utf-8-sig").set_index("item")["value"]
B1 = pd.read_csv(os.path.join(OUT, "step1b_values.csv"), encoding="utf-8-sig").set_index("item")["value"]
SUB = ["dM0", "dM1", "dM2", "dM3", "dM4", "dM5", "dM6"]

# ---------- 1) Table 2 逐行 ----------
rows = []
def ck(channel, claim, body, got, tol, note=""):
    d = abs(float(got) - float(body))
    rows.append(dict(channel=channel, claim=claim, printed=body, computed=round(float(got), 4),
                     diff=round(float(got) - float(body), 4), tol=tol,
                     status="PASS" if d <= tol else "STALE", note=note))

ck("astrometric binaries (RUWE proxy)", "RUWE<1.1 leaves the decline unchanged: relative change (per cent)",
   2.9, 100.0 * abs(float(V["delta decline ruwe11"])) / float(V["decline base all"]), 0.2,
   "绝对变化 0.46 pp，相对 2.9 per cent")
ck("astrometric binaries (RUWE proxy)", "no subtype moves more than 0.3 pp (dM0-dM4)",
   0.3, float(V["max abs ddecl ruwe12 dm0-dm4"]), 0.02)
ck("RV-variable stars", "removing 11 per cent of the sample: fraction removed",
   11.0, float(B1["removed pct snrz5"]), 0.1)
ck("RV-variable stars", "moves the decline by 0.1 pp",
   0.1, abs(float(B1["delta decline snrz5"])), 0.006)
ck("metallicity via tau_c (Ro at fixed P)", "|Delta beta| <~ 0.35 between [Fe/H] halves",
   0.35, float(V["ro weighted mean dbeta 2sigma"]), 0.01, "0.0104 + 2 x 0.1659")
ck("metallicity via tau_c (Ro at fixed P)", "precision-weighted mean is consistent with zero",
   0.0, float(V["ro weighted mean dbeta"]), 0.02, "chi2 = 1.0 for 3 dof")
rows.append(dict(channel="age (not measured here)", claim="status: not tested", printed="no age column",
                 computed="no age column", diff=0.0, tol=0.0, status="PASS",
                 note="逐星表与周期表都没有年龄列，与本行 'not tested' 一致"))
rows.append(dict(channel="residual selection or detection", claim="status: bounded (threshold band)",
                 printed="Section 4.3", computed="复现流程4.3", diff=0.0, tol=0.0, status="PASS",
                 note="由 复现流程4.3 的阈值带检验覆盖，本节不重复计算"))
T2 = pd.DataFrame(rows)
T2.to_csv(os.path.join(OUT, "step3_table2_check.csv"), index=False, encoding="utf-8-sig")

# ---------- 2) 逐亚型位移表 ----------
sub = []
for sp in SUB:
    sub.append(dict(subtype=sp,
                    decline_base=round(float(V["sub decl base %s" % sp]), 3),
                    decline_ruwe12=round(float(V["sub decl ruwe12 %s" % sp]), 3),
                    ddecl_ruwe12=round(float(V["sub ddecl ruwe12 %s" % sp]), 3),
                    decline_ruwe11=round(float(V["sub decl ruwe11 %s" % sp]), 3),
                    ddecl_ruwe11=round(float(V["sub ddecl ruwe11 %s" % sp]), 3),
                    decline_rvremoved=round(float(B1["sub decl removed %s" % sp]), 3),
                    ddecl_rvremoved=round(float(B1["sub ddecl %s" % sp]), 3)))
sub.append(dict(subtype="all (dM0-dM9)", decline_base=round(float(V["decline base all"]), 3),
                decline_ruwe12=round(float(V["decline ruwe12"]), 3),
                ddecl_ruwe12=round(float(V["delta decline ruwe12"]), 3),
                decline_ruwe11=round(float(V["decline ruwe11"]), 3),
                ddecl_ruwe11=round(float(V["delta decline ruwe11"]), 3),
                decline_rvremoved=round(float(B1["decline removed snrz5"]), 3),
                ddecl_rvremoved=round(float(B1["delta decline snrz5"]), 3)))
S = pd.DataFrame(sub)
S.to_csv(os.path.join(D, "表4.11_亚型位移.csv"), index=False, encoding="utf-8-sig")

# ---------- 3) Ro x [Fe/H] 试探表 ----------
Rt = pd.read_csv(os.path.join(D, "金属度_Ro试点.csv"), encoding="utf-8-sig").set_index("subtype")
rr = []
for sp in ["dM0", "dM1", "dM2", "dM4"]:
    r = Rt.loc[sp]
    rr.append(dict(subtype=sp, n=int(r.n), n_low=int(r.n_low), n_high=int(r.n_high),
                   mh_split=round(float(r.mh_split), 3), mh_low=round(float(r.mh_low), 3),
                   mh_high=round(float(r.mh_high), 3), contrast=round(float(r.mh_high - r.mh_low), 3),
                   dbeta=round(float(r.dbeta), 3), dbeta_se=round(float(r.dbeta_se), 3),
                   sigma_int=round(float(r.sigma), 3),
                   n_fit_low=int(V["fit half lo %s" % sp]), n_fit_high=int(V["fit half hi %s" % sp])))
for sp in ["dM3", "dM5", "dM6"]:
    rr.append(dict(subtype=sp, n=int(V["ro n %s" % sp]) if ("ro n %s" % sp) in V.index else 0,
                   n_low=np.nan, n_high=np.nan, mh_split=np.nan, mh_low=np.nan, mh_high=np.nan,
                   contrast=np.nan, dbeta=np.nan, dbeta_se=np.nan, sigma_int=np.nan,
                   n_fit_low=int(V["unsaturated %s" % sp]), n_fit_high=0))
R = pd.DataFrame(rr)
R.to_csv(os.path.join(D, "表4.11_Ro金属度试探.csv"), index=False, encoding="utf-8-sig")

# ---------- 4) 小结 ----------
summ = pd.DataFrame([
    dict(item="Table 2 rows checked", value=len(T2)),
    dict(item="Table 2 PASS", value=int((T2.status == "PASS").sum())),
    dict(item="subtype decline rows written", value=len(S)),
    dict(item="Ro test rows written", value=len(R)),
])
summ.to_csv(os.path.join(OUT, "step3_summary.csv"), index=False, encoding="utf-8-sig")

print("Table 2:", len(T2), T2.status.value_counts().to_dict())
print("subtype table:", S.shape, "| Ro table:", R.shape)
for _, r in T2[T2.status != "PASS"].iterrows():
    print("  %-6s %-32s printed=%-8s computed=%-8s %s" % (r.status, r.channel, r.printed, r.computed, r.note))
