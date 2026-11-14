
# -*- coding: utf-8 -*-
"""step2_check_411.py — §4.11 的逐项比对（正文数 vs step1 重算/归档表）。"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
V = pd.read_csv(os.path.join(OUT, "step1_values.csv"), encoding="utf-8-sig").set_index("item")["value"]
try:
    B1 = pd.read_csv(os.path.join(OUT, "step1b_values.csv"), encoding="utf-8-sig").set_index("item")["value"]
except Exception:
    B1 = None
rows = []
def ck(item, body, got, tol=0.006, hard=0.05, note=""):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 4),
                     diff=round(float(got) - float(body), 4), status=st, note=note))
def ck_int(item, body, got, tol=0, note=""):
    ok = abs(int(got) - int(body)) <= tol
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got) - int(body),
                     status="PASS" if ok else "STALE", note=note))
def ck_le(item, bound, got, note=""):
    ok = float(got) <= float(bound) + 1e-9
    rows.append(dict(item=item, body=bound, computed=round(float(got), 4),
                     diff=round(float(got) - float(bound), 4), status="PASS" if ok else "STALE", note=note))

# ---- 0) 样本与基线 ----
ck_int("sample size (637,973)", 637973, V["n all"])
ck_int("dM0-dM6 subsample (637,352)", 637352, V["n dm0-dm6"])
ck("baseline decline, all subtypes (16.03 pp)", 16.03, V["decline base all"], tol=0.01)

# ---- 1) RUWE 通道 ----
ck("RUWE<1.2 removes 4.3 per cent", 4.3, V["removed pct ruwe12"], tol=0.05)
ck("RUWE<1.1 removes 13.5 per cent", 13.5, V["removed pct ruwe11"], tol=0.1)
ck("decline RUWE<1.2 (15.79 pp)", 15.79, V["decline ruwe12"], tol=0.006)
ck("decline RUWE<1.1 (15.57 pp)", 15.57, V["decline ruwe11"], tol=0.006)
ck("decline moves by 0.23 pp at RUWE<1.2", 0.234, abs(float(V["delta decline ruwe12"])), tol=0.01)
ck("decline moves by 0.46 pp at RUWE<1.1", 0.457, abs(float(V["delta decline ruwe11"])), tol=0.01)
ck("dM0-dM4 max move at RUWE<1.2 (0.28 pp)", 0.28, V["max abs ddecl ruwe12 dm0-dm4"], tol=0.006)
ck("dM0-dM4 max move at RUWE<1.1 (0.47 pp)", 0.47, V["max abs ddecl ruwe11 dm0-dm4"], tol=0.006)
ck("no subtype outside dM5/dM6 moves >0.28 pp at RUWE<1.2", 0.28, V["max abs ddecl ruwe12 dm0-dm4"], tol=0.006)
ck_le("no subtype outside dM6 moves >0.60 pp at RUWE<1.2", 0.60, V["max abs ddecl ruwe12 excl dm6"], "dM5 0.595")
ck("largest subtype move at RUWE<1.1 (dM5, 1.7 pp)", 1.7, V["max abs ddecl ruwe11 excl dm6"], tol=0.02)
ck_int("dM6 stars in the 200-400 pc bin", 114, V["dm6 far n"])
ck_int("dM6 stars below 400 pc", 3647, V["dm6 below400 n"])
ck_int("dM6 stars between 200 and 400 pc (body: 114)", 114, V["dm6 far n"])
ck_int("dM6 in that bin after RUWE<1.1 (body: 94)", 94, V["dm6 far n ruwe11"])
ck_int("dM6 in the 200-400 pc bin after RUWE<1.1", 94, V["dm6 far n ruwe11"])

# ---- 2) RV 变星通道 ----
ck_int("multi-epoch stars (131,331 in text)", 131331,
       B1["snrz5 multi-epoch"] if B1 is not None else float("nan"), note="原始星表重算；归档表 131,183 是去重后口径")
ck_int("archive multi-epoch (131,183, post-dedup)", 131183, V["rv multi-epoch archived"])
ck("multi-epoch fraction (20.6 per cent)", 20.6, 100.0 * V["rv multi-epoch archived"] / V["n all"], tol=0.1)
ck("chi2 criterion flags 53 per cent", 53.0, B1["snrz5 chi2 pct of multi"] if B1 is not None else V["flag chi2 pct"], tol=0.5,
   note="原始重算 70,068 / 131,331 = 53.35%")
ck("ptp criterion flags 26 per cent", 26.0, B1["snrz5 ptp pct of multi"] if B1 is not None else V["flag ptp pct"], tol=0.5,
   note="原始重算 33,674 / 131,331 = 25.64%")
ck_int("union listed (70,660 in text)", 70660, V["rv union rows"])
if B1 is not None:
    ck_int("union rebuilt from raw (70,660)", 70660, B1["snrz5 union"])
    ck_int("70,633 of the union fall inside the sample (in text)", 70633, B1["union snrz5 in sample"])
    ck("decline after removal, rebuilt from raw (15.93)", 15.93, B1["decline removed snrz5"], tol=0.006)
    ck("removal moves the decline by 0.096 pp (raw)", 0.096, abs(float(B1["delta decline snrz5"])), tol=0.006)
    ck("dM3 11.33 after removal (raw)", 11.33, B1["sub decl removed dM3"], tol=0.006)
    ck("dM4 10.86 after removal (raw)", 10.86, B1["sub decl removed dM4"], tol=0.006)
    ck("dM6 22.07 after removal (raw)", 22.07, B1["sub decl removed dM6"], tol=0.006)
    ck_le("no subtype moves >0.7 pp (raw)", 0.7, B1["max abs sub ddecl"])
    ck("RV variable rate among active (15.2, raw)", 15.2, B1["rate among active"], tol=0.05)
    ck("RV variable rate among inactive (10.4, raw)", 10.4, B1["rate among inactive"], tol=0.06)
    ck("200-400 pc active rate (16.9, raw)", 16.9, B1["far rate among active"], tol=0.05)
    ck("200-400 pc inactive rate (10.2, raw)", 10.2, B1["far rate among inactive"], tol=0.05)
    ck("median RV error of the multi-epoch sample (2.8 km/s)", 2.8, B1["median rv err of multi-epoch"], tol=0.3)
ck("union fraction (11.1 per cent)", 11.1, V["union pct of sample"], tol=0.1)
ck_int("stars removed by the RV cut (70,633)", 70633, V["rv removed count"])
ck("decline after RV removal (15.93 pp)", 15.93, V["table decline rvrem"], tol=0.006)
ck("RV removal moves the decline by 0.096 pp", 0.096, abs(float(V["delta decline rv removal"])), tol=0.006)
ck("dM3 11.94 before removal", 11.94, V["sub decl base dM3"], tol=0.006)
ck("dM3 11.33 after removal", 11.33, V["sub decl rvrem dM3"], tol=0.006)
ck("dM4 10.60 before removal", 10.60, V["sub decl base dM4"], tol=0.006)
ck("dM4 10.86 after removal", 10.86, V["sub decl rvrem dM4"], tol=0.006)
ck("dM6 22.18 before removal", 22.18, V["sub decl base dM6"], tol=0.006)
ck("dM6 22.07 after removal", 22.07, V["sub decl rvrem dM6"], tol=0.006)
ck_le("no subtype moves >0.7 pp under the RV cut", 0.7, V["max abs ddecl rvrem"], "最大为 dM3 0.619")
ck("RV variable rate among active (15.2 per cent)", 15.2, V["flag pct among active all"], tol=0.05)
ck("RV variable rate among inactive (10.4 per cent)", 10.4, V["flag pct among inactive all"], tol=0.06)
ck("200-400 pc: rate among active (16.9 per cent)", 16.9, V["flag pct among active far"], tol=0.05)
ck("200-400 pc: rate among inactive (10.2 per cent)", 10.2, V["flag pct among inactive far"], tol=0.05)

# ---- 3) Ro x [Fe/H] 通道 ----
ck_int("periodic stars (3,340)", 3340, V["period sample n"])
ck_int("periodic stars with usable [Fe/H] (2,984)", 2984, V["period sample with mh"])
ck("periodic [Fe/H] coverage (89.3 per cent)", 89.34, V["period sample with mh pct"], tol=0.06)
ck("median [Fe/H] uncertainty of the periodic sample (0.11 dex)", 0.11, V["period sample median mh err"], tol=0.006)
ck("dM0 dbeta (-0.08)", -0.08, V["ro dbeta dM0"], tol=0.006)
ck("dM0 dbeta error (0.25)", 0.25, V["ro dbeta se dM0"], tol=0.006)
ck_int("dM0 halves (768 stars)", 768, V["ro n dM0"])
ck("dM1 dbeta (+0.22)", 0.22, V["ro dbeta dM1"], tol=0.006)
ck("dM1 dbeta error (0.31)", 0.31, V["ro dbeta se dM1"], tol=0.006)
ck_int("dM1 halves (540 stars)", 540, V["ro n dM1"])
ck("dM2 dbeta (-0.27)", -0.27, V["ro dbeta dM2"], tol=0.006)
ck("dM2 dbeta error (0.46)", 0.46, V["ro dbeta se dM2"], tol=0.006)
ck_int("dM2 halves (515 stars)", 515, V["ro n dM2"])
ck("dM4 dbeta (+0.11)", 0.11, V["ro dbeta dM4"], tol=0.006)
ck("dM4 dbeta error (0.42)", 0.42, V["ro dbeta se dM4"], tol=0.006)
ck_int("dM4 halves (643 stars)", 643, V["ro n dM4"])
ck("precision-weighted mean dbeta (0.01)", 0.01, V["ro weighted mean dbeta"], tol=0.006)
ck("precision-weighted mean error (0.17)", 0.17, V["ro weighted mean dbeta se"], tol=0.006)
ck("two-sigma bound |dbeta| (0.35)", 0.35, V["ro weighted mean dbeta 2sigma"], tol=0.01)
ck("chi2 of the four dbeta about their weighted mean (1.0 for 3 dof)", 1.0, V["ro chi2 about mean"], tol=0.05)
ck_int("unsaturated dM3 (body: 20)", 20, V["unsaturated dM3"])
ck_int("unsaturated dM5 (body: 28)", 28, V["unsaturated dM5"])
ck_int("unsaturated dM6 (body: 3)", 3, V["unsaturated dM6"])
ck("smallest half of a fit (body: 49)", 49, V["fit half min"], tol=1)
ck("largest half of a fit (body: 209)", 209, V["fit half max"], tol=1)
ck("widest [Fe/H] half contrast (body: 0.27 dex)", 0.27, max(float(V["mh contrast dM0"]), float(V["mh contrast dM1"]),
   float(V["mh contrast dM2"]), float(V["mh contrast dM4"])), tol=0.006)
ck("narrowest [Fe/H] half contrast (body: 0.18 dex)", 0.18, min(float(V["mh contrast dM0"]), float(V["mh contrast dM1"]),
   float(V["mh contrast dM2"]), float(V["mh contrast dM4"])), tol=0.006)
ck_int("smallest half (body: 257 stars)", 257, V["ro halves min n"])
ck_int("largest half (body: 385 stars)", 385, V["ro halves max n"])

# ---- 4) 归档表与逐星重算互校 ----
for sp in ("dM0", "dM1", "dM2", "dM3", "dM4", "dM5", "dM6"):
    ck("subtype decline, archive vs recompute %s" % sp, float(V["table sub decl base %s" % sp]),
       float(V["sub decl base %s" % sp]), tol=1e-4)
    for tag in ("12", "11"):
        ck("RUWE<%s subtype decline, archive vs recompute %s" % ("1.2" if tag == "12" else "1.1", sp),
           float(V["table sub decl %s %s" % (tag, sp)]), float(V["sub decl ruwe%s %s" % (tag, sp)]), tol=1e-4)
ck("dM0-dM6 baseline, archive vs recompute", float(V["table decline base dm0-dm6"]),
   float(V["decline base dm0-dm6"]), tol=1e-4)
for lab in ("base", "ruwe12", "ruwe11", "rvrem"):
    ck("archive table decline %s (rounded to 3 dp)" % lab, round(float(V["table decline %s" % lab]), 3),
       float(V["table decline %s" % lab]), tol=1e-3)
ck("RV removal from archive vs recompute", float(V["table decline rvrem"]), float(V["decline after rv removal"]), tol=0.002)

# ---- 5) 两处口径差（都已能从原始星表重跑，记录的是"两个不同口径"而非"不可复现"）
arch = [
    dict(item="multi-epoch: raw rebuild vs archived table", value=int(B1["snrz5 multi-epoch"]) if B1 is not None else 131331,
         source="step1b_rv_chain.py vs 数据/RV变星_剔除检验.csv",
         note="原始重跑得 131,331（正文取此值）；归档表 131,183 为去重后口径，差 148 = 0.11 per cent"),
    dict(item="periodic stars with usable [Fe/H]: chosen cut vs old value", value=int(V["period sample with mh"]),
         source="数据/图5_逐星样本_严格3340行_⑤.csv + 金属度_清洁样本_637973.csv",
         note="按 §4.10 判据(|err|<=0.25, -2<=m_h<=0.6)得 2,984；正文原写 2,979，5 种邻域变体(2,814-3,145)均不得 2,979"),
    dict(item="union listed vs removed inside the sample", value=27,
         source="数据/RV变星_SB1候选_多历元.csv / step1b",
         note="清单与原始重算并集均为 70,660，其中 27 个 uid 不在 637,973 星样本内，故实际剔除 70,633"),
]
pd.DataFrame(arch).to_csv(os.path.join(OUT, "step2_archive_only.csv"), index=False, encoding="utf-8-sig")
print("archive-only rows:", len(arch))

R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
print(R.status.value_counts().to_dict(), "| total", len(R))
for _, r in R[R.status != "PASS"].iterrows():
    print("  %-6s %-56s body=%-10s got=%-10s %s" % (r.status, r["item"], r.body, r.computed, r.note))
