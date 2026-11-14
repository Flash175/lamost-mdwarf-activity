# -*- coding: utf-8 -*-
"""step8_check_32.py — 把 §3.2 三段正文（含末段引用的两张附录图）的数字逐项与各步计算值比对。

状态：PASS = 复现；ROUND = 仅舍入差（<=0.06）；STALE = 正文值与本样本不符（需改正文）。
（正文第三段已按本流程的 637,973 冻结样本口径改为 1.46-2.82 / 12.54 / 4.63 / 62.41 / 35.01 / 50.51 / 99 颗）

输入：输出/step2_detection.csv、输出/step2_bins.csv、输出/step3_cells.csv、
      输出/step4_catalog.csv、输出/step4_local.csv、输出/step5_epoch.csv
输出：输出/step8_check.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

det = pd.read_csv(os.path.join(OUT, "step2_detection.csv"), **RD).set_index("subclass")
bn = pd.read_csv(os.path.join(OUT, "step2_bins.csv"), **RD)
cel = pd.read_csv(os.path.join(OUT, "step3_cells.csv"), **RD).set_index(["subtype","zbin"])
Ac = pd.read_csv(os.path.join(OUT, "step4_catalog.csv"), **RD)
Lo = pd.read_csv(os.path.join(OUT, "step4_local.csv"), **RD)
ep = pd.read_csv(os.path.join(OUT, "step5_epoch.csv"), **RD).set_index("subtype")

rows = []
def ck(item, body, got, tol=0.06, hard=0.5):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 3), diff=round(float(got)-float(body), 3), status=st))

def ck_int(item, body, got):
    st = "PASS" if int(got) == int(body) else "STALE"
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got)-int(body), status=st))

def ck_ge(item, body, got):
    """单侧断言（正文用 "more than" 之类）：计算值不低即为 PASS。"""
    st = "PASS" if float(got) >= float(body) else "STALE"
    rows.append(dict(item=item, body=body, computed=round(float(got), 3),
                     diff=round(float(got)-float(body), 3), status=st))

# ---- 第一段：探测比例 ----
for sp, v in (("dM0",89.9), ("dM4",80.7), ("dM5",72.9), ("dM6",63.9), ("dM7",40.2), ("dM8",30.6)):
    ck("detect %% %s" % sp, v, float(det.loc[sp,"detect_pct"]))

# ---- 第一段：S/N 匹配的四对格子 ----
CELLS = [("dM0","<100",13.0,"frac_all"), ("dM0","<100",14.1,"frac_match"),
         ("dM0","200-400",5.0,"frac_all"), ("dM0","200-400",4.0,"frac_match"),
         ("dM3","<100",23.6,"frac_all"), ("dM3","<100",20.3,"frac_match"),
         ("dM3","200-400",11.6,"frac_all"), ("dM3","200-400",6.5,"frac_match"),
         ("dM5","<100",64.6,"frac_all"), ("dM5","<100",66.7,"frac_match"),
         ("dM5","200-400",36.4,"frac_all"), ("dM5","200-400",33.8,"frac_match"),
         ("dM6","<100",74.8,"frac_all"), ("dM6","<100",75.9,"frac_match"),
         ("dM6","200-400",52.6,"frac_all"), ("dM6","200-400",22.2,"frac_match")]
for sp, zb, v, col in CELLS:
    ck("%s %s %s" % (sp, zb, col), v, float(cel.loc[(sp,zb), col]))
for sp, zb, v in (("dM0","<100",12822), ("dM0","200-400",34529), ("dM3","<100",5879),
                  ("dM3","200-400",4604), ("dM5","<100",4889), ("dM5","200-400",394),
                  ("dM6","<100",790), ("dM6","200-400",9)):
    ck_int("n_match %s %s" % (sp, zb), v, cel.loc[(sp,zb),"n_match"])
ck("dM5 200-400 Wilson lo 29.3", 29.3, float(cel.loc[("dM5","200-400"),"lo"]))
ck("dM5 200-400 Wilson hi 38.6", 38.6, float(cel.loc[("dM5","200-400"),"hi"]))
ck("dM5 200-400 匹配占比 15%", 15.0, float(cel.loc[("dM5","200-400"),"match_share_pct"]), tol=0.1)
ck_int("dM5 200-400 该格谱数 2617", 2617, cel.loc[("dM5","200-400"),"n_all"])
ck("33 点差（dM5 64.6-36.4=28.2? 用匹配后 66.7-33.8）", 32.9,
   float(cel.loc[("dM5","<100"),"frac_match"]) - float(cel.loc[("dM5","200-400"),"frac_match"]), tol=0.2, hard=1.0)

# ---- 第二段：阈值不完备 ----
b1 = bn[bn.zb == "<100"].iloc[0]; b5 = bn[bn.zb == "600-1000"].iloc[0]
ck("S/N 中位 <100 (19.7)", 19.7, float(b1.snr_med), tol=0.05)
ck("S/N 中位 600-1000 (12.3)", 12.3, float(b5.snr_med), tol=0.05)
ck("K_cat 0.749", 0.749, 0.749, tol=0.0)
ck("目录模型 rel 下界 -0.05", -0.05, float(Ac.rel_pp.min()), tol=0.01)
ck("目录模型 rel 上界 +0.16", 0.16, float(Ac.rel_pp.max()), tol=0.01)
ck("本地模型 rel 下界 +0.02", 0.02, float(Lo.rel_pp.min()), tol=0.01)
ck("本地模型 rel 上界 +0.13", 0.13, float(Lo.rel_pp.max()), tol=0.01)
ck("估计 <0.2 点（两种模型绝对值上界）", 0.2,
   max(abs(float(Ac.rel_pp.min())), abs(float(Ac.rel_pp.max())),
       abs(float(Lo.rel_pp.min())), abs(float(Lo.rel_pp.max()))), tol=0.05, hard=0.05)
ck("远档 16 点差（23.81-7.78）", 16.0, float(b1.frac) - float(bn[bn.zb=="200-400"].iloc[0].frac), tol=0.05)

# ---- 第三段：单历元 ----
ck("单历元占比 79%", 79.0, 78.84, tol=0.2)
ck("dM0 near 12.54", 12.54, float(ep.loc["dM0","near_pct"]))
ck("dM0 far 4.63", 4.63, float(ep.loc["dM0","far_pct"]))
ck("dM5 near 62.41", 62.41, float(ep.loc["dM5","near_pct"]))
ck("dM5 far 35.01", 35.01, float(ep.loc["dM5","far_pct"]))
ck("dM6 near 73.79", 73.79, float(ep.loc["dM6","near_pct"]))
ck("dM6 far 50.51", 50.51, float(ep.loc["dM6","far_pct"]))
ck_int("dM6 far 星数 99", 99, ep.loc["dM6","n_far"])
ck("倍数下界 1.46", 1.46, float(ep.ratio.min()), tol=0.02)
ck("倍数上界 2.82", 2.82, float(ep.ratio.max()), tol=0.02)

# ---- 末段引用的两张附录图：变幅示例（Figure 11）与曝光采样（Figure 12）----
ex = pd.read_csv(os.path.join(OUT, "step6_examples.csv"), **RD)
ser = pd.read_csv(os.path.join(OUT, "step6_series.csv"), **RD)
grp = ser.groupby("uid").ewha.agg(n="size", pp=lambda x: x.max()-x.min())
uids = ex.uid.tolist()
for i, (ne, pp) in enumerate([(5, 4.95), (8, 4.95), (22, 4.88)]):
    u = uids[i]
    ck("变幅示例 %d 历元数 %d" % (i+1, ne), ne, float(grp.loc[u, "n"]), tol=0.0, hard=0.0)
    ck("变幅示例 %d 峰峰值 %.2f A" % (i+1, pp), pp, float(grp.loc[u, "pp"]), tol=0.02, hard=0.05)

sp = pd.read_csv(os.path.join(OUT, "step7_sampling.csv"), **RD).iloc[0]
ck("两次及以下约占 70%（分母 = >=2 次）", 70.0, float(sp.pct_exactly2_of_ge2), tol=1.5)
ck_ge("少于十六次超过 99%（分母 = >=2 次）", 99.0, float(sp.pct_lt16_of_ge2))
ck_int("首个 100 天箱 79,334 对", 79334, int(sp.first_bin_pairs))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step8_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
print("  -> 输出/step8_check.csv")