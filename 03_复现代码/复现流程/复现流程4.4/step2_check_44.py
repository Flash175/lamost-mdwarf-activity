# -*- coding: utf-8 -*-
"""step2_check_44.py — §4.4 的逐项比对（26 项）。"""
import os
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "输出")
V = pd.read_csv(os.path.join(OUT, "step1_values.csv"), encoding="utf-8-sig").set_index("item")["value"]
rows = []
def ck(item, body, got, tol=0.006, hard=0.05):
    d = abs(float(got) - float(body))
    st = "PASS" if d <= tol else ("ROUND" if d <= hard else "STALE")
    rows.append(dict(item=item, body=body, computed=round(float(got), 4),
                     diff=round(float(got)-float(body), 4), status=st))
def ck_int(item, body, got, tol=0):
    ok = abs(int(got) - int(body)) <= tol
    rows.append(dict(item=item, body=int(body), computed=int(got), diff=int(got)-int(body),
                     status="PASS" if ok else "STALE"))

# 内盘指数
ck("exp a (-1.2262)", -1.2262, float(V["exp a"]), tol=0.0002)
ck("exp b (-0.004920)", -0.004920, float(V["exp b"]), tol=0.000002)
ck("exp H (203 pc)", 203.0, float(V["exp H"]), tol=0.5)
ck("exp A (0.2934)", 0.2934, float(V["exp A"]), tol=0.0002)
ck("extrap 400-600 (2.84)", 2.84, float(V["extrap 400-600"]), tol=0.01)
ck("extrap 600-1000 (1.00)", 1.00, float(V["extrap 600-1000"]), tol=0.01)
ck("ratio 400-600 (2.35)", 2.35, float(V["ratio 400-600"]), tol=0.02)
ck("ratio 600-1000 (7.79)", 7.79, float(V["ratio 600-1000"]), tol=0.02)

# 二项加权
ck("binom a (-1.1783)", -1.1783, float(V["binom a"]), tol=0.0002)
ck("binom b (-0.005198)", -0.005198, float(V["binom b"]), tol=0.000002)
ck("binom H (192)", 192.0, float(V["binom H"]), tol=0.5)
ck("binom H_se (1.5)", 1.45, float(V["binom H_se"]), tol=0.1)

# S/N>=15
ck("S/N>=15 400-600 (3.84)", 3.84, float(V["hi 400-600"]), tol=0.01)
ck("S/N>=15 600-1000 (5.45)", 5.45, float(V["hi 600-1000"]), tol=0.01)

# 分类检验
for sp, v in (("dM0",2.43),("dM1",2.75),("dM2",2.27),("dM3",1.83),("dM4",1.58),("dM5",1.13)):
    ck("C1 %s" % sp, v, float(V["C1 %s" % sp]), tol=0.006)
ck_int("C1 dM5 far (672)", 672, V["C1 dM5 nfar"])
for sp, v in (("dM0",2.55),("dM1",2.32),("dM2",1.94),("dM3",1.67),("dM4",1.39),("dM5",1.08)):
    ck("C2 %s" % sp, v, float(V["C2 %s" % sp]), tol=0.006)

# 星等限制
ck_int("G<=16.5 n 400-600 (5,834)", 5834, V["G16.5 n 400-600"])
ck("G<=16.5 pct 400-600 (8.59)", 8.59, float(V["G16.5 pct 400-600"]), tol=0.006)
ck_int("G<=16.5 n 600-1000 (350)", 350, V["G16.5 n 600-1000"])
ck("G<=16.5 pct 600-1000 (18.00)", 18.00, float(V["G16.5 pct 600-1000"]), tol=0.006)
ck_int("M_G 6-7 n (1,070)", 1070, V["MG6-7 n"])
ck("M_G 6-7 share (4.2%)", 4.2, float(V["MG6-7 share"]), tol=0.05)
ck("M_G 6-7 pct (16.82)", 16.82, float(V["MG6-7 pct"]), tol=0.02)
ck_int("M_G>=7 n (24,703)", 24703, V["MGge7 n"])
ck("M_G>=7 pct (7.37)", 7.37, float(V["MGge7 pct"]), tol=0.02)
sh = [float(V["MGge7 shift %s" % k]) for k in ("in100","m200","m400","m600")]
ck("M_G>=7 shift most negative (-0.39)", -0.39, min(sh), tol=0.02, hard=0.05)
ck("M_G>=7 shift least negative (-0.05)", -0.05, max(sh), tol=0.02, hard=0.05)
ck("M_G>=7 shift >1000 (-3.13)", -3.13, float(V["MGge7 shift far"]), tol=0.02)
ck("floor 400-600 after M_G>=7 (6.54)", 6.54, 6.651 + float(V["MGge7 shift m400"]), tol=0.02)

# Wilson 区间（外两档）
for zb, lo, hi in (('400-600', 6.47, 6.84), ('600-1000', 7.44, 8.10)):
    ck('wilson lo %s' % zb, lo, float(V['wilson lo %s' % zb]), tol=0.006)
    ck('wilson hi %s' % zb, hi, float(V['wilson hi %s' % zb]), tol=0.006)

# 逐亚型外推倍数
ck('subtype ratio min 400-600 (1.54)', 1.54, float(V['subtype ratio min 400-600']), tol=0.006)
ck('subtype ratio max 400-600 (2.81)', 2.81, float(V['subtype ratio max 400-600']), tol=0.006)
ck('subtype ratio min 600-1000 (2.67)', 2.67, float(V['subtype ratio min 600-1000']), tol=0.006)
ck('subtype ratio max 600-1000 (9.67)', 9.67, float(V['subtype ratio max 600-1000']), tol=0.006)
ck_int('subtype count per bin (6)', 6, V['subtype ratio n 400-600'])

# 七结点分段线性 vs 单指数
ck('7-knot RMS (0.43)', 0.43, float(V['L2 rms pp']), tol=0.006)
ck('7-knot worst bin (1.01)', 1.01, float(V['L2 max pp']), tol=0.006)
ck('7-knot n_par (8)', 8.0, float(V['L2 n_par']))
ck('exponential RMS (5.96)', 5.96, float(V['L1 rms pp']), tol=0.006)
ck('exponential worst bin (12.6)', 12.6, float(V['L1 max pp']), tol=0.05)

# 判别力
ck('AUC lower (0.67)', 0.67, float(V['AUC min']), tol=0.006)
ck('AUC upper (0.72)', 0.72, float(V['AUC max']), tol=0.006)
ck('Tjur R2 upper (0.12)', 0.12, float(V['tjur max']), tol=0.006)

# 内盘尺度高
ck('inner per-star H (202)', 202.0, float(V['inner H exp']), tol=0.45)
ck('inner per-star H se (2)', 2.0, float(V['inner H exp se']), tol=0.45)
ck('H rel diff pct (0.8)', 0.8, float(V['H rel diff pct']), tol=0.06)
ck('inner subtype-controlled H (287)', 287.0, float(V['inner H subtypeFE']), tol=0.5)
ck('inner subtype-controlled H se (3.5)', 3.5, float(V['inner H subtypeFE se']), tol=0.05)

# 外推口径
ck('starlevel extrap 400-600 (2.78)', 2.78, float(V['starlevel pred 400-600']), tol=0.006)
ck('starlevel extrap 600-1000 (0.94)', 0.94, float(V['starlevel pred 600-1000']), tol=0.006)
ck('median extrap 400-600 (2.79)', 2.79, float(V['median pred 400-600']), tol=0.006)
ck('median extrap 600-1000 (0.97)', 0.97, float(V['median pred 600-1000']), tol=0.006)
ck('holdout A 400-600 (6.12)', 6.12, float(V['holdoutA pred 400-600']), tol=0.006)
ck('holdout A 600-1000 (4.69)', 4.69, float(V['holdoutA pred 600-1000']), tol=0.006)

# 星等限：maglimit 子样本（正文数字的来源）

# 星等限：AIC 比较
ck('AIC single G (-290)', -290.0, float(V['AIC dG']), tol=0.5)
ck('AIC single M_G (-925)', -925.0, float(V['AIC dMG']), tol=0.5)
ck('AIC single ln d (-1,395)', -1395.0, float(V['AIC dlogd']), tol=0.5)
ck('AIC pair (-1,599)', -1599.0, float(V['AIC pair dMG']), tol=0.5)
ck('AIC pair identical (spread 0)', 0.0, float(V['AIC pair spread']), tol=0.01)

# 地板
ck('floor 400-600 after M_G>=7 (6.54)', 6.54, float(V['MGge7 pct m400']), tol=0.006)

# Hβ 3 万谱分层
ck_int('Hbeta stratified total (30,224)', 30224, V['hb30k total'])
for zb, n in (('<100', 8236), ('100-200', 7935), ('200-400', 9036), ('400-600', 3452), ('600-1000', 1265), ('>1000', 300)):
    ck_int('Hbeta n %s (%d)' % (zb, n), n, V['hb30k n %s' % zb])
for zb, c, a, b in (('<100', 30.09, 29.15, 20.07), ('100-200', 15.70, 14.32, 9.39),
                    ('200-400', 8.36, 6.96, 4.73), ('400-600', 6.84, 5.74, 3.42),
                    ('600-1000', 7.67, 5.53, 3.56), ('>1000', 13.33, 12.00, 4.33)):
    ck('Hbeta cat %s (%s)' % (zb, c), c, float(V['hb30k cat %s' % zb]), tol=0.006)
    ck('Hbeta local Ha %s (%s)' % (zb, a), a, float(V['hb30k locA %s' % zb]), tol=0.006)
    ck('Hbeta local Hb %s (%s)' % (zb, b), b, float(V['hb30k locB %s' % zb]), tol=0.006)
ck('Hbeta/Halpha ratio low (0.60)', 0.60, float(V['hb30k ratio min']), tol=0.006)
ck('Hbeta/Halpha ratio high (0.69)', 0.69, float(V['hb30k ratio max']), tol=0.006)
ck('Hbeta/Halpha ratio >1000 (0.36)', 0.36, float(V['hb30k ratio far']), tol=0.006)
ck('Hbeta in100 with S/N>=20 (29.12)', 29.12, float(V['hb20k locB in100']), tol=0.006)
# 外两档活动率本身
ck("bin fraction 400-600 (6.65)", 6.65, float(V["bin f 400-600"]), tol=0.006)
ck("bin fraction 600-1000 (7.76)", 7.76, float(V["bin f 600-1000"]), tol=0.006)

# 打字漂移能吃掉的份额（C1 -> C2 与 原始 -> C2 两种口径的包络）
ck("drift share lower (20)", 20.0, float(V["drift share C1C2 min"]), tol=0.5)
ck("drift share upper (90)", 90.0, float(V["drift share raw max"]), tol=0.5)
ck("drift share reweighting only dM2 (26)", 26.0, float(V["drift share C1C2 dM2"]), tol=0.6)
ck("drift share reweighting only dM5 (41)", 41.0, float(V["drift share C1C2 dM5"]), tol=0.6)


def ck_gt(item, floor, got):
    ok = float(got) > float(floor)
    rows.append(dict(item=item, body=floor, computed=round(float(got), 4), diff=np.nan,
                     status="PASS" if ok else "STALE"))

# 加噪位移（正文 at most 0.16 percentage points）
ck("noise-induced shift catalogue (0.16)", 0.16, float(V["noise shift catalogue max pp"]), tol=0.006)
ck("noise-induced shift local (0.13)", 0.13, float(V["noise shift local max pp"]), tol=0.006)

# 指数只在 |Z| < 400 pc 内成立：加外档后 H 与 RMS 都变差
ck_gt("H grows when 400-600 pc is added", 0.0, float(V["exp H 4bins"]) - float(V["exp H 3bins"]))
ck_gt("H grows further with all six bins", 0.0, float(V["exp H 6bins"]) - float(V["exp H 4bins"]))
ck_gt("RMS grows when 400-600 pc is added", 0.0, float(V["exp rms 4bins"]) - float(V["exp rms 3bins"]))
ck_gt("RMS grows further with all six bins", 0.0, float(V["exp rms 6bins"]) - float(V["exp rms 4bins"]))

# Hbeta/Halpha 逐档（正文 0.60--0.69）
for zb, v in (("<100", 0.69), ("100-200", 0.66), ("200-400", 0.68), ("400-600", 0.60), ("600-1000", 0.64)):
    ck("Hbeta/Halpha %s (%s)" % (zb, v), v, float(V["hb ratio %s" % zb]), tol=0.006)

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "step2_check.csv"), index=False, encoding="utf-8-sig")
n = res.status.value_counts()
print(res.to_string(index=False))
print("\n  合计 %d 项：%s" % (len(res), "，".join("%s %d" % (k, v) for k, v in n.items())))
