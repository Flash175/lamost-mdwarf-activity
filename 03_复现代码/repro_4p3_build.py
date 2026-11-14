
# -*- coding: utf-8 -*-
"""step1_build.py — 4.3 节的量：六量分布位移分解（线性概率模型 + 2^6 子集 Shapley）与稳健性。

输入：数据/分解样本_snrr_470548.csv（470,548 颗最终分析样本）
      数据/多因素_六因素贡献_snrr版_带误差.csv（bootstrap 16--84 区间）
      数据/多因素_渠道对齐_snrr版_降幅口径.csv（与两渠道分解对齐）
      数据/新图1_数据_2D分箱.csv（全样本降幅 10.90 / 16.03）
输出：输出/step1_values.csv、数据/多因素_稳健性_snrr版.csv
"""
import os, io, itertools, math
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(os.path.dirname(os.path.dirname(HERE)), "数据")
OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)

def load(name):
    p = os.path.join(D, name)
    lines = [l for l in open(p, encoding="utf-8-sig", errors="replace").read().splitlines() if not l.startswith("#")]
    return pd.read_csv(io.StringIO("\n".join(lines)))

vals = []
def put(k, v):
    vals.append(dict(item=k, value=float(v)))

SIX = ["F1", "F2", "F3", "F4", "F5", "F6"]
LAB = {1: "100-200", 2: "200-400"}

def factors(d):
    d = d.copy()
    d["F2"] = (pd.to_numeric(d.teff) - pd.to_numeric(d.teff).mean())/1000.0
    d["F3"] = pd.to_numeric(d.logg) - pd.to_numeric(d.logg).mean()
    d["F4"] = pd.to_numeric(d.m_h) - pd.to_numeric(d.m_h).mean()
    d["F5"] = np.log10(pd.to_numeric(d.snrr)) - np.log10(pd.to_numeric(d.snrr)).mean()
    d["F6"] = np.log10(pd.to_numeric(d.dist_pc)) - np.log10(pd.to_numeric(d.dist_pc)).mean()
    return d

def design(d):
    dum = pd.get_dummies(d.subclass, prefix="s", drop_first=True).astype(float)
    X = pd.concat([dum, d[["F2", "F3", "F4", "F5", "F6"]]], axis=1)
    blocks = {"F1": list(dum.columns)}
    for k in ("F2", "F3", "F4", "F5", "F6"):
        blocks[k] = [k]
    return X, blocks

def fit_lpm(A, y):
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    return A @ beta

def fit_logit(A, y, iters=25):
    beta = np.zeros(A.shape[1])
    for _ in range(iters):
        p = 1.0/(1.0 + np.exp(-(A @ beta)))
        w = p*(1.0 - p)
        step = np.linalg.solve((A*w[:, None]).T @ A, A.T @ (y - p))
        beta = beta + step
        if np.max(np.abs(step)) < 1e-8:
            break
    return 1.0/(1.0 + np.exp(-(A @ beta)))

def decompose(d, edges=(100, 200, 400), keys=SIX, link="lpm"):
    X, blocks = design(d)
    y = d.active.values.astype(float)
    z = d.absz.values
    b0, b1, b2 = edges
    bid = np.where(z < b0, 0, np.where(z < b1, 1, 2))
    masks = {b: (bid == b) for b in (0, 1, 2)}
    frac = {b: float((bid == b).mean()) for b in (0, 1, 2)}
    obs = {b: 100.0*(y[masks[0]].mean() - y[masks[b]].mean()) for b in (1, 2)}
    def delta(S):
        cols = [c for k in S for c in blocks[k]]
        A = np.column_stack([np.ones(len(X))] + [X[c].values for c in cols])
        pred = fit_lpm(A, y) if link == "lpm" else fit_logit(A, y)
        v = {b: pred[m].mean() for b, m in masks.items()}
        return {b: 100.0*(v[b] - v[0]) for b in (1, 2)}
    cache = {}
    for r in range(len(keys) + 1):
        for S in itertools.combinations(keys, r):
            cache[frozenset(S)] = delta(S)
    K = len(keys)
    phi = {b: {k: 0.0 for k in keys} for b in (1, 2)}
    for k in keys:
        others = [x for x in keys if x != k]
        for r in range(len(others) + 1):
            for S in itertools.combinations(others, r):
                wgt = math.factorial(len(S))*math.factorial(K - len(S) - 1)/math.factorial(K)
                dS = cache[frozenset(S)]
                dSk = cache[frozenset(tuple(sorted(S + (k,))))]
                for b in phi:
                    phi[b][k] += wgt*(dSk[b] - dS[b])
    full = cache[frozenset(keys)]
    return dict(obs=obs, phi=phi, full=full, frac=frac, n=len(y))

d0 = load("分解样本_snrr_470548.csv")
put("n sample", len(d0))
d = factors(d0)

sub = d0.subclass.value_counts()
for s in ("dM0", "dM1", "dM2", "dM3", "dM4", "dM5", "dM6"):
    put("subtype n %s" % s, int(sub.get(s, 0)))
put("dM0-dM6 share", 100.0*float(sum(sub.get(s, 0) for s in ("dM0","dM1","dM2","dM3","dM4","dM5","dM6")))/len(d0))

M = decompose(d)
for b in (0, 1, 2):
    lab = "<100" if b == 0 else LAB[b]
    put("n %s" % lab, M["frac"][b]*len(d))
    put("pct %s" % lab, 100.0*M["frac"][b])
for b in (1, 2):
    lab = LAB[b]
    put("obs decline %s" % lab, M["obs"][b])
    for k in SIX:
        put("contrib %s %s" % (lab, k), -M["phi"][b][k])
    put("meas %s" % lab, -(M["phi"][b]["F5"] + M["phi"][b]["F6"]))
    expl = abs(sum(M["phi"][b].values()))
    put("explained %s" % lab, expl)
    put("explained share %s" % lab, 100.0*expl/M["obs"][b])
    put("unexplained %s" % lab, M["obs"][b] - expl)
    put("unexplained share %s" % lab, 100.0*(M["obs"][b] - expl)/M["obs"][b])

ferr = max(abs(sum(M["phi"][b].values()) - M["full"][b]) for b in (1, 2))
put("shapley identity err", ferr)

B = load("新图1_数据_2D分箱.csv").set_index("zbin")
put("full decline 100-200", float(B.loc["<100", "frac"]) - float(B.loc["100-200", "frac"]))
put("full decline 200-400", float(B.loc["<100", "frac"]) - float(B.loc["200-400", "frac"]))

E = load("多因素_六因素贡献_snrr版_带误差.csv")
half = {}
for q in SIX + ["EXPLAINED_SUM", "UNEXPLAINED"]:
    g = E[E.quantity == q]
    half[q] = float(((g.hi84 - g.lo16)/2.0).max())
put("boot half F max", max(half[k] for k in SIX))
put("boot half explained max", half["EXPLAINED_SUM"])
put("boot half unexplained max", half["UNEXPLAINED"])

err = 0.0
for b in (1, 2):
    lab = LAB[b]
    for k in SIX:
        a = float(E[(E.bin == lab) & (E.quantity == k)].pp.iloc[0])
        err = max(err, abs(a - (-M["phi"][b][k])))
put("max |computed - archive| pp", err)

R2 = decompose(d, keys=[k for k in SIX if k != "F6"])
R5 = decompose(d, edges=(125, 250, 400))
for tag, res in (("R2", R2), ("R5", R5)):
    for b in (1, 2):
        put("%s explained %s" % (tag, LAB[b]), abs(sum(res["phi"][b].values())))
        put("%s obs %s" % (tag, LAB[b]), res["obs"][b])
    put("%s F1 %s" % (tag, LAB[1]), -res["phi"][1]["F1"])
    put("%s F1 %s" % (tag, LAB[2]), -res["phi"][2]["F1"])

hi = d[(pd.to_numeric(d.snrr) > 50) & (pd.to_numeric(d.m_h_err) < 0.2)]
R8 = decompose(hi)
put("R8 n", len(hi))
for b in (1, 2):
    put("R8 explained %s" % LAB[b], abs(sum(R8["phi"][b].values())))
    put("R8 obs %s" % LAB[b], R8["obs"][b])
    put("R8 F5 %s" % LAB[b], -R8["phi"][b]["F5"])
    put("R8 F6 %s" % LAB[b], -R8["phi"][b]["F6"])

R3own = decompose(d, link="logit")
for b in (1, 2):
    put("R3own explained %s" % LAB[b], abs(sum(R3own["phi"][b].values())))

R3_MEMO = {"100-200": 3.27, "200-400": 6.11}
put("R3 memo explained 100-200", R3_MEMO["100-200"])
put("R3 memo explained 200-400", R3_MEMO["200-400"])
put("R3 own vs memo diff 100-200", abs(abs(sum(R3own["phi"][1].values())) - R3_MEMO["100-200"]))
put("R3 own vs memo diff 200-400", abs(abs(sum(R3own["phi"][2].values())) - R3_MEMO["200-400"]))
put("logit max factor change (memo)", 0.38)

base = [abs(sum(M["phi"][1].values())), abs(sum(M["phi"][2].values()))]
chg = [abs(R3_MEMO["100-200"] - base[0]), abs(R3_MEMO["200-400"] - base[1])]
for res in (R2, R5):
    for b in (1, 2):
        chg.append(abs(abs(sum(res["phi"][b].values())) - base[b - 1]))
put("variant change logit min pp", min(chg[:2]))
put("variant change all min pp", min(chg))
put("variant change max pp", max(chg))
put("R5 F1 change max", max(abs(-R5["phi"][1]["F1"] - (-M["phi"][1]["F1"])),
                            abs(-R5["phi"][2]["F1"] - (-M["phi"][2]["F1"]))))

A = load("多因素_渠道对齐_snrr版_降幅口径.csv").set_index("bin")
for lab in ("100-200", "200-400"):
    put("within twochannel %s" % lab, float(A.loc[lab, "within_twochannel"]))
    put("within explained %s" % lab, float(A.loc[lab, "within_explained"]))
    put("comp twochannel %s" % lab, float(A.loc[lab, "comp_twochannel"]))

def expl_of(res, b):
    return abs(sum(res["phi"][b].values()))

rows = []
for b in (1, 2):
    lab = LAB[b]
    rows.append(dict(test="R3_logit", bin=lab, n=np.nan, explained_pp=R3_MEMO[lab],
                     observed_pp=M["obs"][b], f5_pp=np.nan, f6_pp=np.nan,
                     source="memo (s635): 六因素分解_第二轮结果与审核备忘_v2.md 第五节"))
    rows.append(dict(test="R3_logit_own", bin=lab, n=len(d), explained_pp=expl_of(R3own, b),
                     observed_pp=M["obs"][b], f5_pp=np.nan, f6_pp=np.nan,
                     source="本流程独立重算（逐星 IRLS logit，64 子集）"))
    rows.append(dict(test="R2_no_F6", bin=lab, n=len(d), explained_pp=expl_of(R2, b),
                     observed_pp=R2["obs"][b], f5_pp=np.nan, f6_pp=np.nan, source="本流程重算"))
    rows.append(dict(test="R5_edges_125_250", bin=lab, n=len(d), explained_pp=expl_of(R5, b),
                     observed_pp=R5["obs"][b], f5_pp=np.nan, f6_pp=np.nan, source="本流程重算"))
    rows.append(dict(test="R8_highSN", bin=lab, n=len(hi), explained_pp=expl_of(R8, b),
                     observed_pp=R8["obs"][b], f5_pp=-R8["phi"][b]["F5"], f6_pp=-R8["phi"][b]["F6"],
                     source="本流程重算（snrr>50 且 m_h_err<0.2）"))
R = pd.DataFrame(rows)[["test", "bin", "n", "explained_pp", "observed_pp", "f5_pp", "f6_pp", "source"]].round(4)
hdr = ["# 4.3 robustness: R3 logit (archived memo value + independent recomputation), R2 (no distance term),",
       "#               R5 (125/250 pc edges), R8 (high S/N subsample snrr>50 and m_h_err<0.2)",
       "# generated by 复现流程/复现流程4.3/step1_build.py (代码/repro_4p3_build.py)"]
with open(os.path.join(D, "多因素_稳健性_snrr版.csv"), "w", encoding="utf-8-sig", newline="") as f:
    f.write("\n".join(hdr) + "\n")
    R.to_csv(f, index=False)

V = pd.DataFrame(vals).round(6)
V.to_csv(os.path.join(OUT, "step1_values.csv"), index=False, encoding="utf-8-sig")
print(V.to_string(index=False)[:1800])
print("  -> 输出/step1_values.csv（%d 项）" % len(V))
