# -*- coding: utf-8 -*-
"""step3_subtypes.py — §3.1 的逐型近/远比例、倍数与星数。

输入：../../数据/subtype_with_giant_cut.csv、../../数据/subtype_retype.csv
输出：输出/step3_subtypes.csv
"""
import os
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(os.path.dirname(HERE))
D = os.path.join(ROOT, "数据"); OUT = os.path.join(HERE, "输出"); os.makedirs(OUT, exist_ok=True)
RD = dict(low_memory=False, encoding="utf-8-sig", comment="#")

cut = pd.read_csv(os.path.join(D, "subtype_with_giant_cut.csv"), **RD).set_index("subtype")
rt = pd.read_csv(os.path.join(D, "subtype_retype.csv"), **RD).set_index("subtype")
subs = ["dM0","dM1","dM2","dM3","dM4","dM5","dM6"]
out = pd.DataFrame({
    "subtype": subs,
    "near_cut_pct": [round(float(cut.loc[s,"near_cut"]),3) for s in subs],
    "far_cut_pct":  [round(float(cut.loc[s,"far_cut"]),3) for s in subs],
    "ratio_cut":    [round(float(cut.loc[s,"ratio_cut"]),3) for s in subs],
    "n_inner": [int(rt.loc[s,"n_inner_label"]) for s in subs],
    "n_far":   [int(rt.loc[s,"n_far_label"]) for s in subs],
})
out.to_csv(os.path.join(OUT, "step3_subtypes.csv"), index=False, encoding="utf-8-sig")
print(out.to_string(index=False))
print("  ratio range: %.2f - %.2f" % (out.ratio_cut.min(), out.ratio_cut.max()))
print("  -> 输出/step3_subtypes.csv")
