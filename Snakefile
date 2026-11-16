# -*- coding: utf-8 -*-
"""Snakefile — 一条命令复现论文全部可机读数值。

    snakemake -n              # 只看计划（dry-run）
    snakemake -j 4            # 执行全部（主流程 5 步 + 22 个子流程 + 字典校验）
    snakemake -j 4 subflow    # 只跑子流程
    snakemake --dag | dot -Tpdf > dag.pdf

设计取舍：
- 所有 Python 用 config.python 指定的解释器（py311），不用 PATH 里的 python；
- 每条规则都通过 tools/run_one.py 执行，统一日志与 rc 语义；
- 输出是 logs/*.done 元信息 JSON，而不是仓库里已存在的 输出/*.csv——
  否则 Snakemake 看到旧结果就判定“已完成”而跳过，等于没跑。
- 想强制重跑：snakemake -j 4 --forcerun step1
"""
import json
import os
import re
from pathlib import Path

configfile: "config.yaml"

PYTHON = str(config["python"])
CODE = "03_复现代码"
FLOW = os.path.join(CODE, "复现流程")
DATA = os.path.join(CODE, "数据")
LOGS = "logs"
DICT_CSV = os.path.join(CODE, "data_dictionary.csv")
RUNNER = "tools/run_one.py"

# 字典要覆盖的表（解压 Dataset 后这些文件才会出现）
TABLES = [
    os.path.join(DATA, "金属度_清洁样本_637973.csv"),
    os.path.join(DATA, "threshold_band_by_Z.csv"),
    os.path.join(DATA, "多因素_六因素贡献_snrr版_降幅口径.csv"),
    os.path.join(DATA, "shift_share.csv"),
    os.path.join(DATA, "原始输入", "lamost", "dr13", "dr13_mstellar_dM.csv"),
    os.path.join(DATA, "原始输入", "gaia", "gaia_dr3_mstars.csv"),
    os.path.join(DATA, "原始输入", "processed", "dr13", "fig5_sample_strict_v4.csv"),
]

# 主流程五步：输出 CSV 只为给人看，Snakemake 的完成标记是 .done
STEPS = [
    ("step1", "step1_sample_chain.py",   "step1_chain.csv"),
    ("step2", "step2_bins_and_band.py",  "step2_table2.csv"),
    ("step3", "step3_decompositions.py", "step3_decomp.json"),
    ("step4", "step4_break_scan.py",     "step4_rss_profile.csv"),
    ("step5", "step5_check_all.py",      "step5_check.csv"),
]


def _natural(s):
    """把 '复现流程4.10' 排到 4.11 之前，而不是字典序 4.10 < 4.2。"""
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", s) if x != ""]


# 自动发现所有子流程目录（每个都有 run_all.py）
SUBFLOWS = sorted(
    [
        p.parent.name
        for p in Path(FLOW).glob("复现流程*/run_all.py")
        if p.parent.name != "复现流程"
    ],
    key=_natural,
)


rule all:
    input:
        f"{LOGS}/main.done",
        expand(f"{LOGS}/subflows/{{sf}}.done", sf=SUBFLOWS),
        f"{LOGS}/dictionary_check.done",


# --------------------------------------------------------------- 主流程
rule step1:
    input: DATA
    output: f"{LOGS}/steps/step1.done"
    params:
        py=PYTHON, script=os.path.join(FLOW, "step1_sample_chain.py"),
        cwd=FLOW, tag="step1"
    script: RUNNER


rule step2:
    input: f"{LOGS}/steps/step1.done", DATA
    output: f"{LOGS}/steps/step2.done"
    params:
        py=PYTHON, script=os.path.join(FLOW, "step2_bins_and_band.py"),
        cwd=FLOW, tag="step2"
    script: RUNNER


rule step3:
    input: f"{LOGS}/steps/step2.done", DATA
    output: f"{LOGS}/steps/step3.done"
    params:
        py=PYTHON, script=os.path.join(FLOW, "step3_decompositions.py"),
        cwd=FLOW, tag="step3"
    script: RUNNER


rule step4:
    input: f"{LOGS}/steps/step3.done", DATA
    output: f"{LOGS}/steps/step4.done"
    params:
        py=PYTHON, script=os.path.join(FLOW, "step4_break_scan.py"),
        cwd=FLOW, tag="step4"
    script: RUNNER


rule step5:
    """汇总比对（61 项 PASS/FAIL）。rc=2 表示有 FAIL，Snakemake 会停在这里。"""
    input:
        f"{LOGS}/steps/step1.done", f"{LOGS}/steps/step2.done",
        f"{LOGS}/steps/step3.done", f"{LOGS}/steps/step4.done",
        DATA
    output: f"{LOGS}/steps/step5.done"
    params:
        py=PYTHON, script=os.path.join(FLOW, "step5_check_all.py"),
        cwd=FLOW, tag="step5"
    script: RUNNER


# --------------------------------------------------------------- 子流程
rule subflow:
    """22 个子流程（2.3–4.11），彼此独立，可并行。"""
    input: f"{FLOW}/{{sf}}/run_all.py"
    output: f"{LOGS}/subflows/{{sf}}.done"
    params:
        py=PYTHON,
        script=f"{FLOW}/{{sf}}/run_all.py",
        cwd=f"{FLOW}/{{sf}}",
        tag="sub_{sf}"
    script: RUNNER


# --------------------------------------------------------------- 数据字典
rule data_dictionary:
    input: DATA
    output: f"{LOGS}/steps/data_dictionary.done"
    params:
        py=PYTHON,
        script="tools/make_data_dictionary.py",
        cwd=".",
        tag="data_dictionary",
        extra=["--tables"] + TABLES + ["--out", DICT_CSV]
    script: RUNNER


rule dictionary_check:
    """字典必须覆盖所有字段且无 TODO，才算归档可以更新 Dataset 版本。"""
    input: f"{LOGS}/steps/data_dictionary.done"
    output: f"{LOGS}/dictionary_check.done"
    params:
        py=PYTHON,
        script="tools/check_data_dictionary.py",
        cwd=".",
        tag="dictionary_check",
        extra=["--tables"] + TABLES + ["--dict", DICT_CSV]
    script: RUNNER


# --------------------------------------------------------------- 汇总
rule main_gate:
    input:
        f"{LOGS}/steps/step1.done", f"{LOGS}/steps/step2.done",
        f"{LOGS}/steps/step3.done", f"{LOGS}/steps/step4.done",
        f"{LOGS}/steps/step5.done"
    output: f"{LOGS}/main.done"
    run:
        metas = [json.loads(Path(f).read_text(encoding="utf-8")) for f in input]
        ok = all(m["rc"] == 0 for m in metas)
        Path(output[0]).write_text(
            json.dumps({"status": "PASS" if ok else "FAIL", "steps": metas},
                       ensure_ascii=False, indent=2),
            encoding="utf-8",
        )
        if not ok:
            raise SystemExit("主流程未全部通过，详见 logs/raw/*.log")
