# -*- coding: utf-8 -*-
"""check_data_dictionary.py — 校验数据字典与真实表头是否一致。

检查四件事：
1. 表里出现的字段，字典里必须有；
2. 字典里不许留 TODO（description / unit / source 不能空）；
3. 表里已不存在的字段标了 STALE（只警告，不失败）；
4. 单位不许空，无量纲必须显式写 "-"。

约定：发现问题 -> 退出码 1（Snakemake 会停在这里）；全部干净 -> 0。
"""
import argparse
import pathlib
import sys

import pandas as pd


def header_of(path: pathlib.Path):
    if path.suffix.lower() in {".parquet", ".pq"}:
        import pyarrow.parquet as pq
        return list(pq.read_schema(path).names)
    return list(pd.read_csv(path, nrows=200, low_memory=False,
                            encoding="utf-8-sig", comment="#").columns)


def check(tables, dict_path: pathlib.Path):
    problems, warnings = [], []

    if not dict_path.exists():
        print(f"[FAIL] 字典不存在: {dict_path}（先跑 data_dictionary 规则）", file=sys.stderr)
        return 1

    d = pd.read_csv(dict_path, encoding="utf-8-sig", dtype=str).fillna("")
    for col in ("table", "field", "unit", "description", "status"):
        if col not in d.columns:
            print(f"[FAIL] 字典缺少列 {col}", file=sys.stderr)
            return 1

    # ---- 逐表核对
    for t in tables:
        p = pathlib.Path(t)
        if not p.exists():
            warnings.append(f"表缺失，跳过（解压 Dataset 后再校验）：{p}")
            continue
        cols = header_of(p)
        have = set(d.loc[d["table"] == p.name, "field"])
        for c in cols:
            if c not in have:
                problems.append(f"[MISSING] {p.name}: 字段 {c} 不在字典里")
        for f in have - set(cols):
            warnings.append(f"[STALE] {p.name}: 字典里的 {f} 已不在表中")

    # ---- 内容完整性
    for _, r in d.iterrows():
        if r["status"] == "STALE":
            continue
        for k in ("description", "unit", "source"):
            v = str(r.get(k, "")).strip()
            if v == "" or v == "TODO":
                problems.append(f"[TODO] {r['table']}:{r['field']} 的 {k} 未填写")
        if str(r.get("valid_range", "")).strip() in ("", "TODO"):
            problems.append(f"[TODO] {r['table']}:{r['field']} 的 valid_range 未填写")

    # ---- 报告
    print(f"字典：{dict_path}（{len(d)} 行）")
    for w in warnings:
        print("  " + w)
    if problems:
        print(f"[FAIL] 共 {len(problems)} 处问题：")
        for p in problems:
            print("  " + p)
        return 1
    print("[PASS] 字典与表头一致，无 TODO 字段")
    return 0


def main():
    try:
        sm = snakemake  # noqa: F821
        tables = [str(t) for t in sm.params.tables]
        dict_path = pathlib.Path(str(sm.input[0]))
    except NameError:
        ap = argparse.ArgumentParser(description="校验数据字典")
        ap.add_argument("--tables", nargs="+", required=True)
        ap.add_argument("--dict", default="03_复现代码/data_dictionary.csv")
        a = ap.parse_args()
        tables, dict_path = a.tables, pathlib.Path(a.dict)
    sys.exit(check(tables, dict_path))


if __name__ == "__main__":
    main()
