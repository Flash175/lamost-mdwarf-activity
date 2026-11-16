# -*- coding: utf-8 -*-
"""make_data_dictionary.py — 生成/更新 data_dictionary.csv（字段级数据字典）。

为什么要它：审稿人和复现者拿到的是一堆列，没有单位、没有口径、没有阈值，
就会把 EW(Hα) 当成连续谱归一化通量、把 zb 当成金属度。字典把这些钉死。

行为：
1. 扫描 config 里列出的每张表，读出表头与 dtype；
2. 用内置 SEED（已从复现脚本实测得到的字段说明）填充；SEED 未覆盖的列留 TODO 待人工补；
3. 已存在的字典会被保留：你手填过的 description/unit/notes 不会被覆盖，
   只有新出现的列才追加；
4. 表里已不存在的字段保留但标 STALE，提醒你确认。

手工运行：
    python tools/make_data_dictionary.py --tables a.csv b.csv --out 03_复现代码/data_dictionary.csv
"""
import argparse
import csv
import os
import pathlib
import sys

import pandas as pd

COLS = ["table", "field", "dtype", "unit", "description", "source",
        "valid_range", "missing", "notes", "status"]

TODO = "TODO"

# ---------------------------------------------------------------- SEED
# 说明统一写法：description 说“是什么”，notes 说“容易踩的坑”。
COMMON = {
    "uid": dict(unit="-", description="LAMOST 恒星唯一标识（每个恒星可有多条光谱）",
                source="LAMOST DR13 mstellar", valid_range="非空字符串",
                missing="NaN=缺失", notes="去重键：每星取 EW(Hα) 最大的一条谱"),
    "ewha": dict(unit="Å", description="Hα 等值宽，LAMOST 谱线测量值（正=发射）",
                  source="LAMOST DR13 mstellar", valid_range="(-20, 50]",
                  missing="NaN=未测", notes="活动判据：>=0.75 且 >3*ewha_err"),
    "ewha_err": dict(unit="Å", description="EW(Hα) 的 1σ 测量误差",
                     source="LAMOST DR13 mstellar", valid_range="(0, 5]",
                     missing="NaN=未测", notes="显著性判据用 3*ewha_err"),
    "ra": dict(unit="deg", description="赤经（ICRS，J2000）",
               source="LAMOST DR13 mstellar", valid_range="[0, 360)",
               missing="NaN=缺失", notes="-"),
    "dec": dict(unit="deg", description="赤纬（ICRS，J2000）",
                source="LAMOST DR13 mstellar", valid_range="[-90, 90]",
                missing="NaN=缺失", notes="-"),
    "gaia_source_id": dict(unit="-", description="Gaia DR3 source_id（19 位整数）",
                           source="LAMOST DR13 交叉", valid_range="19 位整数",
                           missing="NaN=未匹配",
                           notes="按字符串读入再转数值，直接按 float 读会丢精度、匹配率从 ~99.6% 掉到 ~25%"),
    "gaia_g_mean_mag": dict(unit="mag", description="Gaia G 波段平均星等",
                            source="Gaia DR3", valid_range="实数", missing="NaN=缺失",
                            notes="用于去掉亮天体：M_G < 6 剔除"),
    "source_id": dict(unit="-", description="Gaia DR3 source_id（19 位整数）",
                      source="Gaia DR3", valid_range="19 位整数", missing="NaN=缺失",
                      notes="与 LAMOST 表 gaia_source_id 对应"),
    "parallax": dict(unit="mas", description="视差（未去零视差偏差）",
                     source="Gaia DR3", valid_range=">0.2（本样本）", missing="NaN=缺失",
                     notes="距离取 1000/parallax [pc]；不做视差零点修正"),
    "parallax_error": dict(unit="mas", description="视差 1σ 误差",
                           source="Gaia DR3", valid_range=">0", missing="NaN=缺失",
                           notes="显著性判据：parallax/parallax_error > 5"),
    "ruwe": dict(unit="-", description="Gaia 重归一化单位权重误差（天体测量质量指标）",
                 source="Gaia DR3", valid_range="<1.4（本样本）", missing="NaN=缺失",
                 notes="ruwe>=1.4 多为双星/不良解，已剔除"),
    "radial_velocity": dict(unit="km s^-1", description="Gaia DR3 视向速度",
                            source="Gaia DR3", valid_range="|RV|<500（本样本）",
                            missing="NaN=未提供",
                            notes="6D 子样本要求该值有限；NaN 的星被排除"),
}

PER_TABLE = {
    "金属度_清洁样本_637973.csv": {
        "zb": dict(unit="pc", description="日心 |Z| 距离分档（字符串区间）",
                   source="derived: 1000/parallax × 银纬", valid_range="<100 / 100-200 / 200-400 / 400-600 / 600-1000 / >1000",
                   missing="NaN=无", notes="与 Table 2 的六档完全对应；边界左闭右开"),
        "active": dict(unit="-", description="活动判据布尔（1=活动星）",
                       source="derived: EW(Hα)>=0.75 Å 且 >3σ", valid_range="0 / 1",
                       missing="NaN=未判定", notes="阈值 0.75 Å 与 3σ 同时满足才算活动"),
    },
    "threshold_band_by_Z.csv": {
        "zbin": dict(unit="pc", description="距离分档标签", source="derived",
                     valid_range="六档标签", missing="-", notes="与主表 zb 同口径"),
        "band_070_080": dict(unit="pp", description="0.70–0.80 Å 压线带宽度（占该档比例）",
                             source="derived: threshold band scan", valid_range="[0,100]",
                             missing="-", notes="论文 Table 2 的 band 列"),
        "bias_max_pp": dict(unit="pp", description="压线带导致的最大偏差（半宽）",
                            source="derived", valid_range="[0,100]", missing="-",
                            notes="论文 Table 2 的 half 列；400-600 档脚本用 2.53，论文表格印 2.54 属舍入笔误"),
        "band_050_100": dict(unit="pp", description="0.50–1.00 Å 宽压线带比例",
                             source="derived", valid_range="[0,100]", missing="-",
                             notes="论文 Table 2 的 wide 列"),
        "bias_wide_pp": dict(unit="pp", description="宽压线带对应的偏差",
                             source="derived", valid_range="[0,100]", missing="-",
                             notes="论文 Table 2 的 wb 列"),
    },
}

# 不进入字典的表（日志/中间垃圾）
SKIP_NAMES = {"data_dictionary.csv"}


# ---------------------------------------------------------------- 读取表头
def read_header(path: pathlib.Path):
    """返回 [(field, dtype)]。大表只读前 200 行，parquet 读 schema。"""
    if not path.exists():
        raise FileNotFoundError(f"表不存在: {path}")
    if path.suffix.lower() in {".parquet", ".pq"}:
        import pyarrow.parquet as pq
        sch = pq.read_schema(path)
        return [(n, str(sch.field(n).type)) for n in sch.names]
    df = pd.read_csv(path, nrows=200, low_memory=False,
                     encoding="utf-8-sig", comment="#")
    return [(c, str(df[c].dtype)) for c in df.columns]


def lookup(table_name: str, field: str):
    """先查本表 SEED，再查通用 SEED；都没有返回 None。"""
    per = PER_TABLE.get(table_name, {})
    if field in per:
        return per[field], f"SEED:{table_name}"
    if field in COMMON:
        return COMMON[field], "SEED:common"
    return None, None


# ---------------------------------------------------------------- 主逻辑
def build(tables, out_path: pathlib.Path):
    old = {}
    if out_path.exists():
        df_old = pd.read_csv(out_path, encoding="utf-8-sig", dtype=str).fillna("")
        for _, r in df_old.iterrows():
            old[(r["table"], r["field"])] = {c: r.get(c, "") for c in COLS}

    rows, seen = [], set()
    for t in tables:
        p = pathlib.Path(t)
        if p.name in SKIP_NAMES:
            continue
        name = p.name
        try:
            header = read_header(p)
        except FileNotFoundError as e:
            print(f"[WARN] {e}（跳过该表，解压 Dataset 后重跑即可）", file=sys.stderr)
            continue
        for field, dtype in header:
            seen.add((name, field))
            seed, origin = lookup(name, field)
            prev = old.get((name, field), {})
            row = dict(table=name, field=field, dtype=dtype)
            # 人工填过的内容优先保留
            prev_filled = prev and prev.get("description", "") not in ("", TODO)

            if prev_filled:
                row["unit"] = prev.get("unit", "") or "-"
                row["description"] = prev["description"]
                row["source"] = prev.get("source", "") or (origin or "")
                row["valid_range"] = prev.get("valid_range", "") or "-"
                row["missing"] = prev.get("missing", "") or "-"
                row["notes"] = prev.get("notes", "") or "-"
                # 保留人工写入的 status（OK(measured) / OK(partial) 等），不要一律降级成 OK(hand)
                row["status"] = prev.get("status", "") or "OK(hand)"
            elif seed:
                row.update(
                    unit=seed["unit"], description=seed["description"],
                    source=seed["source"], valid_range=seed["valid_range"],
                    missing=seed["missing"], notes=seed["notes"], status="OK(seed)",
                )
            else:
                row.update(unit=TODO, description=TODO, source=TODO,
                           valid_range=TODO, missing=TODO, notes=TODO, status=TODO)
            rows.append(row)

    # 旧字典里存在、但新表头里没有的字段 -> STALE
    for (name, field), rec in old.items():
        if (name, field) not in seen:
            rec = dict(rec)
            rec["status"] = "STALE"
            rows.append(rec)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=COLS).to_csv(out_path, index=False, encoding="utf-8-sig")

    n_ok, n_todo, n_stale = (
        sum(r["status"].startswith("OK") for r in rows),
        sum(r["status"] == TODO for r in rows),
        sum(r["status"] == "STALE" for r in rows),
    )
    print(f"字典 -> {out_path}")
    print(f"  字段总数 {len(rows)} | 已说明 {n_ok} | 待补 {n_todo} | 失效 {n_stale}")
    if n_todo:
        print("  待补字段：", ", ".join(f"{r['table']}:{r['field']}" for r in rows if r["status"] == TODO))


def main():
    try:
        sm = snakemake  # noqa: F821
        tables = [str(t) for t in sm.params.tables]
        out = pathlib.Path(str(sm.output[0]))
    except NameError:
        ap = argparse.ArgumentParser(description="生成字段级数据字典")
        ap.add_argument("--tables", nargs="+", required=True)
        ap.add_argument("--out", default="03_复现代码/data_dictionary.csv")
        a = ap.parse_args()
        tables, out = a.tables, pathlib.Path(a.out)
    build(tables, out)


if __name__ == "__main__":
    main()
