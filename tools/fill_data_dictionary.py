#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""把人工核定的字段定义写回 data_dictionary.csv（只填 status=TODO 的行，已填内容一律不动）。

用法:
    python tools/fill_data_dictionary.py                 # 就地补全
    python tools/fill_data_dictionary.py --dry-run       # 只报告将要填哪些字段

内容来源分三档，写在 status 里:
    OK(measured)  取值范围来自对 03_复现代码/数据 下真实表的前 20 万行统计（已剔除 -9999/-999 哨兵）
    OK(partial)   列名/表名可确认，但官方物理定义未在仓库内找到依据；已标注"本流程未使用"，
                  不编造物理含义，后续可替换为 LAMOST/Gaia 官方数据字典的原文
    OK(seed)      已有内容，脚本不覆盖
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DICT_CSV = REPO / "03_复现代码" / "data_dictionary.csv"

MAIN = "金属度_清洁样本_637973.csv"
THR = "threshold_band_by_Z.csv"
MULTI = "多因素_六因素贡献_snrr版_降幅口径.csv"
SHIFT = "shift_share.csv"
LAMOST = "dr13_mstellar_dM.csv"
GAIA = "gaia_dr3_mstars.csv"
FIG5 = "fig5_sample_strict_v4.csv"

# (table, field, unit, description, source, valid_range, notes, status)
FILLS = [
    # ---------- 主表 金属度_清洁样本_637973.csv ----------
    (MAIN, "subclass", "-",
     "LAMOST 光谱亚型（dM0–dM9）", "LAMOST DR13 mstellar（光谱分类）",
     "dM0–dM9（10 类）", "主表含 dM7–dM9；3.4 冷端分析只取 dM5–dM8", "OK(measured)"),
    (MAIN, "teff", "K",
     "LAMOST 有效温度（恒星参数流水线）", "LAMOST DR13 mstellar",
     "[3114, 4138] K", "Teff 轴与亚型匹配检查使用；-9999 为缺失", "OK(measured)"),
    (MAIN, "logg", "dex",
     "LAMOST 表面重力 log g", "LAMOST DR13 mstellar",
     "[-0.12, 5.09] dex", "低端（<3）为拟合失败值，本流程不据此筛选", "OK(measured)"),
    (MAIN, "m_h", "dex",
     "LAMOST 金属度 [Fe/H]", "LAMOST DR13 mstellar",
     "[-0.83, 0.58] dex", "阈值分档按 m_h；可用性由 mh_ok 标记", "OK(measured)"),
    (MAIN, "m_h_err", "dex",
     "m_h 的 1σ 测量误差", "LAMOST DR13 mstellar",
     "[0.03, 0.70] dex", "mh_ok 判据参考该项", "OK(measured)"),
    (MAIN, "snru", "-",
     "u 波段信噪比", "LAMOST DR13 mstellar",
     "[0, 999]", "压线带分析用 log S/N；-9999 为缺失", "OK(measured)"),
    (MAIN, "MG", "mag",
     "绝对星等 M_G = G - 5log10(d) + 5", "derived: Gaia G + 1000/parallax",
     ">= 6（实测 6.00–16.55 mag）", "M_G < 6 的亮天体已剔除", "OK(measured)"),
    (MAIN, "dist_pc", "pc",
     "日心距离 d = 1000/parallax（不做零点修正）", "derived: Gaia DR3 parallax",
     "(0, 3000]（实测 10.2–2195 pc）", "parallax > 0.2 mas 且 parallax/parallax_error > 5", "OK(measured)"),
    (MAIN, "absz", "pc",
     "日心银盘高度 |Z| = d × |sin b|", "derived: d 与银纬",
     "[0, 3000]（实测 0.003–1594 pc）", "用于六档 |Z| 分箱", "OK(measured)"),
    (MAIN, "mh_ok", "-",
     "金属度可用标志（True=该星 m_h 可进入阈值分档）", "derived: m_h 与 m_h_err 有效性",
     "True / False", "阈值带分析排除 False", "OK(measured)"),

    # ---------- threshold_band_by_Z.csv ----------
    (THR, "n", "-",
     "该 |Z| 档内的星数", "derived: 主表按 zb 计数",
     ">0（实测 1,231–191,179）", "六档左闭右开", "OK(measured)"),
    (THR, "act_pct", "%",
     "该档活动星占比（活动判据见 ewha / active）", "derived: 主表 active 按 zb 汇总",
     "[0, 100]（实测 6.65–23.81）", "论文 Table 2 的 act 列", "OK(measured)"),
    (THR, "band_065_085", "%",
     "0.65–0.85 Å 压线带宽度（占该档比例，阈值敏感性）", "derived: threshold band scan",
     "[0, 100]（实测 1.86–3.08）", "与 0.70–0.80 带联读", "OK(measured)"),

    # ---------- 多因素_六因素贡献_snrr版_降幅口径.csv ----------
    (MULTI, "bin", "pc",
     "距离档（相邻两档对比的目标档）", "derived: 主表 zb",
     "100-200 / 200-400", "-", "OK(measured)"),
    (MULTI, "factor", "-",
     "六因素编号：F1 亚型、F2 Teff、F3 log g、F4 [Fe/H]、F5 log S/N、F6 log d",
     "derived: 4.3 多因素分解", "F1–F6",
     "降幅口径按 snrr 版", "OK(measured)"),
    (MULTI, "contribution_to_decline_pp", "pp",
     "该因素对活动率降幅的贡献（百分点）", "derived: 4.3 分解",
     "实测 -0.84–16.09", "负值表示该因素反向（提升活动率）", "OK(measured)"),
    (MULTI, "share_of_observed_pct", "%",
     "贡献占观测总降幅的百分比", "derived: 4.3 分解",
     "实测 -5.2–100", "可为负", "OK(measured)"),
    (MULTI, "share_of_explained_pct", "%",
     "贡献占已解释降幅的百分比", "derived: 4.3 分解",
     "实测 -13.59–100", "分母为各因素贡献之和", "OK(measured)"),

    # ---------- shift_share.csv ----------
    (SHIFT, "transition", "-",
     "距离档迁移方向（源档>目标档）", "derived: 主表 zb",
     "如 200-400>400-600", "-", "OK(measured)"),
    (SHIFT, "subtype", "-",
     "亚型（dM0–dM6）", "derived: LAMOST subclass",
     "dM0–dM6", "按亚型分解率/组成效应", "OK(measured)"),
    (SHIFT, "rate_pp", "pp",
     "亚型内活动率变化（率效应）", "derived: Simpson 分解",
     "实测 -0.12–1.03", "率效应项", "OK(measured)"),
    (SHIFT, "comp_pp", "pp",
     "亚型构成变化（组成效应）", "derived: Simpson 分解",
     "实测 -1.08–1.00", "组成效应项", "OK(measured)"),
    (SHIFT, "inter_pp", "pp",
     "率与组成的交互项", "derived: Simpson 分解",
     "实测 -0.38–0.25", "交互项", "OK(measured)"),
    (SHIFT, "w1_pct", "%",
     "该亚型在源档的占比（权重）", "derived: 分解中间量",
     "实测 0.004–54.05", "字段名直译 w=weight", "OK(partial)"),
    (SHIFT, "w2_pct", "%",
     "该亚型在目标档的占比（权重）", "derived: 分解中间量",
     "实测 0.004–71.76", "字段名直译 w=weight", "OK(partial)"),
    (SHIFT, "f1_pct", "%",
     "该亚型在源档的活动率", "derived: 分解中间量",
     "实测 4.70–52.63", "字段名直译 f=fraction", "OK(partial)"),
    (SHIFT, "f2_pct", "%",
     "该亚型在目标档的活动率", "derived: 分解中间量",
     "实测 0–50", "字段名直译 f=fraction", "OK(partial)"),
    (SHIFT, "n1", "-",
     "源档该亚型星数", "derived: 主表计数",
     "实测 3–65,497", "-", "OK(measured)"),
    (SHIFT, "n2", "-",
     "目标档该亚型星数", "derived: 主表计数",
     "实测 2–38,742", "-", "OK(measured)"),

    # ---------- dr13_mstellar_dM.csv ----------
    (LAMOST, "obsid", "-",
     "LAMOST 观测 ID（一次曝光中的一条谱）", "LAMOST DR13 mstellar",
     "正整数（实测 101,007–219,111,000）", "与 uid 组合定位单条谱", "OK(measured)"),
    (LAMOST, "gp_id", "-",
     "general plan ID（19 位整数）", "LAMOST DR13 mstellar",
     "19 位整数（实测 1.03e5–6.92e18）", "必须按字符串读：按 float 读会丢精度", "OK(measured)"),
    (LAMOST, "designation", "-",
     "LAMOST 目标名（J2000 坐标命名）", "LAMOST DR13 mstellar",
     "JHHMMSS.ss±DDMMSS.s", "-", "OK(measured)"),
    (LAMOST, "obsdate", "-",
     "观测日期（UTC）", "LAMOST DR13 mstellar",
     "YYYY-MM-DD（实测 328 个不同日期）", "-", "OK(measured)"),
    (LAMOST, "lmjd", "d",
     "本地修正儒略日（LAMOST 时间基准）", "LAMOST DR13 mstellar",
     "实测 55,859–56,718", "lmjd = mjd + 1", "OK(measured)"),
    (LAMOST, "mjd", "d",
     "儒略日（UTC）", "LAMOST DR13 mstellar",
     "实测 55,858–56,717", "lmjd = mjd + 1", "OK(measured)"),
    (LAMOST, "planid", "-",
     "观测计划 ID", "LAMOST DR13 mstellar",
     "如 B5591606 / F5902", "-", "OK(measured)"),
    (LAMOST, "spid", "-",
     "光谱 ID（同一天体内编号）", "LAMOST DR13 mstellar",
     "1–16", "-", "OK(measured)"),
    (LAMOST, "fiberid", "-",
     "光纤编号", "LAMOST DR13 mstellar",
     "1–250", "-", "OK(measured)"),
    (LAMOST, "ra_obs", "deg",
     "观测指向赤经（J2000）", "LAMOST DR13 mstellar",
     "[0, 360)", "与 ra 基本一致（差为指向抖动）", "OK(measured)"),
    (LAMOST, "dec_obs", "deg",
     "观测指向赤纬（J2000）", "LAMOST DR13 mstellar",
     "[-90, 90]（实测 -9.64–62.59）", "-", "OK(measured)"),
    (LAMOST, "snru", "-",
     "u 波段信噪比", "LAMOST DR13 mstellar",
     "[0, 999]；-9999 为缺失", "压线带分析用 log S/N", "OK(measured)"),
    (LAMOST, "snrg", "-",
     "g 波段信噪比", "LAMOST DR13 mstellar",
     "[0, 428]；-9999 为缺失", "-", "OK(measured)"),
    (LAMOST, "snrr", "-",
     "r 波段信噪比", "LAMOST DR13 mstellar",
     "[0, 735]；-9999 为缺失", "论文降幅口径用 snrr 版", "OK(measured)"),
    (LAMOST, "snri", "-",
     "i 波段信噪比", "LAMOST DR13 mstellar",
     "[0, 999]；-9999 为缺失", "-", "OK(measured)"),
    (LAMOST, "snrz", "-",
     "z 波段信噪比", "LAMOST DR13 mstellar",
     "[0, 999]；-9999 为缺失", "-", "OK(measured)"),
    (LAMOST, "class", "-",
     "天体分类（LAMOST 流水线）", "LAMOST DR13 mstellar",
     "STAR（本表全为恒星）", "-", "OK(measured)"),
    (LAMOST, "subclass", "-",
     "光谱亚型", "LAMOST DR13 mstellar",
     "dM0–dM9", "-", "OK(measured)"),
    (LAMOST, "z", "-",
     "红移（对恒星≈日心视向速度对应的红移）", "LAMOST DR13 mstellar",
     "实测 -0.0050–0.0044；-9999 为缺失", "非本流程主变量", "OK(measured)"),
    (LAMOST, "z_err", "-",
     "红移 1σ 误差", "LAMOST DR13 mstellar",
     "实测 1e-8–378（含异常大值）", "未参与筛选", "OK(measured)"),
    (LAMOST, "ps_id", "-",
     "测光星表源 ID（PS1，19 位整数）", "PS1（经 LAMOST 交叉）",
     "19 位整数（实测 9.64e16–1.83e17）", "必须按字符串读", "OK(measured)"),
    (LAMOST, "mag_ps_g", "mag",
     "PS1 g 波段星等", "PS1（经 LAMOST 交叉）",
     "实测 10.36–26.93；-999 为缺失", "-", "OK(measured)"),
    (LAMOST, "mag_ps_r", "mag",
     "PS1 r 波段星等", "PS1（经 LAMOST 交叉）",
     "实测 9.79–22.28；-999 为缺失", "-", "OK(measured)"),
    (LAMOST, "mag_ps_i", "mag",
     "PS1 i 波段星等", "PS1（经 LAMOST 交叉）",
     "实测 7.00–24.95；-999 为缺失", "-", "OK(measured)"),
    (LAMOST, "mag_ps_z", "mag",
     "PS1 z 波段星等", "PS1（经 LAMOST 交叉）",
     "实测 7.36–26.27；-999 为缺失", "-", "OK(measured)"),
    (LAMOST, "mag_ps_y", "mag",
     "PS1 y 波段星等", "PS1（经 LAMOST 交叉）",
     "实测 4.34–24.31；-999 为缺失", "-", "OK(measured)"),
    (LAMOST, "fibertype", "-",
     "光纤类型（Obj 目标星 / F-Std 流量标准星）", "LAMOST DR13 mstellar",
     "Obj / F-Std", "本流程只用 Obj", "OK(measured)"),
    (LAMOST, "offsets", "-",
     "光纤偏移标志（0/1）", "LAMOST DR13 mstellar",
     "0 / 1", "语义见 LAMOST DR13 数据说明；本流程未使用", "OK(partial)"),
    (LAMOST, "offsets_v", "-",
     "光纤偏移量（0–99）", "LAMOST DR13 mstellar",
     "0–99", "语义见 LAMOST DR13 数据说明；本流程未使用", "OK(partial)"),
    (LAMOST, "teff", "K",
     "LAMOST 有效温度", "LAMOST DR13 mstellar",
     "[3114, 4138] K；-9999 为缺失", "-", "OK(measured)"),
    (LAMOST, "teff_err", "K",
     "teff 的 1σ 误差", "LAMOST DR13 mstellar",
     "[11.45, 460] K；-9999 为缺失", "-", "OK(measured)"),
    (LAMOST, "logg", "dex",
     "LAMOST 表面重力 log g", "LAMOST DR13 mstellar",
     "[-0.25, 5.09] dex；-9999 为缺失", "低端为拟合失败值", "OK(measured)"),
    (LAMOST, "logg_err", "dex",
     "logg 的 1σ 误差", "LAMOST DR13 mstellar",
     "[0.02, 0.47] dex", "-", "OK(measured)"),
    (LAMOST, "m_h", "dex",
     "金属度 [Fe/H]", "LAMOST DR13 mstellar",
     "[-0.83, 0.59] dex；-9999 为缺失", "-", "OK(measured)"),
    (LAMOST, "m_h_err", "dex",
     "m_h 的 1σ 误差", "LAMOST DR13 mstellar",
     "[0.03, 0.70] dex", "-", "OK(measured)"),
    (LAMOST, "zeta", "-",
     "LAMOST ζ 参数（恒星参数拟合输出）", "LAMOST DR13 mstellar",
     "实测 -43.18–121.62（主值 0.87–1.40）",
     "具体物理定义见 LAMOST DR13 数据说明；本流程未使用该列", "OK(partial)"),
    (LAMOST, "zeta_err", "-",
     "ζ 参数的 1σ 误差", "LAMOST DR13 mstellar",
     "实测 0–6,436", "本流程未使用", "OK(partial)"),
    (LAMOST, "type", "-",
     "标志位（0/1）", "LAMOST DR13 mstellar",
     "0 / 1", "语义见 LAMOST DR13 数据说明；本流程未使用", "OK(partial)"),
    (LAMOST, "na", "-",
     "Na 指数（LAMOST 谱线指数）", "LAMOST DR13 mstellar",
     "实测 -37.27–1790.53", "本流程未使用该列", "OK(partial)"),

    # ---------- gaia_dr3_mstars.csv ----------
    (GAIA, "pmra", "mas yr^-1",
     "赤经方向自行（Gaia DR3）", "Gaia DR3",
     "实测 -294.6–304.7", "本流程未使用", "OK(measured)"),
    (GAIA, "pmdec", "mas yr^-1",
     "赤纬方向自行（Gaia DR3）", "Gaia DR3",
     "实测 -283.4–223.2", "本流程未使用", "OK(measured)"),
    (GAIA, "phot_g_mean_mag", "mag",
     "Gaia G 波段平均星等", "Gaia DR3",
     "实测 8.10–21.27", "与 LAMOST 表 gaia_g_mean_mag 同源", "OK(measured)"),
    (GAIA, "bp_rp", "mag",
     "Gaia BP–RP 颜色", "Gaia DR3",
     "实测 0.46–6.33", "本流程未使用", "OK(measured)"),
    (GAIA, "teff_gspphot", "K",
     "Gaia GSP-Phot 有效温度", "Gaia DR3",
     "实测 2628–33627（含非 M 星）", "本表为交叉全集，不用于 Teff 筛选", "OK(measured)"),

    # ---------- fig5_sample_strict_v4.csv ----------
    (FIG5, "designation", "-",
     "目标名（LAMOST J 坐标命名）", "LAMOST DR13 mstellar",
     "JHHMMSS.ss±DDMMSS.s", "-", "OK(measured)"),
    (FIG5, "spt", "-",
     "光谱型（与 subclass 同源）", "LAMOST DR13 mstellar",
     "dM0–dM6", "3.4 分组用 dM0-dM3 / dM4-dM6", "OK(measured)"),
    (FIG5, "teff", "K",
     "有效温度", "LAMOST DR13 mstellar",
     "[3114, 4100] K", "-", "OK(measured)"),
    (FIG5, "P", "d",
     "自转周期", "VSX / Fetherolf 等周期目录",
     "实测 0.114–123 d", "来源由 src 列标明", "OK(measured)"),
    (FIG5, "src", "-",
     "周期来源", "VSX / Fetherolf",
     "VSX / Fetherolf", "-", "OK(measured)"),
    (FIG5, "non_single_star", "-",
     "非单星标记（1=双星/多星）", "derived: 上游星表标志",
     "0（严格样本已全部排除）", "恒为 0 即无此类星", "OK(measured)"),
    (FIG5, "alias_flag", "-",
     "周期别名标记（True=该周期可能是别名）", "derived: 上游周期筛选",
     "False（严格样本已全部排除）", "恒为 False", "OK(measured)"),
    (FIG5, "border34", "-",
     "dM3/dM4 分类边界标记", "derived: 上游分类",
     "False（严格样本已全部排除）", "边界星不进严格样本", "OK(measured)"),
    (FIG5, "M", "M_sun",
     "恒星质量", "derived: 质量-光度关系 / 模型",
     "实测 0.276–0.668", "-", "OK(measured)"),
    (FIG5, "chi", "-",
     "周期样本的质量统计量 χ（由上游周期目录给出）", "上游周期目录",
     "实测 4e-6–3.2e-5", "本流程未用该列做筛选；具体定义见上游目录", "OK(partial)"),
    (FIG5, "Ro", "-",
     "Rossby 数 Ro = P/τ_c", "derived: P 与对流周转时标 τ_c",
     "实测 0.002–2.397", "-", "OK(measured)"),
    (FIG5, "sRo", "-",
     "Ro 的 1σ 误差", "derived: 误差传递",
     "实测 0.0003–0.30", "-", "OK(measured)"),
    (FIG5, "RHa", "s^-1",
     "R'_Ha 发射率（Hα 光度归一化到恒星表面）", "derived: EW(Hα) 与 Teff",
     "实测 0–4.54e-4", "-", "OK(measured)"),
    (FIG5, "sRHa", "s^-1",
     "R'_Ha 的 1σ 误差", "derived: 误差传递",
     "实测 0–6.6e-5", "-", "OK(measured)"),
]

BANDS = [
    ("tio5", "TiO5", "实测 -13.76–101.76"),
    ("cah2", "CaH2", "实测 -14.89–77.24"),
    ("cah3", "CaH3", "实测 -7.59–44.92"),
    ("tio1", "TiO1", "实测 -105.65–644.18"),
    ("tio2", "TiO2", "实测 -31.36–28.45"),
    ("tio3", "TiO3", "实测 -820.67–92.30"),
    ("tio4", "TiO4", "实测 -355.26–67.47"),
    ("cah1", "CaH1", "实测 -135.48–272.40"),
    ("caoh", "CaOH", "实测 -305.59–133.91"),
]
BAND_ERRS = {
    "tio5_err": "实测 0–246.06", "cah2_err": "实测 0–354.62", "cah3_err": "实测 0–156.24",
    "tio1_err": "实测 0–285,909", "tio2_err": "实测 0–193.20", "tio3_err": "实测 0–968,015",
    "tio4_err": "实测 0–153,251", "cah1_err": "实测 0–1,000,000", "caoh_err": "实测 0–156.24",
}

for _field, _label, _range in BANDS:
    FILLS.append((LAMOST, _field, "-",
                  f"{_label} 分子带指数（LAMOST 恒星参数拟合用，无量纲）",
                  "LAMOST DR13 mstellar", _range,
                  "定义见 LAMOST DR13 数据说明与 Reid et al. 1995 系列；本流程未用该列筛选",
                  "OK(measured)"))
    FILLS.append((LAMOST, _field + "_err", "-",
                  f"{_label} 指数的 1σ 误差",
                  "LAMOST DR13 mstellar", BAND_ERRS[_field + "_err"],
                  "本流程未使用", "OK(measured)"))


_LAMOST_NO_SENTINEL = {
    "obsid", "gp_id", "designation", "obsdate", "lmjd", "mjd", "planid", "spid",
    "fiberid", "ra_obs", "dec_obs", "class", "subclass", "ps_id", "fibertype",
    "offsets", "offsets_v", "ra", "dec", "uid", "gaia_source_id",
}
_LAMOST_MINUS_999 = {"mag_ps_g", "mag_ps_r", "mag_ps_i", "mag_ps_z", "mag_ps_y"}


def infer_missing(table: str, field: str) -> str:
    """按实测表推缺失表示法，不臆造。"""
    if table == LAMOST:
        if field in _LAMOST_NO_SENTINEL:
            return "无"
        if field in _LAMOST_MINUS_999:
            return "-999（哨兵值）"
        return "-9999（哨兵值）"
    if table == GAIA:
        return "无（交叉全集未过滤）"
    return "无（清洗样本无缺失）"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--dict", default=str(DICT_CSV))
    args = ap.parse_args()

    path = Path(args.dict)
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh))
        fieldnames = list(rows[0].keys())

    patch = {(t, f): (u, d, s, v, n, st) for t, f, u, d, s, v, n, st in FILLS}
    if len(patch) != len(FILLS):
        dup = len(FILLS) - len(patch)
        print(f"WARN 重复键 {dup} 条（后者覆盖前者）", file=sys.stderr)

    filled, skipped, unknown = 0, 0, []
    missing_fixed = 0
    for r in rows:
        if str(r.get("missing", "")).strip() in ("", "TODO"):
            r["missing"] = infer_missing(r["table"], r["field"])
            missing_fixed += 1
        key = (r["table"], r["field"])
        if r["status"].strip() and not r["status"].strip().upper().startswith("TODO"):
            skipped += 1
            continue
        if key not in patch:
            unknown.append(key)
            continue
        unit, desc, source, valid, notes, status = patch[key]
        r["unit"], r["description"], r["source"] = unit, desc, source
        r["valid_range"], r["notes"], r["status"] = valid, notes, status
        filled += 1

    print(f"rows={len(rows)} 已填保持={skipped} 本次填={filled} 未匹配={len(unknown)} missing 补={missing_fixed}")
    for k in unknown:
        print(f"  UNMAPPED {k}", file=sys.stderr)
    if args.dry_run:
        print("dry-run，不写盘")
        return 0 if not unknown else 1

    with path.open("w", encoding="utf-8-sig", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    print(f"已写回 {path}")
    return 0 if not unknown else 1


if __name__ == "__main__":
    raise SystemExit(main())
