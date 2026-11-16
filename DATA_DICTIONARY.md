# Data Dictionary — 逐星表与中间表的字段口径

> 本文与 `03_复现代码/data_dictionary.csv`（机器可读）同源，后者由两道工序产出：
> ① `python tools/make_data_dictionary.py --tables <表...> --out 03_复现代码/data_dictionary.csv`
> 扫表头生成骨架（内置 SEED，未知列标 TODO）；
> ② `python tools/fill_data_dictionary.py` 把人工核定的 105 条字段定义写回（只填 TODO 行，已填内容一律不覆盖）。
> 修改判据时两处必须同步，否则 `snakemake dictionary_check` 会失败。

## 数据集与代码

| 项 | 值 |
|---|---|
| 代码仓库 | <https://github.com/Flash175/lamost-mdwarf-activity> |
| 代码 Software concept DOI | [`10.5281/zenodo.24223745`](https://doi.org/10.5281/zenodo.24223745) |
| 大表 Dataset concept DOI | [`10.5281/zenodo.24219386`](https://doi.org/10.5281/zenodo.24219386) |
| 第三方原始数据 | LAMOST DR13 <https://www.lamost.org/dr13/>；Gaia DR3 <https://gea.esac.esa.int/archive/> |

## 关键判据（论文正文写作时必须与这里一致）

| 判据 | 取值 | 出处 |
|---|---|---|
| 活动星 | `EW(Hα) >= 0.75 Å` **且** `EW(Hα) > 3 × ewha_err` | §2，Table 1 step 3 |
| EW 有效区间 | `-20 < ewha <= 50` Å；`0 < ewha_err <= 5` Å | Table 1 steps 2–3 |
| 每星一条谱 | 按 `ewha` 降序取 `uid` 首个（最大 EW） | Table 1 step 6 |
| Gaia 交叉 | 按 `gaia_source_id ↔ source_id` 精确匹配（非锥形搜索） | Table 1 step 4 |
| 视差质量 | `parallax > 0.2 mas`，`parallax/parallax_error > 5`，`ruwe < 1.4` | Table 1 step 5 |
| 距离 | `d = 1000 / parallax` pc（不做视差零点修正） | §2 |
| 日心 \|Z\| | `absz = d × |sin b|`，保留 `absz <= 3000 pc` | Table 1 step 7 |
| 亮天体剔除 | `M_G = G - 5log10(d) + 5`，剔除 `M_G < 6` | Table 1 step 8 |
| 六维子样本 | Gaia RV 有限且 `|RV| < 500 km/s`、`|Z_gc| < 3000 pc` | Table 1 step 9 |
| 距离分档 | `<100 / 100–200 / 200–400 / 400–600 / 600–1000 / >1000` pc（左闭右开） | Table 2 |

## 样本链（Table 1，实测 61/61 PASS）

```
1,022,559 → 1,018,728 → 1,016,104 → 997,156 → 842,150 → 639,711 → 639,702 → 637,973 → 82,518
```

## 主表字段（`金属度_清洁样本_637973.csv`）

完整字段以 `data_dictionary.csv` 为准；下表是关键列，SEED 已内置。

| 字段 | 单位 | 含义 | 备注 |
|---|---|---|---|
| `uid` | – | LAMOST 恒星唯一标识 | 同一恒星可有多条光谱 |
| `ewha` | Å | Hα 等值宽（正=发射） | 活动判据主量 |
| `ewha_err` | Å | EW(Hα) 的 1σ 误差 | 显著性判据用 3× |
| `ra`, `dec` | deg | 赤经/赤纬（ICRS J2000） | – |
| `gaia_source_id` | – | Gaia DR3 source_id（19 位） | **必须按字符串读**，否则匹配率从 ~99.6% 掉到 ~25% |
| `gaia_g_mean_mag` | mag | Gaia G 星等 | 用于 M_G 剔除 |
| `parallax`, `parallax_error` | mas | 视差及其误差 | 距离来源 |
| `ruwe` | – | 天体测量质量指标 | `ruwe < 1.4` 剔除双星/不良解 |
| `radial_velocity` | km s⁻¹ | Gaia DR3 视向速度 | 6D 子样本要求有限 |
| `zb` | pc | 日心 \|Z\| 距离分档（字符串） | 与 Table 2 六档一一对应 |
| `active` | – | 活动布尔（0/1） | 0.75 Å 与 3σ 同时满足 |

## 字段覆盖状态（2026-11-16 全部落定）

`03_复现代码/data_dictionary.csv` 共 **136 行 / 7 张表**，无 TODO。status 分三档：

| status | 条数 | 含义 |
|---|---|---|
| `OK(seed)` | 31 | 生成器内置 SEED，来自复现脚本实测字段口径 |
| `OK(measured)` | 94 | 按真实表统计（前 20 万行、剔除 -9999/-999 哨兵）核定取值范围后填写 |
| `OK(partial)` | 11 | 列名与来源可确认，但官方物理定义未在仓库内找到依据；已注明"本流程未使用"，不编造物理含义 |

`OK(partial)` 的 11 条：`zeta`、`zeta_err`、`type`、`na`（LAMOST 参数表）、`offsets`、`offsets_v`（光纤偏移）、`chi`（上游周期目录统计量）、`w1_pct`、`w2_pct`、`f1_pct`、`f2_pct`（4.3 分解中间量）。拿到 LAMOST DR13 / Gaia DR3 官方数据字典或 4.3 脚本原文后，可直接替换这 11 条并升级 status。

校验与回归：

```powershell
snakemake -j 1 -- logs/dictionary_check.done      # 应为 PASS rc=0
snakemake -j 1 -R data_dictionary -- logs/dictionary_check.done   # 生成器强制重跑, 验证不覆盖手填内容
```

实测：`[dictionary_check] PASS rc=0`；生成器重跑后字典与人工版逐格 **0 差异**、TODO **0 处**。

缺失表示法：LAMOST 表数值列为 `-9999`（PS1 星等与若干测光列为 `-999`），其余表为清洗后无缺失。

## 中间表

| 表 | 关键字段 | 说明 |
|---|---|---|
| `threshold_band_by_Z.csv` | `zbin`, `band_070_080`, `bias_max_pp`, `band_050_100`, `bias_wide_pp` | Table 2 压线带四列；400–600 pc 档脚本用 2.53，论文表格印 2.54 属舍入笔误（把已舍入的 5.07 再除以 2） |
| `fig5_sample_strict_v4.csv` | – | Table 10 / A.6 断点扫描输入（3,340 条谱） |
| `图5_逐星样本_严格3340行_⑤.csv` | – | **同 3,340 行但 τ 口径不同，勿与上表混用**；误用会得 RSS_min = 1132.40 @ 0.1370（正确值 1132.19 @ 0.1322） |

## 已知限制与使用禁忌

1. **空值不是零**：`radial_velocity`/`parallax` 的 NaN 表示未测，未被当作 0 参与统计，分析时必须先剔除。
2. **`zb` 是字符串标签**，不是数值；数值化请用映射表，避免字典序排序错误。
3. **距离假设**：`d = 1000/parallax` 未做视差零点修正、未用贝叶斯距离；远距离（>1 kpc）绝对误差随视差误差放大。
4. **版本一致性**：Dataset 记录若更新内容需发新版本（本版本 DOI 变，concept DOI `10.5281/zenodo.24219386` 不变），论文引用始终用 concept DOI。
5. **第三方数据不分发**：LAMOST 光谱与 Gaia 天体测量需自行从官方归档获取，本数据集只含派生表与逐星结果。

---
更新：2026-11-16
