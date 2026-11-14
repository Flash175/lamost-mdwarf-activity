# M-dwarf chromospheric activity as a statistical tracer of Galactic height and dynamical heating in LAMOST DR13

复现本论文正文全部表格、图件与定量结论的代码仓库（Software 归档）。

- **样本**：637,973 颗 M 矮星，LAMOST DR13 光谱 + Gaia DR3 天体测量。
- **活动判据**：EW(Hα) ≥ 0.75 Å（3σ）。
- **论文**：*M-dwarf chromospheric activity as a statistical tracer of Galactic height and dynamical heating in LAMOST DR13*, submitted to ApJS.

> 本仓库只包含**代码、小汇总表与图件**；逐星样本与中间大表随独立的 Dataset 记录分发（见下文「数据获取」）。

## 目录结构

| 目录 / 文件 | 内容 |
| --- | --- |
| `03_复现代码/` | 复现脚本：根目录 `repro_*.py` 为各流程脚本的独立副本；`复现流程/复现流程2.3 … 4.11/` 为按小节组织的 step 脚本、`输出/` 与流程 README |
| `01_表格/` | 正文各节的小汇总表（`表X.Y_*.csv`），直接对应 §X.Y |
| `04_图件/` | 正文插图源文件（PNG，部分含 PDF） |
| `MANIFEST.csv` | 全量文件清单：类别、路径、大小、sha256 前 16 位、说明 |
| `requirements.txt` | Python 第三方依赖 |
| `LICENSE` | 代码许可（MIT） |

## 环境与安装

需要 Python 3.11+（实测 3.14.6，Windows 11）。

```bash
pip install -r requirements.txt
```

第三方依赖：numpy、pandas、matplotlib、scipy、astropy、Pillow。

## 数据获取

- **本仓库不含逐星样本与中间大表**（体积与二次分发原因），它们打包在独立的 Zenodo Dataset 记录中：
  - Dataset DOI：`10.5281/zenodo.24219386` — <https://doi.org/10.5281/zenodo.24219386>（all versions，解析到最新版本）
  - 下载解压后放到 `03_复现代码/数据/`（脚本以 `__file__` 上两级定位 `数据/`），保持其中 `原始输入/lamost/...`、`原始输入/gaia/...` 等子目录结构；各流程 README 列出了具体输入文件名。
- **LAMOST DR13** 光谱目录（1,022,559 条谱 / 763,853 颗星）与 **Gaia DR3** 天体测量属于第三方数据，本仓库不分发。官方入口：
  - LAMOST DR13：<https://www.lamost.org/dr13/>
  - Gaia DR3：<https://gea.esac.esa.int/archive/>
- 本仓库不含任何第三方论文 PDF（版权原因）。

## 如何复现

每个小节独立运行。例如复现 §4.10（金属度通道）：

```bash
cd 03_复现代码/复现流程/复现流程4.10
python run_all.py
```

每个流程目录的 `README.md` 写明：输入文件、输出文件、比对项数与实测通过数。例如：

- `复现流程4.10`：89 个量、52 项正文比对、Table 5 的 54 格逐格核对，以及与论文图 9 逐像素一致的重绘图（52/52 + 54/54 PASS）。
- `复现流程4.11`：从原始星表重跑的 RV 变星链与 109+8 项比对。

## 主要结果（与论文一致）

- 活动比例随高度：<100 pc 23.81% → 100–200 pc 12.90% → 200–400 pc 7.78% → 400–600 pc 6.65% → 600–1000 pc 7.76%。
- 控制 S/N ≥ 15 后，400 pc 以外两档为 3.84% / 5.45%（主要为选择效应）。
- 82,518 颗六维子样本中，活动星在该高度运动学上更冷（25 pc 内 σ(vz) 低 6.08±0.23 km/s，200 pc 以外符号反转）；611,094 条 LAMOST 视向速度把间隔延伸到 125 pc。
- 零结果：dM6–dM8 无转折；R′Hα–Rossby 斜率无 |Z| 依赖（|Δβ| = 0.037 ± 0.177）；首 kpc [Fe/H] 下降 0.13 dex。

## 许可与引用

- **代码**：MIT License（见 `LICENSE`）。
- **随仓库分发的表与图**：CC BY 4.0；逐星/中间大表的许可见 Dataset 记录。
- 若使用本仓库，请引用论文与以下两条归档记录：
  - 本 Software 存档（Zenodo，由 GitHub release v1.0.0 联动生成）：DOI `10.5281/zenodo.24223745` — <https://doi.org/10.5281/zenodo.24223745>（all versions；本版本 `10.5281/zenodo.24223746` — <https://doi.org/10.5281/zenodo.24223746>）。
  - Dataset 存档（all versions）：`10.5281/zenodo.24219386` — <https://doi.org/10.5281/zenodo.24219386>（本版本 `10.5281/zenodo.24223938` — <https://doi.org/10.5281/zenodo.24223938>）。
  - 论文 DOI：待注册。

## 版本

version 1（与投稿版本一致）。正文 32 页、40 条参考文献、12 张图、10 张表。
