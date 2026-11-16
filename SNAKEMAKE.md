# 用 Snakemake 一条命令复现

本仓库的复现脚本本身已经能跑（`python 03_复现代码/复现流程/run_all.py`），
Snakemake 只是给它们加一层**依赖图 + 失败即停 + 并行 + 统一日志**，
让审稿人看到 `snakemake -j 4` 就能跑完主流程与全部子流程。

## 目录约定

```
Snakefile            # 规则定义（自动发现 复现流程*/run_all.py）
config.yaml          # 解释器路径 + 判据参数
envs/analysis.yml    # conda 环境（含 snakemake-minimal）
tools/run_one.py     # 统一执行器：落日志、写 rc 元信息
tools/make_data_dictionary.py
tools/check_data_dictionary.py
DATA_DICTIONARY.md   # 人读的字段口径
logs/                # 运行日志与完成标记（不入库）
```

## 前置：数据必须先就位

`03_复现代码/数据/` 只随附 README，实际表从 Dataset 记录下载解压：

1. 打开 <https://doi.org/10.5281/zenodo.24219386>（concept DOI，指向最新版本）；
2. 下载 zip 解压到 `03_复现代码/数据/`，保持 README 里写的相对路径；
3. 第三方原始表（LAMOST/Gaia）按 `03_复现代码/数据/README.md` 放 `原始输入/` 下。

数据没到位时，`snakemake -n` 仍能出计划，但执行会在第一条规则上报缺文件——这是预期行为。

## 安装与运行

```powershell
# 方式 A：conda 环境（推荐，含 snakemake）
mamba env create -f envs\analysis.yml
conda activate mdwarf-repro
cd C:\Users\lingz\Desktop\lamost-mdwarf-activity

snakemake -n                 # dry-run，先看要跑什么
snakemake -j 4               # 主流程 + 22 个子流程 + 字典校验
snakemake --dag | dot -Tpdf > dag.pdf   # 依赖图，可作补充材料
```

```powershell
# 方式 B：不装 conda，直接 pip 进 py311
C:\Users\lingz\anaconda3\envs\py311\python.exe -m pip install snakemake-minimal
C:\Users\lingz\anaconda3\envs\py311\Scripts\snakemake.exe -n
```

> Snakemake 在 Windows 原生支持较弱（规则里的 shell 片段依赖 POSIX）。
> 本仓库已把所有规则改成 `script: tools/run_one.py`（纯 Python 调用），
> 因此 Windows 上可直接跑；若仍遇到怪问题，改用 WSL2 或云主机。

## 常用目标

| 命令 | 作用 |
|---|---|
| `snakemake -j 4` | 全部 |
| `snakemake -j 4 main_gate` | 只跑主流程五步 |
| `snakemake -j 4 subflow` | 只跑子流程（22 个，可并行） |
| `snakemake -j 4 dictionary_check` | 只校验数据字典 |
| `snakemake -j 4 --forcerun step1` | 强制重跑某条规则（脚本/数据变了才会自动重跑） |
| `snakemake -j 4 --keep-going` | 某条失败后继续跑其他分支，最后统一报错 |

## 输出怎么读

- 每条规则的原始输出在 `logs/raw/<tag>.log`；
- 每条规则的完成标记是 `logs/.../<tag>.done`（JSON：`rc`、`status`、耗时、日志路径）；
- 主流程汇总在 `logs/main.done`（`PASS` 要求五步 rc 全为 0）；
- 子流程每个都有 `实测 XX/XX PASS` 的打印，rc=0 即全过。

失败定位顺序：`logs/main.done` → 对应 `logs/raw/*.log` 尾部 → 单独复跑该脚本。

## 数据字典

```powershell
# 生成/更新（已填内容不会被覆盖，只追加新列）
python tools\make_data_dictionary.py --tables ... --out 03_复现代码\data_dictionary.csv

# 校验（不许留 TODO，字段必须全覆盖）
python tools\check_data_dictionary.py --tables ... --dict 03_复现代码\data_dictionary.csv
```

`DATA_DICTIONARY.md` 与 csv 必须同源：改了判据就同步两处，否则 `dictionary_check` 会失败。

## 与论文/归档的衔接

- 论文 availability 段可加一句：
  *The full analysis pipeline is provided as a Snakemake workflow (`snakemake -j 4`), with all thresholds exposed in `config.yaml`.*
- 归档包增加 `Snakefile`、`config.yaml`、`envs/analysis.yml`、`tools/`、`DATA_DICTIONARY.md`；
- Dataset 内容有变 → 发新版本（本版本 DOI 变，论文仍引 concept DOI `10.5281/zenodo.24219386`）。
