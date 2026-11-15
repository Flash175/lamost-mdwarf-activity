# 数据目录（不随本仓库分发）

脚本统一从本目录（`03_复现代码/数据/`）读取派生表与逐星表。这些文件体积较大，**未纳入代码仓库**，
请从独立的 Zenodo Dataset 记录下载后解压到此处：

- Dataset DOI：`10.5281/zenodo.24219386` — <https://doi.org/10.5281/zenodo.24219386>（all versions，解析到最新版本）

## 目录结构

解压后请保持脚本所期望的相对路径，例如：

```
数据/
├── 金属度_清洁样本_637973.csv        # 637,973 行最终逐星样本
├── 分解样本_snrr_470548.csv
├── RV变星_SB1候选_多历元.csv
├── hbeta30k_per_spectrum.csv
├── 表4.11_RV变星并集_原始重算.csv
├── 表4.12_距离精度_逐星.csv
├── 原始输入/
│   ├── lamost/dr13/dr13_mstellar_dM.csv
│   ├── gaia/gaia_dr3_mstars.csv
│   └── processed/dr13/...
└── ...（其余中间表，见各流程 README 列出的输入文件名）
```

各流程目录（`复现流程/复现流程X.Y/`）的 README.md 写明了该步骤需要的具体输入文件与输出文件。

## 第三方原始数据

LAMOST DR13 光谱目录与 Gaia DR3 天体测量属于第三方数据，不在 Dataset 内分发，
需从官方归档获取：LAMOST <https://www.lamost.org/dr13/>；Gaia <https://gea.esac.esa.int/archive/>。
