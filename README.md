# 肺腺癌 Hallmark 转录程序的跨平台预后基准研究

本项目评估四类 MSigDB Hallmark 转录程序（增殖、缺氧/糖酵解、EMT/间质、免疫）在肺腺癌（LUAD）中的总生存（OS）关联、对临床及测量因素的敏感性、与经典已发表增殖签名的信息冗余度，以及由 TCGA-LUAD 冻结风险模型向外部队列迁移时的区分度与校准表现。全研究遵循 REMARK 与 TRIPOD/TRIPOD+AI 报告框架，纳入 6 个公开队列共 1,443 名合格患者（其中 1,420 名患者、469 例死亡进入临床完整病例分析）。

## 队列概览

| 队列 | 研究角色 | 检测平台 | 合格 OS | 完整病例 | 完整病例死亡数 | 中位潜在随访（月） |
|---|---|---|---:|---:|---:|---:|
| TCGA-LUAD | 发现与模型训练 | RNA-seq (STAR log2) | 503 | 485 | 177 | 29.0 |
| GSE31210 | 外部验证（GPL570） | Affymetrix 芯片 | 204 | 204 | 30 | 61.0 |
| GSE37745 | 外部验证（GPL570） | Affymetrix 芯片 | 106 | 106 | 77 | 126.4 |
| GSE50081 | 外部验证（GPL570） | Affymetrix 芯片 | 127 | 127 | 51 | 65.9 |
| GSE72094 | 扩展外部验证（HuRSTA） | Rosetta/Merck 芯片 | 398 | 393 | 111 | 29.9 |
| CPTAC-LUAD | 独立外部验证（RNA-seq） | Illumina RNA-seq (RSEM) | 105 | 105 | 23 | 22.8 |
| **合计** | 1 训练集 + 5 外部队列 | 芯片 + RNA-seq | **1,443** | **1,420** | **469** | — |

## 主要结果

- **跨平台关联与外部汇总**：在五队列（*k* = 5，*N* = 935，292 例死亡）两阶段 REML-mKH 随机效应合并中，增殖程序校正分期、年龄和性别后的 HR 为 **1.35（95% CI 1.15–1.59）**，Holm 校正后 *P* = **0.027**（*I*² = 0%）；原三芯片队列与扩展四芯片队列点估计分别为 1.32 与 1.34。在五队列一阶段队列分层 IPD 模型中，增殖与缺氧 HR 分别为 **1.35** 与 **1.40**（Holm 校正后 *P* 均为 **4.4 × 10⁻⁷**），而 EMT/间质与免疫程序未呈现稳定的合并关联。[两阶段与 IPD 汇总](results/tables/04b_ipd_and_pan_external_synthesis.csv) · [芯片 Meta 汇总](results/tables/04_meta_analysis.csv)
- **混杂敏感性与评分尺度**：在 323 名精确匹配 DNA 推导肿瘤纯度（ABSOLUTE）的 TCGA 患者中，校正纯度后增殖 HR 为 **1.33 → 1.34**（Spearman ρ = +0.05），缺氧 HR 为 **1.34 → 1.54**（ρ = −0.29），免疫 HR 为 **0.98 → 1.26**（ρ = −0.70）。在各队列匹配病例上额外校正吸烟史或驱动基因（`EGFR`、`KRAS`、`TP53`、`STK11`）对增殖点估计影响有限。TCGA 均值 gene-wise z-score 的增殖 HR 经样本中位数中心化后由 **1.08** 回升至 **1.24**（tie-aware ssGSEA 为 **1.31**）。[纯度敏感性](results/tables/03c_purity_sensitivity.csv) · [吸烟及驱动基因](results/tables/03d_driver_smoking_adjustment.csv) · [评分尺度敏感性](results/tables/06_sensitivity.csv)
- **冗余性基准对比**：增殖复合评分与基于秩评分的 `CCP-31`（31 基因）和 `meta-PCNA`（131 基因）高度相关（Spearman ρ = **0.927–0.965**），与 `E2F`/`G2M` 相关为 **0.94–0.98**。在全部 5 个外部队列的 10 组配对冻结迁移对比中，复合评分相对秩评分 `CCP-31` 或 `meta-PCNA` 的 Harrell ΔC 置信区间均包含零。[关联基准](results/tables/03e_published_signature_benchmark.csv) · [迁移基准](results/tables/05d_published_signature_transport.csv)
- **冻结迁移区分度与局部基线重校准**：加入增殖复合评分在 CPTAC 的 Uno IPCW ΔC 为 **+0.071（95% CI −0.023 至 0.192）**，在 GSE72094 的 Harrell ΔC 与 Uno ΔC 分别为 **+0.034（0.008 至 0.062）** 与 **+0.051（0.006 至 0.102）**。在低基线死亡率的 GSE31210 中，复合模型 36 个月平均预测死亡率经固定线性预测子的局部基线重校准由 **25.6%** 降至 **8.2%**（观测 KM 为 **8.0%**），IPCW Brier 由 **0.0935** 降至 **0.0675**；重校准风险在用于估计基线的同一目标队列上评估。[CPTAC 迁移](results/tables/05b_rnaseq_transport.csv) · [GEO 迁移](results/tables/05_incremental_value.csv) · [GEO 重校准](results/tables/05c_recalibration.csv)

## 仓库结构

```text
.
├── LICENSE
├── Makefile
├── README.md
├── run_all.sh
├── config/                         # 队列配置、基因集说明与已发表签名列表
├── data/
│   ├── README.md                   # 数据来源、映射与校验说明
│   ├── checksums.csv               # 公开输入文件 MD5 校验清单
│   └── genesets/                   # MSigDB Hallmark GMT 定义
├── R/                              # 可复用核心函数库（不直接产生文件副作用）
├── scripts/                        # 按编号顺序执行的分析与绘图脚本（00–99）
├── docs/                           # 英文手稿、分析方案、方法学、溯源与局限性说明
└── results/
    ├── figures/
    │   ├── main/                   # 主图 Figure 1–5（PDF 与 300 dpi PNG）
    │   └── supplementary/          # 补充图 Figure S1–S12（PDF 与 300 dpi PNG）
    └── tables/                     # 统计结果 CSV 与 30 页签汇总工作簿
```

## 分析脚本与产出对照

| 脚本 | 主要任务 | 核心产出 |
|---|---|---|
| [install_packages.R](scripts/install_packages.R) | 安装缺失的 CRAN 与 Bioconductor 依赖包 | 本地 R 运行环境 |
| [fetch_inputs.R](scripts/fetch_inputs.R) | 下载公开数据集并核对 MD5 与文件大小 | `data/raw/`（本地忽略） |
| [00_setup.R](scripts/00_setup.R) | 检查运行依赖、输入快照与记录运行环境 | `00_params.csv`、`docs/sessionInfo.txt` |
| [01_load_data.R](scripts/01_load_data.R) | 队列清洗、ID 匹配、平台注释与单位统一 | `01_cohort_audit.csv`、`01_tcga_participant_flow.csv` |
| [02_score.R](scripts/02_score.R) | tie-aware ssGSEA 评分与共享迁移基因宇宙构建 | `02_axis_definition.csv`、`02_gene_coverage.csv` |
| [03_tcga_primary.R](scripts/03_tcga_primary.R) | TCGA 主关联、C-index 阶梯与分期分层分析 | `03_cox_*.csv`、`03_stage_stratified.csv` |
| [03b_mutual_adjustment.R](scripts/03b_mutual_adjustment.R) | 四程序联合互调 Cox 模型 | `03b_mutual_adjustment.csv` |
| [03c_benchmark.R](scripts/03c_benchmark.R) | 纯度、吸烟与驱动基因敏感性及已发表签名基准 | `03c_*.csv`、`03d_*.csv`、`03e_*.csv` |
| [04_external_meta.R](scripts/04_external_meta.R) | 外部队列关联估计与芯片两阶段 REML-mKH 合并 | `04_cohort_estimates.csv`、`04_meta_analysis.csv` |
| [05_incremental.R](scripts/05_incremental.R) | GEO 队列冻结模型迁移、Uno C 与基线重校准 | `05_incremental_value.csv`、`05c_recalibration.csv` |
| [05b_rnaseq_validation.R](scripts/05b_rnaseq_validation.R) | CPTAC 独立评估及五队列两阶段/一阶段 IPD 汇总 | `05b_rnaseq_*.csv`、`04b_ipd_and_pan_external_synthesis.csv` |
| [06_negative_control.R](scripts/06_negative_control.R) | 结构匹配随机基因集、标签置换与尺度敏感性 | `06_random_set_null.csv`、`06_sensitivity.csv` |
| [07_figures.R](scripts/07_figures.R) | 绘制主图（Figure 1–5）与补充图（Figure S1–S12） | `results/figures/main/`、`results/figures/supplementary/` |
| [08_tables.R](scripts/08_tables.R) | 汇总报告表与 30 页签 Excel 工作簿 | `table1`–`table9.csv`、`project_tables.xlsx` |
| [99_verify.R](scripts/99_verify.R) | 核验输入校验、数学不变量、输出表与图件完整性 | `results/logs/99_verify.log` |

## 运行环境与复现

所需 R 软件包及已测版本见 [运行环境](docs/sessionInfo.txt)，其中 `org.Hs.eg.db`、`hgu133plus2.db` 与 `AnnotationDbi` 来自 Bioconductor，PDF 导出需要支持 Cairo 的 R 环境。首次运行前可按需安装依赖，然后在仓库根目录执行数据下载与全流程脚本：

```sh
Rscript scripts/install_packages.R
Rscript scripts/fetch_inputs.R
Rscript scripts/run_pipeline.R
```

也可以通过 `make all` 或 `bash run_all.sh` 运行完整流程。`scripts/run_pipeline.R` 支持按数字前缀重跑部分步骤，例如基于本地中间对象重建图表、工作簿并运行校验：

```sh
Rscript scripts/run_pipeline.R 07 08 99
```

若仅运行 TCGA 与三套 GPL570 核心芯片队列（跳过需要额外下载的 GSE72094 与 CPTAC），可在运行前设置环境变量 `LUAD_CORE_ONLY=1`。

## 文档、图件与数据表索引

| 内容 | 文件 |
|---|---|
| 英文手稿（含主图/附图图注与声明） | [Markdown](docs/manuscript.md) · [Word 投稿版](docs/manuscript.docx) |
| 研究设计、实现细节与数据溯源 | [分析方案 v2](docs/analysis_plan_v2.md) · [方法](docs/methods.md) · [数据来源](docs/data_provenance.md) |
| Figure 1 · 队列概况与评分分布 | [PDF](results/figures/main/Figure1_Cohort_and_Scores.pdf) · [PNG](results/figures/main/Figure1_Cohort_and_Scores.png) |
| Figure 2 · 关联分析与 IPD 汇总 | [PDF](results/figures/main/Figure2_Association_and_IPD_Synthesis.pdf) · [PNG](results/figures/main/Figure2_Association_and_IPD_Synthesis.png) |
| Figure 3 · 特异性、纯度与驱动基因 | [PDF](results/figures/main/Figure3_Specificity_Purity_and_Drivers.pdf) · [PNG](results/figures/main/Figure3_Specificity_Purity_and_Drivers.png) |
| Figure 4 · 简单标志物与已发表签名基准 | [PDF](results/figures/main/Figure4_Redundancy_Benchmark.pdf) · [PNG](results/figures/main/Figure4_Redundancy_Benchmark.png) |
| Figure 5 · 冻结风险校准与基线重校准 | [PDF](results/figures/main/Figure5_Calibration_and_Recalibration.pdf) · [PNG](results/figures/main/Figure5_Calibration_and_Recalibration.png) |
| 补充图（Figure S1–S12） | [补充图目录](results/figures/supplementary/) · [补充图图注](docs/manuscript.md#supplementary-figure-legends) |
| 30 页签汇总工作簿与报告索引 | [工作簿](results/tables/project_tables.xlsx) · [报告索引](docs/reporting_crosswalk.md) |
| 解释边界与版本变更记录 | [局限性](docs/limitations.md) · [修订记录](docs/changelog.md) · [归档方案 v1](docs/analysis_plan.md) |

## 数据说明与许可

全部源数据均来自公开数据库（TCGA、NCBI GEO、LinkedOmics CPTAC），本地下载的原始文件与患者级衍生中间对象由 `.gitignore` 排除在版本控制之外，可通过 [输入校验清单](data/checksums.csv) 与上述复现命令完整重建。本研究代码基于 MIT 协议发布，见 [LICENSE](LICENSE)。
