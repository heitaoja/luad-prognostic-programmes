# 肺腺癌 Hallmark 转录程序的跨平台预后基准研究

本项目评估四类 Hallmark 转录程序在肺腺癌中的总生存关联、对临床和测量因素的敏感性，以及冻结的 TCGA 风险模型在外部队列的区分度与校准。研究纳入 6 个队列、1,443 名患者，其中 1,420 名患者和 469 例死亡进入临床完整病例分析。TCGA 用于发现和训练，四个 GEO 芯片队列与 CPTAC RNA-seq 队列提供外部证据。

## 核心结果

- 五队列两阶段 REML-mKH 合并中，增殖 HR 为 **1.35（95% CI 1.15–1.59）**，Holm 校正后 *P* = **0.027**。五队列一阶段分层 IPD 中，增殖和缺氧 HR 分别为 **1.35** 和 **1.40**，两者的 Holm 校正后 *P* 均为 **4.4 × 10⁻⁷**。[两阶段与 IPD 汇总](results/tables/04b_ipd_and_pan_external_synthesis.csv)
- 在同一批 323 名 TCGA 患者中，ABSOLUTE 纯度校正使增殖 HR 从 **1.33** 变为 **1.34**，缺氧 HR 从 **1.34** 变为 **1.54**。可用的吸烟与驱动基因校正对增殖点估计影响较小。TCGA 均值 gene-wise z-score 的增殖 HR 经样本中位数中心化由 **1.08** 回升至 **1.24**。[纯度](results/tables/03c_purity_sensitivity.csv) · [吸烟及驱动基因](results/tables/03d_driver_smoking_adjustment.csv) · [评分尺度](results/tables/06_sensitivity.csv)
- 增殖复合评分与秩评分版本的 CCP-31 和 meta-PCNA 高度相关，Spearman ρ 为 **0.927–0.965**。外部队列中，复合评分相对于两种已发表秩评分的配对 Harrell ΔC 置信区间均包含零。[关联基准](results/tables/03e_published_signature_benchmark.csv) · [迁移基准](results/tables/05d_published_signature_transport.csv)
- CPTAC 的复合评分相对临床模型的 Uno ΔC 为 **+0.071（95% CI −0.023 至 0.192）**。GSE31210 的复合模型 36 个月平均预测死亡率经基线重校准由 **25.6%** 降至 **8.2%**，观测为 **8.0%**，IPCW Brier 由 **0.0935** 降至 **0.0675**。更新后的风险在用于估计基线的同一目标队列上评估。[CPTAC 迁移](results/tables/05b_rnaseq_transport.csv) · [GEO 迁移](results/tables/05_incremental_value.csv)

## 手稿、图件与数据表

| 内容 | 文件 |
|---|---|
| 英文手稿、主图和补充图图注 | [手稿](docs/manuscript.md) |
| 研究设计、实现细节与数据来源 | [分析方案](docs/analysis_plan_v2.md) · [方法](docs/methods.md) · [数据来源](docs/data_provenance.md) |
| Figure 1 · 队列与评分分布 | [PDF](results/figures/Figure1_Cohort_and_Scores.pdf) · [PNG](results/figures/Figure1_Cohort_and_Scores.png) |
| Figure 2 · 关联与 IPD 汇总 | [PDF](results/figures/Figure2_Association_and_IPD_Synthesis.pdf) · [PNG](results/figures/Figure2_Association_and_IPD_Synthesis.png) |
| Figure 3 · 特异性、纯度和驱动基因 | [PDF](results/figures/Figure3_Specificity_Purity_and_Drivers.pdf) · [PNG](results/figures/Figure3_Specificity_Purity_and_Drivers.png) |
| Figure 4 · 简单与已发表签名基准 | [PDF](results/figures/Figure4_Redundancy_Benchmark.pdf) · [PNG](results/figures/Figure4_Redundancy_Benchmark.png) |
| Figure 5 · 冻结风险与基线重校准 | [PDF](results/figures/Figure5_Calibration_and_Recalibration.pdf) · [PNG](results/figures/Figure5_Calibration_and_Recalibration.png) |
| 完整 30 页签工作簿与图表索引 | [工作簿](results/tables/project_tables.xlsx) · [报告索引](docs/reporting_crosswalk.md) |
| 解释边界与版本记录 | [局限性](docs/limitations.md) · [变更记录](docs/changelog.md) |

## 复现

在仓库根目录运行以下命令。所需软件包及已用版本见运行环境记录。org.Hs.eg.db、hgu133plus2.db 和 AnnotationDbi 来自 Bioconductor。PDF 输出需要支持 Cairo 的 R。[运行环境](docs/sessionInfo.txt) · [输入校验清单](data/checksums.csv)

```sh
Rscript R/fetch_inputs.R
Rscript R/run_pipeline.R
```

完整流程运行 00 至 99 共 13 个脚本，并核查输入版本、患者流、统计输出、PDF/PNG 图件、Excel 页签以及手稿结构。指定前缀可重跑相关步骤，例如 `Rscript R/run_pipeline.R 07 08 99` 可用已有分析结果重建图表和执行检查。患者级衍生表与中间对象保存在本地忽略目录。原始设计与后续扩展的时间顺序见分析方案和变更记录。[原分析方案](docs/analysis_plan.md) · [修订分析方案](docs/analysis_plan_v2.md)
