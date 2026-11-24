# 肺腺癌转录程序的预后关联与跨队列模型验证

项目比较肺腺癌中四类 Hallmark 转录程序与总生存的关系，并评估固定的 TCGA 模型在外部队列中的表现。电子稿、分析代码、可追溯的输入清单、汇总表和可发表图件随项目提供。

## 本版结果

- **6 个队列，1,443 名患者。** 主分析包括 TCGA 和四个 GEO 芯片队列，共 1,338 人。四个 GEO 队列有 830 人、269 例死亡进入多因素合并分析。独立 CPTAC RNA-seq 队列有 105 人、23 例死亡进入补充验证。
- **增殖方向一致。** 原三队列 REML + 改良 Knapp–Hartung 合并 HR 为 **1.32（0.93–1.86）**。四队列扩展结果为 **1.34（1.10–1.62）**，四项 Holm 校正后 *p* = **0.069**。
- **DNA 纯度校正后的估计相近。** 在同一批 323 名 TCGA 患者中，增殖 HR 由 **1.33** 变为 **1.34**。
- **固定模型增量有限。** 加入复合增殖评分后，外部 C-index 增加 0.017–0.054，四组中一组区间完全高于零。复合评分对 MKI67、E2F 或 G2/M 的直接优势尚不明确，部分队列的绝对风险校准偏差较大。
- **CPTAC 中的估计精度有限。** 增殖 HR 为 **1.575（1.002–2.476）**。用 CPTAC 可测的 17,593 个基因在 TCGA 中重新拟合后，复合评分的外部 C-index 增量为 **0.010（−0.071 至 0.104）**。

| 内容 | 位置 |
|---|---|
| 英文稿件及图注 | [docs/manuscript.md](docs/manuscript.md) |
| 设计、方法与输入 | [分析方案](docs/analysis_plan_v2.md) · [方法细节](docs/methods.md) · [数据来源](docs/data_provenance.md) |
| 发表用图件 | [四程序森林图](results/figures/fig3_forest_all_cohorts.pdf) · [固定模型 ΔC](results/figures/fig6_incremental.pdf) · [简单指标对照](results/figures/fig9_marker_benchmark.pdf) · [校准](results/figures/fig10_calibration.pdf) · [纯度](results/figures/fig11_purity_sensitivity.pdf) |
| 汇总工作簿 | [results/tables/project_tables.xlsx](results/tables/project_tables.xlsx) |
| 独立 RNA-seq 验证 | [补充图 S1](results/figures/figS1_rnaseq_validation.pdf) · [纳入及测量覆盖](results/tables/05b_rnaseq_audit.csv) · [固定模型表现](results/tables/05b_rnaseq_transport.csv) |
| 输入版本及检查 | [data/checksums.csv](data/checksums.csv) · [运行环境](docs/sessionInfo.txt) |

## 复现

在仓库根目录使用 R，所需包及本版版本见 [运行环境](docs/sessionInfo.txt)。`org.Hs.eg.db`、`hgu133plus2.db` 和 `AnnotationDbi` 来自 Bioconductor，其余依赖来自 CRAN 或 R 自带库。PDF 图件需要 Cairo。

```sh
Rscript R/fetch_inputs.R
Rscript R/run_pipeline.R
```

已有原始文件可按 [清单](data/checksums.csv) 放入相应路径。完整流程校验文件版本，从数据清洗运行到成果核对。可用 `Rscript R/run_pipeline.R 04 05` 复算指定的关联合并和模型评估。原始数据、患者级衍生表和中间对象写入本地数据目录。历史 v1 设计见 [原分析方案](docs/analysis_plan.md) 与 [变更记录](docs/changelog.md)。
