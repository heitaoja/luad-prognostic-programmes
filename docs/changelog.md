# 变更记录（偏离预注册与设计转向的记录）

> 规则：任何偏离 `docs/analysis_plan.md` 的改动都要在这里记一行。
> 记"改了什么、为什么、改之前的结论是什么"——尤其是最后一栏。

| 日期时间 | 版本 | 改动 | 原因 | 改之前的结论 |
|---|---|---|---|---|
| 2026-11-06 17:40 | v0.1 | 立项：原设计为"树突状细胞成熟轴"（沿用已有的单细胞数据） | 先尝试从单细胞延长到 bulk | 无 |
| 2026-11-06 17:52 | v0.2 | **放弃单细胞驱动的 DC 轴**，改为纯 bulk 设计 | 单细胞发现层会拖长周期，且与 bulk 技能不直接对应 | — |
| 2026-11-06 17:56 | v0.3 | 原假设"增殖−免疫平衡轴"（`z(P) − z(I)`）在 TCGA 探索中不成立：免疫成分单独看 HR 0.99（p = 0.91），P×I 交互 p = 0.92 | 数据不支持，拒绝硬写 | 平衡轴 adj HR 1.21（1.04–1.42），但无可解释的交互 |
| 2026-11-06 17:58 | v0.4 | 改为**四个程序的正面对比**（增殖、缺氧、EMT/间质、免疫） | 单轴设计在免疫阴性时没有故事；四程序对比是更完整、更可证伪的问题 | 平衡轴 |
| 2026-11-06 18:05 | v0.9 | 外部队列确认为 GSE31210 / GSE37745 / GSE50081；**排除 GSE72094** | GSE72094 为 HuRSTA 自定义芯片，无可用注释包，探针无法映射 | 原计划 4 个验证队列 |
| 2026-11-06 18:20 | v1.0 | 分析计划定稿（`docs/analysis_plan.md`），外部队列在读入前已确定模型与协变量 | 预注册 | — |

## 约定

- 程序定义（`AXES`）定稿后不再改动；如需改动，版本号 +1 并说明。
- 主终点变更必须明示。
- 队列剔除写进 `docs/data_provenance.md`，不只写在代码注释里。
- 探索性分析单独标注，不进主结论。

## 已知的"事后"决策（诚实披露）

1. **四个程序的选择**：候选池是 Hallmark 全集，但具体选这四类，是在 TCGA 里做过一次
   Hallmark 扫描后才定的。TCGA 的效应量可能受选择偏倚影响，因此正文以外部合并估计为主要证据，
   并把 TCGA 用于发现和模型训练。
2. **分期分层**：是在看到总效应后才做的（探索性），已在 README 中标注。
3. **互调模型**（四程序同时入模）：在看到单变量结果后补充，用于回应"纯度/共线"质疑，标注为探索性。

## 2026-11-24 revision

- Corrected the GMT description parsing, used midranks for expression ties and recalculated all 14 gene-set scores.
- Changed the TCGA unit from specimen to participant and tested an alternate vial rule. Separated median OS and potential follow-up.
- Restored GSE72094 via the versioned original-centre HuRSTA map. Kept the three original external cohorts identifiable and added a four-cohort synthesis.
- Added modified Knapp–Hartung inference, Holm correction, model diagnostics, molecular benchmark comparisons, and exact-specimen purity adjustment.
- Moved patient-level derived scores to ignored processed data. Replaced the previous within-cohort incremental analysis and risk reclassification summaries with frozen cohort transport, calibration and paired C-index comparisons.
- Used a shared assay gene list and within-sample averaging of duplicate features for the final transport analysis. This removes target-cohort mean expression from duplicate-probe preprocessing.
- Added an independent RNA-seq assessment using pinned CPTAC-LUAD inputs, including separate association estimates and evaluation of TCGA models fitted to the 17,593 comparable genes.

## 2026-11-24 完整发布

以下扩展在初版修订分析审查后完成。历史发现和修订设计的时间顺序见分析方案。[修订分析方案](analysis_plan_v2.md)

| 步骤 | 改动 | 原因和原有证据 |
|---|---|---|
| Step 1 | 在相同完整病例上比较吸烟及驱动突变校正，增加驱动基因亚组和样本中心化 z-score | 原分析以分期、年龄、性别为主，TCGA 未中心化的 gene-wise z-score 增殖 HR 为 1.08 |
| Step 2 | 增加分层 IPD 和纳入 CPTAC 的五队列合并，保留原数组合并表 | 原三队列与四队列 REML-mKH 提供芯片证据，CPTAC 原先单独报告 |
| Step 3 | 增加 Uno IPCW C 和配对间隔，增加预测子固定的 36 个月基线重校准 | 原冻结迁移以 Harrell C 与完全冻结的绝对风险为主，GSE31210 存在明显高估 |
| Step 4 | 加入 Wistuba 肺癌 CCP-31 和 Venet meta-PCNA 的关联、条件模型及固定迁移对照 | 原冗余性基准仅包含 MKI67、E2F 和 G2/M |
| Step 5 | 导出五张组合主图及 30 页签工作簿，保留原单图 | 使新旧分析结果在主图、附图和工作表中可追溯 |
| Step 6 | 扩写英文手稿和配套文档，增加稿件规模、结构、图注与参考文献检查 | 汇总完整证据，并区分发现、扩展分析、冻结外部验证与目标队列内重校准 |

[英文稿件](manuscript.md) · [图表索引](reporting_crosswalk.md) · [结果工作簿](../results/tables/project_tables.xlsx)
