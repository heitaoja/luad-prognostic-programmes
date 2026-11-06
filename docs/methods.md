# 方法与参数理由

> 每个参数都要写"为什么是这个值"。报告与复现说明的 Methods 从这里取。

## 1. 数据获取

| 数据 | 方式 | 说明 |
|---|---|---|
| TCGA-LUAD 表达 | `gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_counts.tsv.gz` | UCSC Xena 维护的 GDC 镜像；内容是 GDC STAR counts，已 `log2(count+1)` |
| TCGA-LUAD 生存 | `TCGA-LUAD.survival.tsv.gz` | `OS.time` 原始单位是**天**，已除 30.4375 转为月 |
| TCGA-LUAD 临床 | `TCGA-LUAD.clinical.tsv.gz` | AJCC 分期、年龄、性别、吸烟包年 |
| GEO 队列 | `ftp.ncbi.nlm.nih.gov/geo/series/.../matrix/GSE*_series_matrix.txt.gz` | 矩阵内含表达值与样本注释 |

**为什么不用 Xena 的 Toil 重组矩阵**：直接取 GDC 的 STAR counts，与临床字段同源、版本可追溯。

## 2. 基因注释

| 平台 | 注释包 | 说明 |
|---|---|---|
| TCGA Ensembl ID | `org.Hs.eg.db` | 去掉版本号后缀（`.15`）后映射到 Symbol |
| GPL570（三个芯片队列） | `hgu133plus2.db` | 探针 → Symbol |

**重复基因**：保留平均表达最高的一行，而不是取平均，避免把不同探针的信号混在一起。

## 3. 表达尺度

- TCGA：STAR counts → `log2(count+1)`（数据源已完成）。
- 芯片：GSE31210 的矩阵是**线性强度**（中位数 99.6），先 `log2(x+1)`；GSE37745 / GSE50081 已是 log2。
- 所有分析只在队列内标准化；**跨队列不比较绝对值**，只合并效应量。

## 4. ssGSEA 打分

**算法**：Barbie et al. 2009 的 single-sample GSEA。样本内按表达量降序排名，
命中基因权重 `(N − rank + 1)^0.25`，未命中基因步长 `−1/(N − k)`，富集分数 = 游走曲线积分。

**实现**：`R/lib_score.R` 用闭式求和代替逐位置游走：

```
ES = Σ_hits w_j·(N−j+1) − (Σ_nonhits (N−j+1)) / (N−k)
```

复杂度从 O(N) 降到 O(k)（k = 基因集大小），这是能在 500 次随机基因集上跑零分布的前提。
**正确性**：与逐步游走实现逐元素比较，最大绝对误差 9.7×10⁻¹⁷。

**为什么 α = 0.25**：GSVA 默认值，对基因排序不敏感。

**程序得分**：成员基因集 z 分数平均后再 z 标准化。用平均而不是取第一主成分，
是因为 PC1 的符号不稳定，跨队列不可比。

## 5. 生存模型

- Cox 比例风险模型；连续变量，报告每 1 SD 的 HR。
- 协变量：AJCC 分期（I/II/III/IV 因子）、年龄（标准化）、性别。
- 完整病例分析（TCGA 515 → 497，缺分期或年龄）；缺失情况见 `results/tables/01_cohort_audit.csv`。
- 事件/变量：TCGA 185 / 5 ≈ 37，远高于 EPV ≥ 10 的经验下限。
- KM 图仅用于可视化，不做推断（推断用连续变量模型）。

## 6. 合并

DerSimonian–Laird 随机效应模型，在 log-HR 尺度合并；报告合并 HR、95% CI、I²、τ²、Q。
选随机效应而非固定效应，是因为三个验证队列的人群与随访强度不同（I² 在 0–70% 之间）。

## 7. 增量价值

| 指标 | 实现 | 说明 |
|---|---|---|
| ΔC-index | 两个 Cox 模型的 `concordance` 之差；bootstrap 200 次给 CI | 在同一批完整病例上比较 |
| 似然比检验 | `anova(fit0, fit1)` | 嵌套模型 |
| 连续 NRI | 36 个月预测风险的变化方向，事件与非事件分别计算后相加 | — |
| DCA | IPCW 加权净获益，36 个月，阈值 0.05–0.60 | 删失用 KM 估计 |

## 8. 阴性对照

- 随机基因集：按真实基因集大小逐个匹配，500 次（TCGA）/200 次（GSE31210）。
- 标签置换：1000 次。
- 经验 p 用 `(计数+1)/(B+1)`，避免 p = 0。

## 9. 软件

R 4.6.2；`survival`、`data.table`、`ggplot2`、`patchwork`、`openxlsx`、
`org.Hs.eg.db`、`hgu133plus2.db`、`AnnotationDbi`。环境写入 `docs/sessionInfo.txt`。

## 10. 明确的未做项

- **肿瘤纯度未校正**（ESTIMATE / ABSOLUTE 未运行）。这是最大的方法学缺口，已写进局限。
- **未训练数据驱动模型**（LASSO / 随机森林）。刻意的设计选择：避免"用结局筛基因"的循环论证；
  代价是失去了与常规预后模型的直接比较。
- **未做竞争风险模型**：把非肺癌死亡视为删失，可能高估风险。
