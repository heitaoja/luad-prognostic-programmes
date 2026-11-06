# 数据来源与取舍

## 1. 纳入的数据

### TCGA-LUAD（发现队列）

| 项目 | 内容 |
|---|---|
| 来源 | GDC Data Portal；通过 UCSC Xena 维护的 `gdc-hub` 镜像下载 |
| 文件 | `TCGA-LUAD.star_counts.tsv.gz`（64 MB）、`TCGA-LUAD.survival.tsv.gz`、`TCGA-LUAD.clinical.tsv.gz` |
| 表达 | 60,660 个 Ensembl ID × 589 个样本；值域 0–22.4，为 `log2(count+1)` |
| 样本筛选 | 只保留 barcode `-01`（原发肿瘤）：528 例；排除 `-02`（复发，2 例）与 `-11`（癌旁，58 例） |
| 基因筛选 | 映射到 36,650 个基因符号 |
| 生存 | `OS.time` 原始单位为天（中位 669 天），转月后中位 22.0 月；185 例死亡 |
| 临床 | AJCC 分期（I 280 / II 122 / III 80 / IV 25 / 缺失 8）、年龄（缺失 10）、性别、吸烟包年 |
| 校验 | `data/checksums.csv` |

### GSE31210（验证队列）

| 项目 | 内容 |
|---|---|
| 平台 | Affymetrix HG-U133 Plus 2.0（GPL570） |
| 原始样本 | 246（226 例肿瘤 + 20 例正常肺） |
| 筛选 | `tissue == primary lung tumor` 且有 OS，且未被作者标记 `exclude for prognosis analysis` |
| 最终 | 204 例，30 例死亡，中位 OS 59.6 月 |
| 尺度 | 矩阵为**线性强度**，已 `log2(x+1)` |
| 额外字段 | 分期、年龄、性别、吸烟、EGFR/KRAS/ALK 状态、无复发生存 |

### GSE37745（验证队列）

| 项目 | 内容 |
|---|---|
| 平台 | GPL570 |
| 原始样本 | 196 例 NSCLC（混合组织学） |
| 筛选 | `histology == adeno` 且有 OS → **106 例**，77 例死亡，中位 OS 47.8 月 |
| 尺度 | 已是 log2 |

### GSE50081（验证队列）

| 项目 | 内容 |
|---|---|
| 平台 | GPL570 |
| 原始样本 | 181 例 NSCLC |
| 筛选 | `histology == adenocarcinoma` 且有 OS → **127 例**，51 例死亡，中位 OS 52.4 月 |
| 尺度 | 已是 log2；**`survival time` 原始单位为年**，已 ×12 |

### 基因集来源

Hallmark（50 个基因集）来自 Enrichr 的 `MSigDB_Hallmark_2020` 库，
缓存于 `data/genesets/Hallmark_Enrichr.gmt`（离线可用，不依赖在线接口）。

## 2. 排除的数据及理由

| 数据 | 排除理由 |
|---|---|
| **GSE72094**（442 例 LUAD，vital_status + survival_time_in_days） | 平台是 Rosetta/Merck HuRSTA 自定义芯片（GPL15048），探针名形如 `merck-AA001022_at`；Bioconductor 无对应注释包（`huRSTA2a520709.db` 已从当前版本移除），GEO 的 `/annot/` 目录也不可访问。无法把探针映射到基因符号，**且不做强行映射**（会引入错误注释）。这是本可以最大的验证队列，值得在局限里说明。 |
| GSE65858（270 例 HNSCC） | 癌种不符（头颈鳞癌） |
| GSE41613（97 例 OSCC） | 癌种不符 |
| GSE68465 / GSE30219 / GSE14520 | 与本研究问题或癌种不符，或组织学混合度过高 |

## 3. 已知的取舍

1. **组织学严格限定腺癌**：代价是样本量下降（GSE37745 196→106，GSE50081 181→127），
   好处是四个队列的生物学可比。这个取舍会影响外部效度：结论只对肺腺癌成立。
2. **不合并表达矩阵**：只合并效应量。跨平台矩阵拼接会引入批次效应，且难以在审稿中辩护。
3. **不做生存时间插补**：没有 OS 的样本直接排除，不推断。

## 4. 复现方式

```bash
# TCGA（约 65 MB）
curl -sSL -o data/raw/tcga/TCGA-LUAD.star_counts.tsv.gz \
  https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.star_counts.tsv.gz
curl -sSL -o data/raw/tcga/TCGA-LUAD.survival.tsv.gz \
  https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.survival.tsv.gz
curl -sSL -o data/raw/tcga/TCGA-LUAD.clinical.tsv.gz \
  https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LUAD.clinical.tsv.gz

# GEO（每个 20–60 MB）
for g in GSE31210 GSE37745 GSE50081; do
  pre=$(echo $g | sed 's/...$//')
  curl -sSL -o data/raw/geo/${g}_series_matrix.txt.gz \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/${pre}nnn/${g}/matrix/${g}_series_matrix.txt.gz"
done

# 基因集
curl -sSL -o data/genesets/Hallmark_Enrichr.gmt \
  "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020"
```
