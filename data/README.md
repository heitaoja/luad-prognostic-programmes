# data/

本目录只放数据与基因集，不放代码与结果。

```
genesets/     Hallmark GMT 缓存（离线可用）
checksums.csv 原始文件的 md5（模板，下载后由 01_load_data.R 追加）
raw/          下载的原始文件（不入库，见 .gitignore）
processed/    中间对象（.rds，不入库，可由 01/02 重建）
```

## 已纳入的数据

| 数据集 | 来源 | 用途 | 规模 |
|---|---|---|---|
| TCGA-LUAD | GDC Data Portal（经 gdc-hub 镜像） | 发现队列 | 515 例有 OS / 185 例死亡 |
| GSE31210 | GEO | 验证队列 | 204 例 / 30 例死亡 |
| GSE37745 | GEO | 验证队列（仅腺癌） | 106 例 / 77 例死亡 |
| GSE50081 | GEO | 验证队列（仅腺癌） | 127 例 / 51 例死亡 |
| Hallmark 基因集 | MSigDB（Enrichr 分发） | 程序定义 | 50 个基因集 |

## 规则

1. `raw/` 只读，任何清洗写到 `processed/`，并记录理由。
2. 下载后立刻算 md5 写进 `checksums.csv`，重下时比对。
3. 大文件不入 git（见 `.gitignore`）。
4. 每个数据集的来源、字段、取舍见 [`../docs/data_provenance.md`](../docs/data_provenance.md)。

公开数据不需要伦理审批，但方法里写清来源与使用条款。
