# DeepLoc 2.0 protocol benchmark / DeepLoc 2.0 协议基准

## 中文

本基准在同一套 DeepLoc 2.0 数据和评测协议下比较两个模型：

- **DeepLoc 2.0**：官方 ProtT5 五折 checkpoint；
- **ProteinDigger Nucleus Specialist**：ESM2-3B layer-36 残基表示、label-query head 与三查询 Nucleus specialist。

共有 28,300 条有效蛋白。官方 partition 0–4 分别作为外层测试折；第 *i* 个 checkpoint 只预测第 *i* 折，不做模型 ensemble。每折阈值使用其余四折按 DeepLoc 2.0 原始 per-label maximum-MCC 算法确定。ProteinDigger checkpoint 只按 inner-validation focal loss 选择，外层折标签不参与模型选择。

| 指标 | DeepLoc 2.0 | ProteinDigger Nucleus Specialist | 差值 |
|---|---:|---:|---:|
| Overall micro-F1 | 0.7253 ± 0.0140 | **0.7462 ± 0.0133** | **+2.09 pp** |
| Overall macro-F1 | 0.6579 ± 0.0080 | **0.6647 ± 0.0094** | **+0.68 pp** |
| Nucleus precision | 0.8177 ± 0.0328 | **0.8383 ± 0.0204** | **+2.06 pp** |
| Nucleus F1 | 0.7894 ± 0.0135 | **0.7972 ± 0.0200** | **+0.78 pp** |
| Nucleus MCC | 0.6861 ± 0.0172 | **0.7008 ± 0.0244** | **+1.47 pp** |
| Nucleus ROC-AUC | 0.9281 ± 0.0065 | **0.9330 ± 0.0068** | **+0.49 pp** |
| Nucleus PR-AUC | 0.8843 ± 0.0155 | **0.8959 ± 0.0151** | **+1.16 pp** |

数值是五个官方外层折的平均值 ± 标准差。完整机器可读值、数据锁和源结果哈希见 [`metrics.json`](metrics.json)。该表对应官方五折 OOF 基准，不包含额外的 HPA 外部测试集。

## English

This benchmark compares two models under the same DeepLoc 2.0 dataset and evaluation protocol:

- **DeepLoc 2.0:** official five-fold ProtT5 checkpoints;
- **ProteinDigger Nucleus Specialist:** ESM2-3B layer-36 residue representations, a label-query head, and a three-query Nucleus specialist.

The logical dataset contains 28,300 proteins. Official partitions 0–4 serve as outer test folds; checkpoint *i* predicts partition *i* only, with no model ensemble. Thresholds are fitted on the other four partitions using the original DeepLoc 2.0 per-label maximum-MCC procedure. ProteinDigger checkpoints are selected only by inner-validation focal loss, without using outer-fold targets.

| Metric | DeepLoc 2.0 | ProteinDigger Nucleus Specialist | Difference |
|---|---:|---:|---:|
| Overall micro-F1 | 0.7253 ± 0.0140 | **0.7462 ± 0.0133** | **+2.09 pp** |
| Overall macro-F1 | 0.6579 ± 0.0080 | **0.6647 ± 0.0094** | **+0.68 pp** |
| Nucleus precision | 0.8177 ± 0.0328 | **0.8383 ± 0.0204** | **+2.06 pp** |
| Nucleus F1 | 0.7894 ± 0.0135 | **0.7972 ± 0.0200** | **+0.78 pp** |
| Nucleus MCC | 0.6861 ± 0.0172 | **0.7008 ± 0.0244** | **+1.47 pp** |
| Nucleus ROC-AUC | 0.9281 ± 0.0065 | **0.9330 ± 0.0068** | **+0.49 pp** |
| Nucleus PR-AUC | 0.8843 ± 0.0155 | **0.8959 ± 0.0151** | **+1.16 pp** |

Values are mean ± standard deviation across the five official outer folds. See [`metrics.json`](metrics.json) for full-precision values, the dataset lock, and source-result hashes. This table reports the official five-fold OOF benchmark and does not include the separate external HPA test set.
