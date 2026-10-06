#!/usr/bin/env python3
"""Write the concise scientific README for the clean Stage1-B1 release."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--release-root", type=Path, required=True)
    args = parser.parse_args()

    counts_frame = pd.read_csv(args.release_root / "03_audit/b1_counts.tsv", sep="\t")
    counts = dict(zip(counts_frame["metric"], counts_frame["value"]))
    formal = int(counts["corrected_B1_positive"])
    legacy = int(counts["legacy_B1_positive"])
    intersection = int(counts["corrected_and_legacy_positive"])
    corrected_only = int(counts["corrected_only_positive"])
    legacy_only = int(counts["legacy_only_positive"])
    rerun = int(counts["corrected_full_rerun_B1_positive"])
    rerun_diff = int(counts["corrected_formal_vs_rerun_set_difference"])
    flips = int(counts["corrected_recompute_nls_threshold_flips_14556"])
    deep_positive = int(counts["deeploc_gt_0p5"])
    corrected_nls_positive = int(counts["corrected_nls_gt_0p5"])
    legacy_nls_positive = int(counts["legacy_nls_gt_0p5"])

    text = f"""# Stage1-B1：DeepLoc–NLSExplorer 双模型核定位筛选

状态：`COMPLETE / CURRENT_STAGE0_UNIVERSE`
方法版本：`stage1_B1_complete_current_universe_20260904`

## 1. 目标与输入

B1 用两个互补的序列模型筛选没有适用于 canonical/displayed sequence 的严格 UniProt `Nucleus` 注释、且不属于当前 anchor 的蛋白：

- DeepLoc 2.0 评估整条蛋白的核定位倾向；
- NLSExplorer 评估序列窗口中的核输入信号证据。

唯一输入来自 Stage0：

```text
../00_common_input/02_non_anchor/
route_B_uniprot_without_canonical_nucleus_non_anchor_14625.tsv
route_B_uniprot_without_canonical_nucleus_non_anchor_14625.fasta
```

该表包含 **14,625 条 reviewed human canonical proteins**。B1 对这 14,625 条蛋白使用同一模型、同一阈值统一判定；不存在“保留某个旧候选集合再追加新候选”的集合并联。

## 2. 模型与推理口径

### 2.1 DeepLoc 2.0 Accurate

- 编码器：ProtT5-XL-UniRef50；
- 分类器：5 个官方 DeepLoc 2.0 ProtT5 heads 的概率均值；
- B1 使用 `Nucleus` 概率；
- 长度超过 4,000 aa 时沿用既定的中间截断口径；
- 当前 14,625 条均具有连续 DeepLoc 分数。

### 2.2 NLSExplorer 的两个保留版本

两个版本使用相同的 ESM-1b 编码器、NLSExplorer 权重、1,022 aa 窗口、200 aa 重叠和 max-window 蛋白汇总规则。差别只发生在 NLSExplorer 分类头的 batching。

| 版本 | 分类头输入 | 本项目地位 |
|---|---|---|
| `CORRECTED_EXACT_LENGTH` | 仅将真实长度完全相同的窗口组成 batch，因此分类头不含人为补零 | **正式 B1** |
| `LEGACY_BATCH_PADDED` | 不同长度窗口补零后直接送入无 padding mask 的分类头 | 方法敏感性对照 |

legacy 版本的分数会受 batch 组成影响，因此保留其完整分数和候选结果用于审计，但不与 corrected 结果取并集，也不向正式 B1 补入 legacy-only 蛋白。

## 3. 正式筛选规则

```text
DeepLoc_positive = deeploc_nucleus_prob > 0.5
NLS_positive     = nlsexplorer_nls_prob_max > 0.5
B1_positive      = DeepLoc_positive AND NLS_positive
```

阈值均为严格大于 0.5。候选按 `joint_min_prob`、`joint_mean_prob` 降序和 accession 升序稳定排序；排序不改变阳性定义。

正式 corrected 连续分数表由两部分构成：已完成一致性验证的 14,556 条 corrected 分数，以及对 Stage0 新纳入 69 条蛋白使用完全相同方法得到的 corrected 分数。合并后先验证为完整、无重复的 14,625 条连续分数，再统一应用上述双模型阈值。

## 4. 结果

### 4.1 正式 corrected 结果

| 层级 | proteins |
|---|---:|
| B1 输入全集 | 14,625 |
| DeepLoc `Nucleus > 0.5` | {deep_positive:,} |
| corrected NLSExplorer `max-window > 0.5` | {corrected_nls_positive:,} |
| 两模型均 `> 0.5`，正式 B1 | **{formal:,}** |

### 4.2 legacy 方法对照

| 比较项 | proteins |
|---|---:|
| legacy NLSExplorer `max-window > 0.5` | {legacy_nls_positive:,} |
| legacy 双模型阳性 | {legacy:,} |
| corrected 与 legacy 共同阳性 | {intersection:,} |
| corrected-only | {corrected_only:,} |
| legacy-only（不进入正式 B1） | {legacy_only:,} |

### 4.3 corrected 全量复算审计

为验证从输入到输出的完整可复现性，14,625 条 corrected NLSExplorer 另行全量复算一次。该复算产生 {rerun:,} 条双模型阳性，与正式集合的对称差为 {rerun_diff:,} 条；旧 14,556 条在 NLS `>0.5` 边界上的翻转为 {flips:,} 条。细微分差来自 GPU 半精度复算，正式数值仍按预先冻结的“已验证 14,556 + 同法新增 69”接口生成。

## 5. 数据接口

正式 B1：

```text
01_model_scores/b1_corrected_formal_scores_14625.tsv
02_candidates/b1_corrected_positive_both_gt0p5.tsv
02_candidates/b1_corrected_positive_both_gt0p5.accessions.txt
02_candidates/b1_corrected_positive_both_gt0p5.fasta
```

legacy 对照：

```text
01_model_scores/b1_legacy_batch_padded_scores_14625.tsv
02_candidates/legacy_batch_padded_comparator/b1_legacy_positive_both_gt0p5.tsv
02_candidates/legacy_batch_padded_comparator/b1_legacy_positive_both_gt0p5.accessions.txt
02_candidates/legacy_batch_padded_comparator/b1_legacy_positive_both_gt0p5.fasta
```

审计：

```text
03_audit/b1_counts.tsv
03_audit/b1_version_comparison.tsv
03_audit/b1_candidate_membership_comparison.tsv
03_audit/b1_corrected_recompute_consistency_14556.tsv
03_audit/b1_corrected_full_rerun_scores_14625.tsv
03_audit/b1_manifest.tsv
```

原始推理日志及窗口级分数保存在 `04_inference_runs/`。`b1_manifest.tsv` 记录输入、模型权重、结果文件的 SHA256。

## 6. 下游边界

- B2 必须从本版正式 corrected B1 中排除 direct-access 蛋白后重新计算；旧 B2 不再是当前 Stage0/B1 的正式下游。
- Stage2-B1 必须针对本版正式 B1 重新建立 HI-union / RF2-PPI 结果；历史 Stage2 与“旧候选 + 新增批次”的接口仅保留在归档中。
- legacy-only 蛋白不进入 B2 排除集合，也不进入正式 Stage2-B1。

## 7. 复现

```bash
python3 scripts/build_b1_14625.py [冻结输入参数]
python3 scripts/validate_b1_14625.py [冻结输入参数]
```

验证覆盖 14,625 = 14,556 + 69 的不重叠并集、两版连续分数覆盖、严格阈值、候选排序、accession/FASTA 接口、corrected 值来源、manifest SHA256 和未完成 `.partial` 文件。
"""

    path = args.release_root / "README.md"
    partial = path.with_name(path.name + ".partial")
    if path.exists() or partial.exists():
        raise FileExistsError(f"Refusing to overwrite README: {path}")
    partial.write_text(text)
    os.replace(partial, path)
    print(path)


if __name__ == "__main__":
    main()
