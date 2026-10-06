# Source layout

The source tree mirrors the current analysis stages.

| Directory | Purpose |
|---|---|
| `anchors/` | Validate the merged regulatory-anchor collections |
| `stage0/` | Build and validate the isoform-aware canonical A/B partition |
| `stage1/route_a/` | Validate the annotation-first A route |
| `stage1/route_b1/` | Build the DeepLoc–NLSExplorer consensus route |
| `stage1/route_b2/` | Build network features and rank the signal-to-nucleus route |
| `stage1/route_b3/` | Prepare alternative isoforms, run nuclear models and construct the B3 union |
| `stage2/route_a/`, `route_b1/`, `route_b2/` | Join canonical candidates to HI-union and RF2-PPI evidence |
| `stage2/route_b3/` | Run PPLM-PPI scoring for isoform–anchor shards |
| `stage25/` | Stage 2.5 v3 sequence/domain/function extension and validation |

Most build scripts expose their data paths through `argparse`. A small number of frozen orchestration/manifest scripts preserve the directory constants used in the validated run; adapt those constants or call their component scripts with explicit paths when reconstructing the workflow elsewhere.

The scripts are designed to fail closed on unexpected row counts, duplicate accessions, incomplete model coverage, threshold mismatches or broken subset relations. Large upstream resources and model weights are not distributed in this repository.
