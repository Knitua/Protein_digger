# Current result release

## Formal routes

- `candidates_1580.tsv`: one compact record per formal route-level candidate gene. Columns are route, candidate gene, UniProt accession and primary protein name.
- `candidate_anchor_pairs_5957.tsv`: one row per retained candidate–anchor protein pair, with route, candidate identifiers, regulatory-anchor identifiers and evidence class.

| Route | Candidate records | Pair rows |
|---|---:|---:|
| A | 954 | 3,530 |
| B1 | 429 | 1,593 |
| B2 | 85 | 344 |
| B3 | 112 | 490 |
| **Total** | **1,580** | **5,957** |

B3 pair rows collectively retain 152 alternative isoform accessions for 112 genes. Route-level counts are therefore not interchangeable with unique accession counts.

## Stage 2.5 v3

`stage25/` contains the nine extension candidates, their passing seed edges, the machine-readable core summary and validation receipt, the positive-control threshold matrix and the post-freeze RAS case audit.

Stage 2.5 results are a separate extension set and are not added to the 1,580 formal records or 5,957 candidate–anchor pairs.
