# Stage 2.5 v3 protocol

Stage 2.5 v3 uses the current A1/B1/B2 Stage 2 universe while keeping the previously frozen sequence, domain and direct-function rules.

## Input lock

| Route | Stage 2 screened | final80 positive seeds | final80 failure pool |
|---|---:|---:|---:|
| A1 | 2,877 | 1,348 | 1,529 |
| B1 | 1,350 | 429 | 921 |
| B2 | 266 | 85 | 181 |
| **Total** | **4,493** | **1,862** | **2,631** |

## Frozen selection rules

1. Run a general k-mer/length prefilter followed by BLOSUM62 local alignment.
2. Require identity ≥70%, seed coverage ≥80% and candidate coverage ≥80%.
3. Use Pfam first and InterPro as fallback; require identical domain order and copy number plus normalized boundary minimum IoU ≥0.80.
4. Accept GOA BP/MF direct annotations only for EXP, IDA, IMP and IGI evidence; reject `NOT`, `contributes_to` and `colocalizes_with` qualifiers.
5. Accept only direct human, non-disease, non-inferred Reactome catalyst/regulator roles.
6. Require at least one of:
   - Reactome direct Jaccard ≥0.50 with at least one shared event;
   - GO-MF direct Jaccard ≥0.50 with at least two shared terms;
   - GO BP+MF direct IDF-Jaccard ≥0.25 with at least three shared terms.
7. Freeze the general result before positive-control sensitivity analysis and the separate RAS case audit.

## Validation contract

The selected candidates must be unique primary Stage 2 failures, all passing edges must satisfy the sequence, domain-architecture and direct-function rules, and no case-specific logic may participate in core selection. The machine-readable validation receipt is published with the results.
