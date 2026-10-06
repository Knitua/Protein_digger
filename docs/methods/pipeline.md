# ProteinDigger current pipeline

## 1. Reference universe and anchors

The canonical reference universe contains 20,416 reviewed human UniProt proteins. A frozen working set of 2,368 canonical transcriptional and epigenetic regulatory proteins is used as the anchor set.

Stage 0 parses the scope of every UniProt `Nucleus` annotation instead of propagating an alternative-isoform annotation to its canonical parent. A canonical/displayed sequence enters A when it has an unqualified Nucleus annotation, or when an isoform-specific Nucleus annotation applies to the displayed isoform. All other canonical proteins enter B.

| Set | Before anchor exclusion | After accession/gene exclusion |
|---|---:|---:|
| A: canonical/displayed Nucleus annotation | 5,669 | 3,418 |
| B: no strict canonical/displayed Nucleus annotation | 14,747 | 14,625 |

The A/B sets are disjoint and reconstruct all 20,416 canonical records. Anchor removal uses either an exact canonical accession match or an atomic-gene match.

## 2. Stage 1 candidate generation

### Route A: annotation-first discovery

The 3,418 non-anchor A proteins are classified by curated regulatory evidence, direct functional evidence, complex association and research/function profile. Three nested candidate sets are constructed:

```text
A1_WIDE    = all 3,418 minus 541 KNOWN_CURATED                    = 2,877
A2_NOVELTY = A1 minus 236 KNOWN_FUNCTIONAL_DIRECT                = 2,641
A3_FOCUSED = 440 complex-associated + 924 under-characterized
             + 1,093 mixed/moonlighting                          = 2,457
```

The public delivery uses A3_FOCUSED and the stricter Stage 2 RF2-PPI final90 background.

### Route B1: direct nuclear-access consensus

B1 applies two complementary sequence models to all 14,625 non-anchor B proteins:

```text
DeepLoc 2.0 Nucleus probability > 0.5
AND
corrected exact-length NLSExplorer max-window probability > 0.5
```

NLSExplorer uses 1,022-aa windows with 200-aa overlap and max-window protein aggregation. Windows reaching the classifier head have equal true length, preventing padding from changing the head output. The Stage 1 rule yields 1,350 B1 proteins.

### Route B2: signal-to-nucleus network

B2 is localization-agnostic. It removes the 1,350 B1 positives from the 14,625-protein B universe, leaving 13,275 eligible proteins. Candidates and anchors are connected through GO, UniProt Keyword and Reactome annotation terms plus the directed OmniPath network.

Four eligible-set percentiles receive equal weight:

1. annotation-graph distance to an anchor seed;
2. shared-seed IDF evidence;
3. directed OmniPath distance from candidate to anchor;
4. reachable anchors within three hops.

Stable sorting and `ceil(2% × 13,275)` select 266 Stage 1 B2 proteins.

### Route B3: alternative isoforms

B3 keeps the alternative isoform accession and true sequence as the unit of analysis. Nuclear evidence combines isoform-specific UniProt annotation with DeepLoc/NLSExplorer consensus. Candidate–anchor sequence interactions are then screened by PLM-interact and PPLM-PPI; the delivery endpoint requires PLM-interact `>0.99` without length truncation and PPLM-PPI `>0.9`.

The current delivery contains 112 B3 genes represented by 152 alternative isoforms and 490 retained isoform–anchor pairs.

## 3. Stage 2 interaction evidence

For canonical routes, Stage 2 connects each Stage 1 candidate to the 2,368 anchors using:

- HI-union experimental binary interactions, resolved at gene level;
- RF2-PPI official Data S4 (`final80`) and Data S3 (`final90`) sets, resolved primarily by UniProt accession while retaining match provenance.

`final80` and `final90` name the published RF2-PPI precision sets; they are not thresholds applied to a local probability column. The public delivery uses final90 for A and final80 for B1/B2.

| Route | Stage 1 input | Delivery pass | Retained pairs |
|---|---:|---:|---:|
| A / A3_FOCUSED | 2,457 | 954 | 3,530 |
| B1 | 1,350 | 429 | 1,593 |
| B2 | 266 | 85 | 344 |
| B3 | alternative isoforms | 112 genes / 152 isoforms | 490 |

HI-union support is gene-level and therefore does not identify a specific isoform. RF2-PPI support is a computational interaction prediction. Neither evidence type by itself constitutes a new biochemical binding experiment.

## 4. Stage 2.5 extension

Stage 2.5 reviews primary Stage 2 failures against final80-positive seeds from the same canonical input universe. An edge must jointly satisfy:

- BLOSUM62 local-alignment identity ≥70%;
- seed and candidate coverage ≥80%;
- identical domain order and copy number, with normalized boundary minimum IoU ≥0.80;
- a strict direct-function rule based on filtered GOA experimental evidence and/or direct human Reactome roles.

The v3 run screened 4,493 canonical proteins (1,862 positive seeds and 2,631 failures) and retained nine extension candidates. The extension preserves each candidate's primary Stage 2 PPI failure label; it supplies orthogonal homology/domain/function support rather than a candidate–anchor PPI claim.

## 5. Reproducibility boundaries

The repository includes current analysis, inference and validation code plus compact result artifacts. Raw database snapshots, licensed or externally hosted model weights, embeddings, training checkpoints and multi-gigabyte intermediate tables are intentionally excluded. Reproduction therefore requires obtaining the corresponding upstream resources and passing their local paths to the stage scripts.
