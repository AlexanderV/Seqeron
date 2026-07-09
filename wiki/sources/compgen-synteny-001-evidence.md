---
type: source
title: "Evidence: COMPGEN-SYNTENY-001 (Synteny / collinearity block detection — MCScanX DP scoring model)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-SYNTENY-001-Evidence.md
source_commit: 326807a15f6c06f4b66ac9411bb2e965713b6f25
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-SYNTENY-001

The validation-evidence artifact for test unit **COMPGEN-SYNTENY-001** — synteny / collinearity
block detection between two whole genomes, following the **MCScanX collinearity model** (Wang et
al. 2012). This is a **Comparative-genomics** family Evidence file and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. It is the
comparative-genomics counterpart to the chromosome-scale [[chrom-synt-001-evidence]], and it
reuses the shared synteny anchor [[synteny-and-rearrangement-detection]] rather than re-deriving
the syntenic-block definition — this file adds the concrete **MCScanX dynamic-programming scoring
parameters** behind that anchor's `FindSyntenyBlocks`. Its sibling COMPGEN units are
[[average-nucleotide-identity]], [[conserved-gene-clusters-common-intervals]],
[[ortholog-detection-reciprocal-best-hits]], [[genome-rearrangement-breakpoint-distance]], and the
end-to-end pipeline [[genome-comparison-core-dispensable]] that consumes these blocks for its
`OverallSynteny` fraction. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources (all rank 1 unless noted):**
  - **MCScanX** (Wang et al. 2012, *Nucleic Acids Res.* 40(7):e49, PMC3326336) — the collinearity
    DP model, extracted verbatim:
    - **Scoring recurrence:** `Score(v) = max(MatchScore(v), max(Score(u) + MatchScore(v) +
      GapPenalty × NumberofGaps(u,v)))` — a chain DP that rewards adjacent collinear anchor pairs
      and penalizes the distance between them.
    - **MatchScore = 50** per anchored gene pair; **GapPenalty = −1**.
    - **NumberofGaps(u,v)** = max intervening genes between anchors `u` and `v`; must be **< 25**
      (default `MAX_GAPS = 25`).
    - **Reporting threshold:** non-overlapping chains scoring **over 250** (i.e. **≥ 5 collinear
      gene pairs**, since `5 × 50 = 250`); default minimum block = **5** anchor pairs.
    - **Match collapsing:** consecutive BLASTP matches sharing a gene whose partners are < 5 genes
      apart collapse to the representative pair with the smallest E-value; **anchor E-value cutoff
      10⁻⁵**.
    - **Directionality:** matches are sorted in **both transcriptional directions**, so both
      **forward** and **inverted (reverse)** collinear blocks are detected.
  - **MCScanX** (Oxford Academic HTML) — the synteny-vs-collinearity distinction: *"Collinearity, a
    more specific form of synteny, requires conserved gene order."* Anchors are homologous
    (ortholog/paralog) gene pairs.
  - **Wikipedia "Synteny"** (rank 4, definitions + to locate primaries) — conserved/shared synteny =
    preserved co-localization of homologous genes across species; syntenic blocks detected by an
    MCScan-family DP over shared homologous genes, accounting for gene loss/gain.
- **Algorithm behaviour:** chain DP over an input **ortholog/anchor map** (anchor *generation* is
  out of scope — delegated to COMPGEN-ORTHO-001); reports non-overlapping forward and inverted
  blocks of ≥ 5 collinear anchors with gaps < 25, each tagged `IsInverted`.
- **Datasets (documented oracles):**
  - *Forward chain* — 5 adjacent anchors, identical order both genomes (0 gaps) → score `5×50 −
    1×0 = 250` → one forward block, `IsInverted = false`, GeneCount 5.
  - *Reverse chain* — 5 anchors with genome2 order reversed (4,3,2,1,0) → one block,
    `IsInverted = true`.
  - *Sub-threshold* — 4 adjacent anchors → score `4×50 = 200 ≤ 250` → **no block**.
  - *Gap over cutoff* — two anchor runs separated by ≥ 25 intervening genes → gap breaks the chain,
    neither sub-run reaches 5 pairs → no block.
  - *Empty genome / empty ortholog map* → no blocks.

## Deviations and assumptions

Two **ASSUMPTIONs**, both source-backed, neither a correctness gap:

1. **Threshold boundary at exactly 250.** The paper says "scores over 250" yet also "at least 5
   collinear gene pairs", and a 5-pair zero-gap chain scores exactly 250. The unit adopts the
   paper's explicit **≥ 5-pair minimum** as the operative rule: report iff score ≥ MinChainScore
   (250) **and** anchors ≥ MinAnchors (5). Source-backed resolution of the wording tension.
2. **Anchors supplied as an ortholog map.** MCScanX derives anchors from BLASTP (E ≤ 1e-5) with
   near-duplicate collapsing; this repo delegates anchor/ortholog identification to a separate unit
   (COMPGEN-ORTHO-001) and accepts an `orthologMap` input. The chaining algorithm under test is
   unchanged; anchor generation is out of scope.

No contradictions among sources — the two MCScanX renderings and Wikipedia give the same DP /
collinearity model. A `VisualizeSynteny` smoke test and an O(n²) per-block invariant (GeneCount ≥ 5,
coordinates within parent gene bounds) round out the recommended coverage.
</content>
