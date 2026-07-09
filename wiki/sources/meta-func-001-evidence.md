---
type: source
title: "Evidence: META-FUNC-001 (functional prediction — BLAST bit-score/E-value homology transfer + hypergeometric ORA)"
tags: [validation, metagenomics]
doc_path: docs/Evidence/META-FUNC-001-Evidence.md
sources:
  - docs/Evidence/META-FUNC-001-Evidence.md
source_commit: ab6548a3b59cfe5105eb5aa7cb252e63ffed5de7
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: META-FUNC-001

The validation-evidence artifact for test unit **META-FUNC-001** — **functional prediction**: assign
biological function to a metagenome's predicted proteins by **homology-based annotation transfer**
(`MetagenomicsAnalyzer.PredictFunctions`) and score which pathways are over-represented among them by the
**hypergeometric ORA** (`FindPathwayEnrichment`). This is the **fifth ingested unit of the Metagenomics
family** and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern. The method is synthesized in its own concept,
[[functional-prediction]]; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (all mutually consistent, no contradictions):**
  - **Altschul et al. — NCBI BLAST tutorial "The Statistics of Sequence Similarity Scores"**
    (accessed 2026-06-13, authority rank 2) — the Karlin-Altschul E-value `E = K·m·n·e^(−λS)`, the
    normalized **bit score** `S' = (λS − ln K)/ln 2`, and the bit-score→E-value identity
    `E = m·n·2^(−S')`.
  - **NCBI `blast_stat.c` — `blosum62_values`** (rank 3, reference implementation) — the **ungapped**
    BLOSUM62 Karlin-Altschul parameters `λ = 0.3176`, `K = 0.134`, `H = 0.4012` (the gap-open =
    gap-extend = INT2_MAX row).
  - **NCBI BLOSUM62 matrix file** (rank 2) — the diagonal (self-match) scores used to score an exact
    signature self-alignment: A4 R5 N6 D6 C9 Q5 E5 G6 H8 I4 L4 K5 M5 F6 P7 S4 T5 **W11** Y7 V4.
  - **PNNL Proteomics Data Analysis tutorial §8.2 (ORA)** (rank 4; the underlying hypergeometric is
    standard) — the right-tail over-representation p-value
    `P(X ≥ x) = 1 − Σ_{i=0}^{x−1} C(M,i)·C(N−M,n−i)/C(N,n)` with N/M/n/x definitions and a worked example.

- **Extracted formulas & identities:** bit score linear in S with slope `λ/ln2 > 0` (strictly
  increasing); E-value `∝ e^(−λS)` (strictly decreasing in S); the algebraic identity
  `K·m·n·e^(−λS) = m·n·2^(−S')` verified to machine precision; hypergeometric right-tail p ∈ [0,1].

- **Documented corner / failure cases:**
  - Karlin-Altschul EVD requires a positive expected maximal score; an exact self-match of a non-empty
    segment over BLOSUM62 has `S > 0` (all diagonals ≥ 4).
  - ORA: `x = 0` (no overlap) ⇒ empty sum ⇒ `P = 1`; degenerate margins `M = 0` or `n = 0` ⇒ `P = 1`;
    **right-tail only** (under-representation is a separate lower-tail test).

- **Datasets (documented oracles):**
  - **Homology-transfer self-consistency:** signature `"WWW"` (3 × W, diagonal 11) → raw `S = 33` →
    bit `S' = 18.0202932787533…`; with `m = n = 3` both E-value forms give
    `3.3852730346546 × 10⁻⁵` (agree to machine precision).
  - **ORA worked example (PNNL §8.2):** `N = 8000, M = 400, n = 100, x = 20` → `P(X ≥ 20) = 7.88 × 10⁻⁸`
    (R: `phyper(19, 400, 7600, 100, lower.tail = FALSE)`).

## Recommended test coverage (from the Evidence file)

MUST: bit score for the `WWW` self-match (18.0202932787533); the two-form E-value equality
(3.3852730346546 × 10⁻⁵); the ORA p-value 7.88 × 10⁻⁸; best-hit selection = lowest E-value / highest bit
score. SHOULD: ORA `p = 1` when `x = 0` or `M = 0` or `n = 0`; null/empty inputs → no annotations / no
enrichment without throwing. COULD: monotonicity — larger raw score ⇒ larger bit score ⇒ smaller E-value.

## Deviations and assumptions

- **ASSUMPTION (ASM-01) — ungapped exact-match scoring model.** `PredictFunctions` transfers function by
  detecting an **exact occurrence** of a database signature in the query protein and scoring that ungapped
  self-match with BLOSUM62, so the source-backed ungapped `λ, K` apply directly and no gap model is
  needed. This affects **which hits are found** (divergent homologs a full Smith-Waterman/BLAST search
  would find are missed), **not** the bit-score/E-value formulas, which are used exactly as cited. The
  unit's canonical complexity `O(n × g)` (per-gene × per-database-entry scan) is consistent with signature
  matching rather than full alignment.

No source contradictions — the Altschul E-value/bit-score forms, the `blast_stat.c` ungapped parameters,
the BLOSUM62 diagonals, and the PNNL hypergeometric formula are mutually consistent.

## Scope note — the ORA half has its own dedicated unit

The hypergeometric ORA (`FindPathwayEnrichment` / `HypergeometricUpperTail`) appears here as **component B**
of functional prediction, but it also has its **own** test unit **META-PATHWAY-001**
(`docs/algorithms/Metagenomics/Pathway_Enrichment_ORA.md`, evidence `META-PATHWAY-001-Evidence.md`, **not
yet ingested**). META-FUNC-001 validates **only** `docs/algorithms/Metagenomics/Functional_Prediction.md`;
the shared ORA statistics are captured on [[functional-prediction]] and will be revisited when
META-PATHWAY-001 is ingested.
</content>
