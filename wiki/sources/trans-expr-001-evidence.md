---
type: source
title: "Evidence: TRANS-EXPR-001 (Expression quantification — TPM / FPKM / RPKM + quantile normalization)"
tags: [validation, transcriptome]
doc_path: docs/Evidence/TRANS-EXPR-001-Evidence.md
sources:
  - docs/Evidence/TRANS-EXPR-001-Evidence.md
source_commit: deb32560df90dfd97221f5218e71f6fd6cf3b2fd
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: TRANS-EXPR-001

The validation-evidence artifact for test unit **TRANS-EXPR-001** — **RNA-seq expression
quantification**: TPM (transcripts per million), FPKM/RPKM (fragments/reads per kilobase per million),
and quantile normalization. A **Transcriptome / RNA-seq family** unit and one instance of the
templated per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. The algorithm is
synthesized in its own concept, [[expression-quantification]] (a distinct quantification/normalization
method from the family anchor [[differential-expression]]); [[test-unit-registry]] tracks the unit.
Impl `TranscriptomeAnalyzer.CalculateTPM` / `CalculateFPKM` / `QuantileNormalize`
(`Seqeron.Genomics.Annotation`; MCP tools `CalculateTpm` / `QuantileNormalize`).

## What this file records

- **Online sources (mutually consistent, no contradictions):**
  - **Wagner, Kin & Lynch (2012)**, *Theory in Biosciences* (rank 1) — proposes **TPM** as a fix for
    RPKM's inconsistency across samples; RPKM measures relative molar concentration but does not
    respect its invariance property. The **average TPM = 10⁶ / (number of annotated transcripts)** is
    a constant → TPM values within a sample **sum to exactly 10⁶**.
  - **Zhao, Ye & Stanton (2020)**, *RNA* (PMC7373998, rank 1) — verbatim formulas:
    `RPKM = 10⁹ · reads / (total reads · length)`;
    `TPM = 10⁶ · (reads/length) / Σ(reads/length)`; `TPM = 10⁶ · RPKM / Σ(RPKM)`. Restates the
    sum-to-a-million invariant, and warns TPM/RPKM are **misused when compared across samples/protocols**.
  - **Pimentel (2014)** blog review (rank 3, corroborating) — `TPM_i = (X_i/l_i)/Σ_j(X_j/l_j)·10⁶` and
    `FPKM_i = X_i/((l_i/10³)(N/10⁶)) = X_i·10⁹/(l_i·N)`; `TPM_i = FPKM_i/Σ_j FPKM_j·10⁶`
    (X = counts, l = effective length, N = total fragments).
  - **Wikipedia "Quantile normalization"** citing **Bolstad et al. 2003** (rank 4) — sort each column,
    set each rank to the mean of the values at that rank across columns, re-place at original
    positions; **tie rule**: tied values get the mean of the rank means they would span.

- **Documented corner cases / failure modes:** **all-zero counts** → Σ(reads/length)=0 → TPM = 0/0
  undefined (impl emits TPM = 0 for all genes); **zero / non-positive gene length** → reads/length
  undefined → FPKM = 0, excluded from RPK; **cross-sample comparison caveat** (interpretation, not an
  output rule); quantile-norm **tied ranks** averaged; **empty matrix** → undefined (no rank means).

- **Datasets (documented oracles):**
  - *TPM three-gene* — A(X=10, l=2000), B(X=20, l=4000), C(X=30, l=1000): RPK 0.005/0.005/0.030,
    ΣRPK=0.04 → **(125000, 125000, 750000)**, sum **1000000** (invariant).
  - *FPKM single-gene* — X=1000, l=2000, N=10⁶ → `1000·10⁹/(2000·10⁶)` = **500**.
  - *Quantile normalization* — C1=(5,2,3,4), C2=(4,1,4,2), C3=(3,4,6,8); rank means 2.0 / 3.0 /
    14/3=4.666… / 17/3=5.666…; tie mean (r2,r3)=31/6=**5.166…** — both tied `4`s in C2 (rows A, C) →
    5.166… (final matrix 5.17); C1 out = (5.67, 2.00, 3.00, 4.67).

- **Test-coverage recommendations:** MUST — TPM three-gene = (125000,125000,750000) summing to 10⁶;
  TPM sums to 10⁶ for any non-degenerate input (invariant); FPKM single-gene = 500; quantile norm
  reproduces the Wikipedia worked-example matrix incl. C2 tie handling. SHOULD — empty → empty (all
  three); non-positive length / total reads → FPKM = 0. COULD — all-zero counts → all TPM = 0.

## Deviations and assumptions

- **ASSUMPTION (all-zero TPM → 0):** no source specifies the output when Σ(reads/length)=0; the
  mathematically defined quantity is 0/0. The impl emits TPM = 0 for all genes — a degenerate-input
  convention that never affects a non-zero-denominator result (every such case is fully defined by the
  cited formula).
- **ASSUMPTION (effective length = annotated length):** the cited formulas use effective length `l̃_i`
  (Pimentel 2014); the public API takes the annotated transcript length directly. `l̃_i = l_i` is the
  standard within-sample substitution — it omits only the optional fragment-length correction, not the
  formula structure.

No source contradictions — Wagner 2012, Zhao/Ye/Stanton 2020, Pimentel 2014, and Bolstad 2003 are
mutually consistent (Mortazavi 2008 cited for the original RPKM). This is the **classic
count-normalization** layer, not an effective-length fragment model (kallisto/salmon) or TMM/
median-of-ratios size factors (edgeR/DESeq2). Research-grade, not for clinical use.
