---
type: concept
title: "Expression quantification (TPM / FPKM / RPKM + quantile normalization)"
tags: [transcriptome, algorithm]
mcp_tools:
  - calculate_tpm
sources:
  - docs/algorithms/Transcriptome/Expression_Quantification.md
  - docs/Evidence/TRANS-EXPR-001-Evidence.md
source_commit: a6c1fb6eec93a97474ade58430d30f5b57a63223
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-expr-001-evidence
      evidence: "Test Unit ID: TRANS-EXPR-001, Algorithm: Expression Quantification (TPM, FPKM/RPKM, Quantile Normalization); Methods TranscriptomeAnalyzer.CalculateTPM / CalculateFPKM / QuantileNormalize"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:differential-expression
      source: trans-expr-001-evidence
      evidence: "Both TRANS-EXPR-001 and TRANS-DIFF-001 are TranscriptomeAnalyzer units of the Transcriptome/RNA-seq family; per-gene normalized expression (TPM/FPKM) or a quantile-normalized matrix is the upstream input a two-group DE test consumes."
      confidence: high
      status: current
---

# Expression quantification (TPM / FPKM / RPKM + quantile normalization)

**Convert raw RNA-seq read counts into comparable expression values.** Raw counts are not
comparable across genes (a longer gene collects more reads at the same expression level) or across
samples (different library sizes / sequencing depth). This unit packages the three standard
corrections: **TPM** and **FPKM/RPKM** normalize a single sample's counts by gene length and library
size, and **quantile normalization** forces a multi-sample matrix onto one common distribution. A
**Transcriptome / RNA-seq family** sibling of the family anchor [[differential-expression]] — a
distinct *quantification / normalization* method, upstream of the two-group DE test. Validated under
test unit **TRANS-EXPR-001**; the record is [[trans-expr-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern.

Impl `TranscriptomeAnalyzer` (`Seqeron.Genomics.Annotation`, the Annotation server):
`CalculateTPM(...)` → per-gene `GeneExpression`, `CalculateFPKM(rawCount, length, totalReads)` →
`double`, and `QuantileNormalize(...)` → the normalized matrix. MCP tools `CalculateTpm` /
`QuantileNormalize`.

## Method contract (algorithm spec)

Genuinely-distinct implementation detail from the canonical PRIMARY spec
(`docs/algorithms/Transcriptome/Expression_Quantification.md`, TRANS-EXPR-001, *Production*), in
`Seqeron.Genomics.Annotation/TranscriptomeAnalyzer.cs`:

| Method | Signature | Returns |
| --- | --- | --- |
| `CalculateTPM` | `(IEnumerable<(string GeneId, double RawCount, int Length)> geneCounts)` | `IEnumerable<GeneExpression>` — per-gene record; TPM sums to 10⁶ |
| `CalculateFPKM` | `(double rawCount, int length, double totalReads)` | `double` — single-gene FPKM |
| `QuantileNormalize` | `(IEnumerable<IEnumerable<double>> samples)` | `IEnumerable<IReadOnlyList<double>>` — one normalized vector per sample, original positions preserved |

- **`CalculateTPM` also fills the `FPKM` field** of each `GeneExpression`, using the **sample's total
  raw count as N** (the FPKM depth normalizer) — so a single `CalculateTPM` pass returns both the TPM
  and its matching FPKM per gene.
- **`QuantileNormalize` uses the first sample's length as the gene count** (all samples must be equal
  length); inputs are materialized to lists and each sample is re-sorted twice (once to build the rank
  means, once to place them back) for clarity rather than caching.
- **Input constraints:** `RawCount ≥ 0`, `Length > 0` for a meaningful rate; FPKM's `length > 0` and
  `totalReads > 0` (else returns 0). No exceptions on degenerate input — empty `geneCounts` / empty
  `samples` / zero genes yield an empty sequence.
- **INV-04 / INV-05 (QN, explicit):** quantile normalization **preserves within-column rank order**
  (rank → rank-mean is non-decreasing), and an untied column maps to a **permutation of the same
  rank-mean multiset** (each rank is assigned exactly one rank mean).
- **Deviation (fixed):** a pre-fix QN implementation ignored tied ranks (tied values got the wrong
  rank means); corrected to the Bolstad tie-average rule.
- **Complexity:** `CalculateTPM` **O(n)** time / **O(n)** space (n genes, one pass for rates, one for
  output); `CalculateFPKM` **O(1)** / **O(1)**; `QuantileNormalize` **O(s·n log n)** / **O(s·n)** for
  s samples of n values (the per-sample sort dominates). No substring/pattern search — the repository
  suffix tree is N/A.
- **Not implemented (defer to DE tooling):** library-size scaling factors (edgeR TMM, DESeq2
  median-of-ratios) — out of scope for this within-sample/distribution-matching unit.

## 1. TPM — transcripts per million

Divide count by length first (**reads per kilobase, RPK**), then rescale so each sample **sums to
exactly 10⁶**:

```
RPK_i = X_i / l_i
TPM_i = (X_i / l_i) / Σ_j (X_j / l_j) · 10⁶
```

`X_i` = read/fragment count for gene i, `l_i` = its (effective) length. Because the length division
happens **before** the per-sample rescaling, the total is a fixed 10⁶ in every sample — the
**sum-to-a-million invariant** (Wagner, Kin & Lynch 2012; Zhao, Ye & Stanton 2020). Equivalently the
**average** TPM = 10⁶ / (number of annotated transcripts), a constant. This is exactly the property
RPKM lacks, which is why TPM is preferred for **within-sample** proportion comparisons.

## 2. FPKM / RPKM — length- and depth-normalized

RPKM (single-end reads) and FPKM (paired-end fragments) are the **same formula** — normalize by both
gene length (per kilobase) and library size (per million mapped reads):

```
FPKM_i = X_i / ( (l_i/10³) · (N/10⁶) ) = X_i · 10⁹ / (l_i · N)
RPKM   = 10⁹ · (reads mapped to transcript) / (total reads · transcript length)   (identical)
```

`N` = total sequenced fragments/reads in the sample. **TPM is FPKM rescaled to a million:**
`TPM_i = FPKM_i / Σ_j FPKM_j · 10⁶` — so unlike TPM, an FPKM sample does **not** sum to a fixed
constant (the flaw Zhao 2020 warns against when comparing FPKM across samples).

## 3. Quantile normalization — cross-sample distribution matching

Force every sample (column) onto one **shared distribution** so their quantiles match, with **no
reference sample** needed (Bolstad et al. 2003):

1. **Sort** each column ascending, remembering original positions.
2. **Rank mean** = arithmetic mean across columns of the values at each rank (highest → mean of the
   column-maxima, etc.).
3. **Re-place** each rank mean back at its original position.

**Tie rule:** values tied within a column receive the **average of the rank means they would
otherwise span** — not each rank separately.

## Invariants and edge cases

- **INV (TPM):** for any non-degenerate input, `Σ_i TPM_i = 10⁶` and every `TPM_i ≥ 0`.
- **TPM = FPKM rescaled:** `TPM_i = FPKM_i / Σ FPKM_j · 10⁶`.
- **All-zero counts** → `Σ(X/l) = 0` → TPM is 0/0 (undefined); the impl emits **TPM = 0 for every
  gene** (ASSUMPTION — a degenerate-input convention, no source specifies it; never affects a
  non-zero-denominator result).
- **Non-positive length or total reads** → `X/l` undefined → **FPKM = 0**, excluded from the RPK sum.
- **Empty input → empty output** for all three methods.
- **Effective length = annotated length** (ASSUMPTION — the cited formulas use effective length
  `l̃_i`; the public API takes the annotated length directly, the standard `l̃_i = l_i` substitution;
  omits only the optional fragment-length correction, not the formula structure).

Worked oracles (from [[trans-expr-001-evidence]]):
- **TPM** three genes A(X=10,l=2000), B(X=20,l=4000), C(X=30,l=1000): RPK 0.005/0.005/0.030,
  ΣRPK=0.04 → TPM **(125000, 125000, 750000)**, sum **1000000** (invariant). Equal-length genes give
  TPM directly proportional to counts.
- **FPKM** X=1000, l=2000, N=10⁶ → `1000·10⁹/(2000·10⁶)` = **500**.
- **Quantile norm** columns C1=(5,2,3,4), C2=(4,1,4,2), C3=(3,4,6,8): rank means 2.0 / 3.0 / 4.666… /
  5.666…; C2's two tied `4`s (rows A, C) span the top two ranks → **both become (4.666…+5.666…)/2 =
  5.166…** (verbatim final matrix shows 5.17 for both). C1 out = (5.67, 2.00, 3.00, 4.67).

## Relationship to the rest of the RNA-seq family

- **[[differential-expression]]** (TRANS-DIFF-001) — the family anchor and the natural **downstream
  consumer**: a normalized expression matrix (or per-sample TPM/FPKM) is what the two-group log2FC +
  Welch-t + BH test operates on. Quantile normalization is one of the cross-sample corrections that
  makes such a comparison valid; TPM/FPKM are within-sample corrections.
- **Cross-sample caveat:** TPM and RPKM/FPKM are **within-sample relative** measures; Zhao, Ye &
  Stanton (2020) document that both are *misused* when compared directly across samples or protocols —
  the very failure mode quantile normalization (or a proper DE model) is meant to address.

## Scope and limitations

A [[research-grade-limitations|research-grade]] correctness reference for TPM, FPKM/RPKM, and
quantile normalization. It is the **classic count-normalization** layer — **not** an
effective-length / fragment-model estimator (kallisto/salmon) nor the TMM/median-of-ratios size
factors (edgeR/DESeq2). Two source-backed assumptions (all-zero TPM → 0; effective length = annotated
length). All formulas match their primary sources (Wagner 2012; Zhao/Ye/Stanton 2020; Pimentel 2014;
Bolstad 2003) with **no source contradictions**. **Not for clinical use.**
