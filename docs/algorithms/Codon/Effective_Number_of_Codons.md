# Effective Number of Codons (ENC / Nc)

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Usage Analysis |
| Test Unit ID | CODON-ENC-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

The effective number of codons (Nc, also ENC) measures synonymous codon-usage bias in a single coding sequence. It answers: "how many codons are effectively in use in this gene?" Nc ranges from 20 (extreme bias — exactly one codon used per amino acid) to 61 (no bias — every synonymous codon used equally) [1][2]. It is a deterministic, count-based statistic computed from one gene, requiring no reference set (unlike CAI). The implementation follows Wright's original 1990 estimator exactly as reproduced in Fuglsang (2004) [2].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Most amino acids are encoded by more than one synonymous codon. Organisms and genes differ in how evenly they use these synonyms. Nc quantifies this evenness using codon "homozygosity" — the probability that two randomly chosen codons for the same amino acid are identical — aggregated across the degeneracy classes of the standard genetic code [1][2].

### 2.2 Core Model

For an amino acid with `k` synonymous codons, total count `n`, and codon frequencies `p_i = n_i / n`, the codon homozygosity is (Wright 1990, Eq. 1 in [2]):

```
F̂ = ( n·Σ_{i=1..k} p_i² − 1 ) / ( n − 1 )
```

The effective number of codons for that amino acid is `N̂c(aa) = 1 / F̂` (Eq. 2 [2]). The gene-level value aggregates class averages `F̂_2, F̂_3, F̂_4, F̂_6` over the standard-code degeneracy classes (Eq. 3 [2]):

```
N̂c = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆
```

The constant `2` is the contribution of the two single-codon amino acids Met (ATG) and Trp (TGG); `9, 1, 5, 3` are the numbers of two-, three-, four- and six-fold degenerate amino acids in the standard (NCBI table 1) genetic code [3]. Stop codons are excluded.

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Codons within an amino acid follow the standard genetic code degeneracy classes (9 two-fold, 1 three-fold, 5 four-fold, 3 six-fold, 2 single) | A non-standard code changes the class numerators/constant and the result is no longer Wright's Nc [3] |
| ASM-02 | Each represented amino acid has n ≥ 2 codons so F̂ is defined | For n ≤ 1 the amino acid is skipped (denominator n−1) and its class falls back to the within-class average (Eq. 4) [2] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 20 ≤ Nc ≤ 61 | Extreme-bias / no-bias limits of Wright's estimator; upper value re-adjusted to 61 per Eq. 3 [2] |
| INV-02 | One codon per amino acid (each used ≥2×) ⇒ Nc = 20 | F̂ = 1 for every class ⇒ each N̂c(aa) = 1, sum = 9+1+5+3+2 = 20 [1][2] |
| INV-03 | Near-uniform usage ⇒ Nc re-adjusted to exactly 61 | Wright's overshoot rule caps Nc at 61 [2] |
| INV-04 | Deterministic | Pure function of codon counts |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | Nc (ENC) | CAI |
|--------|----------|-----|
| Reference set required | No (single gene) | Yes (highly expressed genes) |
| Range | 20–61 | 0–1 |
| Measures | Evenness of synonymous usage | Adaptation to a reference |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence | `string` or `DnaSequence` | required | Coding DNA sequence | Read in frame as consecutive non-overlapping triplets from index 0; non-ACGT codons skipped; case-insensitive |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | Effective number of codons, in [20, 61]; 0 for null/empty string |

### 3.3 Preconditions and Validation

`CalculateEnc(DnaSequence)` throws `ArgumentNullException` for null. `CalculateEnc(string)` returns 0 for null/empty. Input is upper-cased; codons are read 0-based in non-overlapping triplets; a trailing partial codon (length < 3) is ignored; codons containing any non-ACGT character are skipped (consistent with `CountCodons`). Amino acids with total count ≤ 1 are skipped (F̂ undefined).

## 4. Algorithm

### 4.1 High-Level Steps

1. Count valid ACGT codons in frame (reuse `CountCodons`).
2. For each amino acid with degeneracy > 1 and n ≥ 2, compute F̂ by Wright Eq. (1).
3. Average F̂ within each degeneracy class (2, 3, 4, 6) — Eq. (4).
4. If the 3-fold class (isoleucine) is unestimable, use `F̂₃ = (F̂₂ + F̂₄)/2` — Eq. (5a).
5. Aggregate via Eq. (3): `Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`; a class with no estimable F̂ contributes its full codon count.
6. Clamp to [20, 61] (upper re-adjustment per Wright; lower bound structural).

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

Degeneracy-class numerators (standard genetic code [3]): two-fold = 9, three-fold = 1 (Ile), four-fold = 5, six-fold = 3, single = 2 (Met, Trp). These are named constants in the implementation, each cited to Wright/Fuglsang.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateEnc | O(n) | O(1) | n = sequence length; codon table is fixed size (64) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonUsageAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs)

- `CodonUsageAnalyzer.CalculateEnc(string)`: canonical Wright 1990 computation on a raw sequence.
- `CodonUsageAnalyzer.CalculateEnc(DnaSequence)`: delegates to the string overload via `.Sequence`.

### 5.2 Current Behavior

F̂ is computed from frequencies `p_i = n_i/n` (Eq. 1), not from raw counts. Class averages substitute for absent amino acids (Eq. 4); the isoleucine 3-fold fallback (Eq. 5a) applies when no isoleucine is present but the 2- and 4-fold classes are estimable. A degeneracy class with no estimable F̂ at all contributes its full codon count. The result is re-adjusted to ≤ 61 (Wright) and floored at 20. This unit does not perform substring search, so the repository suffix tree is **not** applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Eq. (1) homozygosity `F̂ = (n·Σ p_i² − 1)/(n − 1)` with `p_i = n_i/n` [2].
- Eq. (3) aggregation `Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆` with standard-code class counts [2][3].
- Eq. (4) within-class averaging for absent amino acids [2].
- Eq. (5a) isoleucine fallback `F̂₃ = (F̂₂ + F̂₄)/2` [2].
- Upper re-adjustment of Nc to 61 [2].

**Intentionally simplified:**

- Lower clamp at 20: Wright/Fuglsang only prescribe the upper re-adjustment; **consequence:** a structurally-floored value of 20 is returned for degenerate inputs, consistent with the published range but not an explicit Wright instruction.

**Not implemented:**

- The Fuglsang (2006) N̂c "sampling-without-replacement rounding" variant; **users should rely on:** the standard Wright Nc returned here (the most widely used estimator [2]).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null `DnaSequence` | `ArgumentNullException` | Contract |
| empty / null string | 0 | Degenerate input; no Wright rule for empty gene |
| amino acid with n ≤ 1 | skipped; class uses Eq. 4 average | F̂ undefined (n−1) [2] |
| isoleucine absent | F̂₃ = (F̂₂ + F̂₄)/2 | Eq. 5a [2] |
| near-uniform short gene | re-adjusted to 61 | Eq. 3 overshoot rule [2] |
| non-ACGT codon | skipped | consistent with `CountCodons` |

### 6.2 Limitations

Standard genetic code only (ASM-01). Nc overestimates for very short genes [2]. The result is undefined-but-clamped when no degeneracy class is estimable (returns 20–61 by the full-count fallback).

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
double nc = CodonUsageAnalyzer.CalculateEnc("ATGGCTGCAGCTGCA"); // in [20, 61]
```

**Numerical / biological walk-through:**

Gene with only Phe (TTT×3, TTC×1): n=4, p=(0.75,0.25), Σp²=0.625, F̂=(4·0.625−1)/3=0.5, N̂c(Phe)=2. No other classes estimable ⇒ they contribute their full counts: Nc = 2 + 9/0.5 + 1 + 5 + 3 = 29.0.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonUsageAnalyzer_CalculateEnc_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonUsageAnalyzer_CalculateEnc_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [CODON-ENC-001-Evidence.md](../../../docs/Evidence/CODON-ENC-001-Evidence.md)

## 8. References

1. Wright, F. 1990. The 'effective number of codons' used in a gene. *Gene* 87(1):23–29. https://doi.org/10.1016/0378-1119(90)90491-9
2. Fuglsang, A. 2004. The 'effective number of codons' revisited. *Biochemical and Biophysical Research Communications* 317(3):957–964. https://doi.org/10.1016/j.bbrc.2004.03.138
3. Fuglsang, A. 2006. Estimating the 'effective number of codons': the Wright way of determining codon homozygosity leads to superior estimates. *Genetics* 172(2):1301–1307. https://academic.oup.com/genetics/article/172/2/1301/5923091
