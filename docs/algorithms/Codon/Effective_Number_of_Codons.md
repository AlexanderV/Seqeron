# Effective Number of Codons (ENC / Nc)

| Field | Value |
|-------|-------|
| Algorithm Group | Codon Usage Analysis |
| Test Unit ID | CODON-ENC-001 |
| Related Projects | Seqeron.Genomics.MolTools |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

The effective number of codons (Nc, also ENC) measures synonymous codon-usage bias in a single coding sequence. It answers: "how many codons are effectively in use in this gene?" Nc ranges from 20 (extreme bias — exactly one codon used per amino acid) to 61 (no bias — every synonymous codon used equally) [1][2]. It is a deterministic, count-based statistic computed from one gene, requiring no reference set (unlike CAI). The implementation follows Wright's original 1990 estimator as reproduced in Fuglsang (2004) [2] and as implemented by the de-facto reference tool CodonW 1.4.4 (`enc_out`) [4], including its "Nc not calculated" rule, and is genetic-code aware.

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

The constant `2` is the contribution of the two single-codon amino acids Met (ATG) and Trp (TGG); `9, 1, 5, 3` are the numbers of two-, three-, four- and six-fold degenerate amino acids in the standard (NCBI table 1) genetic code [3]. Stop codons are excluded. For another genetic code the same formula is applied to that code's classes, `Nc = K₁ + Σ_z K_z / F̂_z` (K_z = number of amino acids with z codons; CodonW derives them from the selected code, `-enc -code`) [4]; e.g. table 2 has 12 two-fold, 6 four-fold and 2 six-fold amino acids, table 3 an 8-fold Thr.

### 2.3 Modeling Assumptions (Optional)

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | Synonymous classes are those of the selected genetic code (default NCBI table 1: 9 two-fold, 1 three-fold, 5 four-fold, 3 six-fold, 2 single) | Using the wrong code mis-assigns codons to families [4] |
| ASM-02 | An amino acid is estimable when n ≥ 2 and F̂ > 0 | For n ≤ 1 (denominator n−1) or F̂ = 0 (every observed codon used once; CodonW `bb > 0.0000001`) the amino acid is left out and its class uses the within-class average (Eq. 4) [2][4] |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | 20 ≤ Nc ≤ 61 (table 1) whenever Nc is calculable; otherwise 0 | Extreme-bias / no-bias limits; upper value re-adjusted to 61 [2][4]. For other codes the upper limit is that code's sense-codon count |
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
| sequence | `string` or `DnaSequence` | required | Coding DNA or RNA sequence | Read in frame as consecutive non-overlapping triplets from index 0; U read as T; triplets with other symbols skipped without shifting the frame; case-insensitive |
| code | `GeneticCode` | `GeneticCode.Standard` | Genetic code defining the synonymous classes | non-null; any of the 27 NCBI tables |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | Effective number of codons, in [20, 61] (table 1); **0 when Nc cannot be calculated** (null/empty input, or a synonymous class with no estimable amino acid — CodonW prints `*****`) |

### 3.3 Preconditions and Validation

`CalculateEnc(DnaSequence[, GeneticCode])` throws `ArgumentNullException` for a null sequence; every overload throws it for a null code. `CalculateEnc(string[, GeneticCode])` returns 0 for null/empty. Input is upper-cased and U is read as T; codons are read 0-based in non-overlapping triplets; a trailing partial codon (length < 3) is ignored; codons containing any other character are skipped (consistent with `CountCodons`). Amino acids with total count ≤ 1 or F̂ = 0 are not estimable.

## 4. Algorithm

### 4.1 High-Level Steps

1. Count valid codons in frame (the canonical `CountCodons` core).
2. Group the sense codons of the genetic code into synonymous families; K_z = number of amino acids with z codons.
3. For each amino acid with z > 1 and n ≥ 2, compute F̂ by Wright Eq. (1); keep it only if F̂ > 0.0000001 (CodonW).
4. Average F̂ within each degeneracy class z — Eq. (4).
5. If the class z = 3 has a single amino acid (Ile) and it is unestimable, use `F̂₃ = (F̂₂ + F̂₄)/2` — Eq. (5a).
6. If any other class z > 1 has no estimable amino acid, return 0 (Nc not calculated; CodonW `*****`).
7. Aggregate `Nc = K₁ + Σ_z K_z / F̂_z` (table 1: `2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆`).
8. Re-adjust values above the number of sense codons of the code (61 for table 1) down to it; the lower bound 20 (= number of amino acids) only guards floating-point rounding.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

Degeneracy-class numerators are derived from `GeneticCode.CodonTable` (standard code [3]: two-fold = 9, three-fold = 1 (Ile), four-fold = 5, six-fold = 3, single = 2 (Met, Trp)). Codons that are context-dependent stops in tables 27/28/31 belong to the amino acid they encode (as for RSCU/CAI).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateEnc | O(n) | O(1) | n = sequence length; codon table is fixed size (64) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [CodonUsageAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs)

- `CodonUsageAnalyzer.CalculateEnc(string, GeneticCode)`: canonical Wright 1990 / CodonW computation on a raw sequence.
- `CodonUsageAnalyzer.CalculateEnc(string)`, `CalculateEnc(DnaSequence)`, `CalculateEnc(DnaSequence, GeneticCode)`: same core (standard code by default).
- `CodonUsageAnalyzer.GetStatistics(...)`.`Enc` and MCP `effective_number_of_codons` / `codon_usage_statistics` delegate to this core (standard code).

### 5.2 Current Behavior

F̂ is computed from frequencies `p_i = n_i/n` (Eq. 1), not from raw counts. Class averages substitute for absent amino acids (Eq. 4); F̂ = 0 is not an estimate (CodonW). The isoleucine 3-fold fallback (Eq. 5a) applies when the single 3-fold amino acid is unestimable but the 2- and 4-fold classes are estimable. Any other empty class ⇒ 0 (not calculated). The result is re-adjusted to ≤ the code's sense-codon count (61 for table 1). Cross-checked against the compiled CodonW 1.4.4 binary on 6456 gene × code cases (codes 1,2,3,4,5,6,9,10; 3136 "not calculated" cases all 0; all others equal to 2 dp) and against a Python port of `enc_out` on all 27 tables (21789 cases, exact). This unit does not perform substring search, so the repository suffix tree is **not** applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Eq. (1) homozygosity `F̂ = (n·Σ p_i² − 1)/(n − 1)` with `p_i = n_i/n` [2].
- Eq. (3) aggregation `Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆` with standard-code class counts [2][3].
- Eq. (4) within-class averaging for absent amino acids [2].
- Eq. (5a) isoleucine fallback `F̂₃ = (F̂₂ + F̂₄)/2` [2].
- Upper re-adjustment of Nc to 61 [2][4]; "Nc not calculated" when a synonymous class is empty (except Ile) [4].

**Documented divergence:**

- For genetic codes whose sense-codon count is not 61, the upper re-adjustment uses that count (uniform usage gives Nc = Σ K_z·z), as codonbias 0.5.0 does (`min(len(P), ENC)`); CodonW hard-codes 61 for every code. Identical for tables 1 and 11.

**Not implemented:**

- The Fuglsang (2006) N̂c "sampling-without-replacement rounding" variant; **users should rely on:** the standard Wright Nc returned here (the most widely used estimator [2]).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| null `DnaSequence` | `ArgumentNullException` | Contract |
| empty / null string | 0 | Nc not calculated |
| amino acid with n ≤ 1, or F̂ = 0 | left out; class uses Eq. 4 average | F̂ undefined (n−1) [2]; CodonW `bb > 0.0000001` [4] |
| a synonymous class (other than a lone Ile) with no estimable amino acid | 0 | Nc not calculated (CodonW `*****`, Wright 1990) [4] |
| isoleucine absent | F̂₃ = (F̂₂ + F̂₄)/2 | Eq. 5a [2] |
| near-uniform short gene | re-adjusted to 61 | Eq. 3 overshoot rule [2] |
| non-ACGT(U) codon | skipped, frame kept | consistent with `CountCodons` |
| RNA / lower case | same as DNA upper case | CodonW reads U as T |

### 6.2 Limitations

Nc overestimates for very short genes [2]; the Novembre (2002) background-corrected Nc′ and the Sun, Yang & Xia (2013) variant are not implemented (not claimed). A return value of 0 means "not calculated" and must not be averaged with real Nc values.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
double nc = CodonUsageAnalyzer.CalculateEnc(gene);                               // 20..61, or 0 if not calculable
double mt = CodonUsageAnalyzer.CalculateEnc(gene, GeneticCode.GetByTableNumber(2)); // vertebrate mitochondrial classes
```

**Numerical / biological walk-through:**

Gene M3 (TTT×4 TTC; CTG×3 CTC×2 TTA; ATT×3 ATC×2 ATA; GTG×4 GTC; AGC×3 TCT×2 TCA; CGC×4 CGT×2; GGC×3 GGT×2 GGA): F̂₂ = 0.6, F̂₃ = 0.2667, F̂₄ = 0.4333, F̂₆ = 0.3333 ⇒ Nc = 2 + 15 + 3.75 + 11.538 + 9 = 41.288461538461526 (CodonW 1.4.4: 41.29). Adding His as CAT+CAC (F̂ = 0) leaves Nc unchanged. A gene with only Phe (TTT×3, TTC×1) has empty 3/4/6-fold classes ⇒ not calculated (0; CodonW `*****`).

### 7.3 Related Tests, Evidence, or Documents

- Tests: [CodonUsageAnalyzer_CalculateEnc_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/MolTools/CodonUsageAnalyzer_CalculateEnc_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [CODON-ENC-001-Evidence.md](../../../docs/Evidence/CODON-ENC-001-Evidence.md)

## 8. References

1. Wright, F. 1990. The 'effective number of codons' used in a gene. *Gene* 87(1):23–29. https://doi.org/10.1016/0378-1119(90)90491-9
2. Fuglsang, A. 2004. The 'effective number of codons' revisited. *Biochemical and Biophysical Research Communications* 317(3):957–964. https://doi.org/10.1016/j.bbrc.2004.03.138
3. Fuglsang, A. 2006. Estimating the 'effective number of codons': the Wright way of determining codon homozygosity leads to superior estimates. *Genetics* 172(2):1301–1307. https://academic.oup.com/genetics/article/172/2/1301/5923091
4. Peden, J.F. 1999. *Analysis of codon usage* (PhD thesis, Univ. Nottingham) and CodonW 1.4.4 source, `codon_us.c` `enc_out`; `README_indices.txt` ("When there are no amino acids in a synonymous family, Nc is not calculated …"). https://codonw.sourceforge.net/
