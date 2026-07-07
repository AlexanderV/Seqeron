# Coding Potential Calculation (CPAT Hexamer Usage-Bias Score)

| Field | Value |
|-------|-------|
| Algorithm Group | Extended Annotation |
| Test Unit ID | ANNOT-CODING-001 |
| Related Projects | Seqeron.Genomics.Annotation |
| Implementation Status | Framework |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Coding potential scoring distinguishes protein-coding sequence from non-coding sequence using intrinsic sequence statistics. This implementation realizes the **hexamer usage-bias** measure of CPAT [1]: the mean per-hexamer log-likelihood ratio of in-frame hexamer frequencies between a coding and a non-coding training set [2]. It is a probabilistic, alignment-free measure: positive scores indicate a coding sequence, negative scores a non-coding one [1]. Because it depends on user-supplied (organism-specific) frequency tables, it is a *framework* component — the caller provides the trained tables.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Coding regions exhibit non-random hexamer (6-mer) usage caused by codon-usage bias and dependence between adjacent amino acids [1]. Comparing the frequency of each in-frame hexamer in coding versus non-coding training data yields a discriminative log-likelihood signal; CPAT uses this hexamer score as one feature of its logistic-regression classifier [1]. An alternative composite measure is the Fickett TESTCODE statistic [4]; it is not implemented here (§5.3).

### 2.2 Core Model

For a sequence `S`, extract its in-frame hexamers `H1, H2, … , Hm` (window size 6, step 3, starting at offset 0) [2]. With coding frequency table `F` and non-coding table `F'`, the score is the mean per-hexamer contribution [2]:

- if `F(k) > 0` and `F'(k) > 0`: contribution = `ln( F(k) / F'(k) )` (natural log) [2];
- if `F(k) > 0` and `F'(k) = 0`: contribution = `+1` [2];
- if `F(k) = 0` and `F'(k) > 0`: contribution = `−1` [2];
- if `F(k) = 0` and `F'(k) = 0`: skipped (`continue`) and **not counted** [2];
- a hexamer absent from either table is skipped and not counted [2].

Score `= ( Σ contributions ) / (number of scored hexamers)` [2]. Positive ⇒ coding, negative ⇒ non-coding [1].

### 2.3 Modeling Assumptions

| ID | Assumption | Consequence if Violated |
|----|------------|--------------------------|
| ASM-01 | The supplied tables represent the coding/non-coding populations of the analysed organism [1] | Scores lose discriminative power; sign may mislead |
| ASM-02 | Both tables use the same units (raw counts or proportions) [2] | Score shifts by the additive constant `ln(ΣF / ΣF')` |
| ASM-03 | The reading frame of interest is frame 0 of the supplied sequence [2] | Out-of-frame hexamers are not scored; pass the ORF in frame |

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Score = (Σ per-hexamer contributions) / (count of scored hexamers) | Definition [2] |
| INV-02 | Per-hexamer contribution = `ln(F(k)/F'(k))` when both > 0 | `kmer_ratio` [2] |
| INV-03 | Coding-only hexamer ⇒ +1; non-coding-only ⇒ −1 | `kmer_ratio` branches [2] |
| INV-04 | Coding-biased tables ⇒ score > 0; non-coding-biased ⇒ score < 0 | Sign convention [1] |
| INV-05 | Only in-frame full-length hexamers (offset 0, step 3, length = wordSize) are scored | `word_generator` [2] |
| INV-06 | `sequence.Length < wordSize` ⇒ score = 0 | Guard in `kmer_ratio` [2] |

### 2.5 Comparison with Related Methods

| Aspect | CPAT hexamer score | Fickett TESTCODE [4] |
|--------|--------------------|----------------------|
| Signal | In-frame hexamer log-likelihood vs trained tables | Position + composition asymmetry, no training tables |
| Requires training data | Yes (coding & non-coding tables) | No |
| Output sign | + coding / − non-coding | Probability-like score |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequence` | `string` | required | DNA sequence to score | non-null; case-insensitive (upper-cased internally) |
| `codingHexamerFrequencies` | `IReadOnlyDictionary<string,double>` | required | In-frame hexamer table from coding (CDS) training set | non-null; uppercase A/C/G/T keys; non-negative values |
| `noncodingHexamerFrequencies` | `IReadOnlyDictionary<string,double>` | required | In-frame hexamer table from non-coding (background) set | non-null; same units as coding table |
| `wordSize` | `int` | 6 | Hexamer length [1][2] | > 0 |
| `stepSize` | `int` | 3 | In-frame window step [2] | > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `double` | Mean per-hexamer log-likelihood ratio. > 0 ⇒ coding, < 0 ⇒ non-coding [1]. 0 when no hexamer is scorable. |

### 3.3 Preconditions and Validation

Null `sequence` or either table → `ArgumentNullException`. `wordSize ≤ 0` or `stepSize ≤ 0` → `ArgumentOutOfRangeException`. `sequence.Length < wordSize` (including empty) → returns 0 [2]. Input is upper-cased (`ToUpperInvariant`); only frame 0 is scored. Keys absent from either table are skipped; no exception for hexamers containing non-ACGT characters (they simply are not in the tables) [2].

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; if `sequence.Length < wordSize`, return 0 [2].
2. Upper-case the sequence.
3. For `i = 0, stepSize, 2·stepSize, …` while `i + wordSize ≤ length`, take the hexamer `seq[i, wordSize]` [2].
4. If the hexamer is missing from either table, skip it [2].
5. Otherwise add its contribution (`ln(F/F')`, `+1`, `−1`, or 0) and increment the scored-hexamer count [2].
6. Return `sum / count`, or 0 if `count = 0`.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

The per-hexamer contribution rule (§2.2) is the only decision table; constants `wordSize = 6` and `stepSize = 3` originate from CPAT [1][2]. The frequency tables themselves are user-supplied training data, not hard-coded.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Score one sequence | O(n) | O(1) | n = sequence length; `(n − wordSize)/stepSize + 1` hexamer lookups, each O(1) hash lookup over fixed wordSize |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomeAnnotator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs)

- `GenomeAnnotator.CalculateCodingPotential(string, IReadOnlyDictionary<string,double>, IReadOnlyDictionary<string,double>, int, int)`: computes the CPAT hexamer usage-bias mean log-likelihood (frame 0).

### 5.2 Current Behavior

Implements `FrameKmer.kmer_ratio` frame-0 path verbatim [2]. The `coding == 0 && noncoding == 0` case (a hexamer present in both tables with zero value) is `continue`d — skipped and **not counted** toward the scored-hexamer denominator, matching the reference branch `elif coding[k]==0 and noncoding[k]==0: continue` (verified against canonical CPAT `liguowang/cpat` `src/cpmodule/FrameKmer.py` and lncScore, 2026-06-15). For the degenerate "no scorable hexamer" case the implementation returns 0 (the reference returns −1 there); see §5.4. No substring/search reuse is relevant: scoring is a single linear in-frame scan with O(1) table lookups, so the repository suffix tree does not apply (no occurrence enumeration or multi-query matching).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- In-frame hexamer extraction: start 0, step 3, word size 6, full-length words only [2].
- Per-hexamer contribution `ln(F/F')` (natural log), `+1` (coding-only), `−1` (non-coding-only), skip-if-missing [2].
- Mean over the number of scored hexamers [2].
- Sign convention positive=coding / negative=non-coding [1].

**Intentionally simplified:**

- No-scorable-hexamer result: returns 0 instead of the reference −1 sentinel; **consequence:** for inputs with zero scorable hexamers the score reads 0 ("no information") rather than −1. Both are non-scores.

**Not implemented:**

- The full CPAT logistic-regression classifier (ORF size, ORF coverage, Fickett TESTCODE, hexamer score combined) [1]; **users should rely on:** this method as the hexamer feature only — combine externally for a full coding/non-coding call.
- Fickett TESTCODE statistic [4]; **users should rely on:** EMBOSS `tcode` or a dedicated implementation — no current in-repo alternative.
- Built-in organism frequency tables; **users should rely on:** supplying tables (e.g. generated by CPAT `make_hexamer_tab.py`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | No-scorable-hexamer returns 0, reference returns −1 | Deviation | Only affects inputs with zero scored hexamers; both are non-scores | accepted | ASM/ASSUMPTION-1 |
| 2 | Tables must share units (counts vs proportions) | Assumption | Score shifts by `ln(Σ/Σ)` if mixed | accepted | ASM-02 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| `sequence.Length < wordSize` (incl. empty) | 0 | Guard [2] |
| Hexamer missing from one/both tables | skipped, not counted | `has_key` guard [2] |
| Coding-only hexamer | +1 contribution | `kmer_ratio` branch [2] |
| Non-coding-only hexamer | −1 contribution | `kmer_ratio` branch [2] |
| Both tables value 0 for an in-both hexamer | skipped (`continue`), not counted | branch `coding==0 && noncoding==0: continue` [2] |
| No scorable hexamer at all | 0 | implementation choice (§5.4) |
| null sequence / table | `ArgumentNullException` | validation |
| wordSize ≤ 0 / stepSize ≤ 0 | `ArgumentOutOfRangeException` | validation |

### 6.2 Limitations

Requires user-supplied, organism-appropriate hexamer tables (no built-in defaults). Scores only frame 0 of the supplied sequence — the caller must supply the ORF in frame. The hexamer score alone is one feature; a complete coding/non-coding decision in CPAT also uses ORF size, ORF coverage and Fickett TESTCODE [1]. Mixing table units silently shifts the score (ASM-02).

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var coding = new Dictionary<string, double> { ["ATGAAA"] = 8, ["AAACCC"] = 2 };
var noncoding = new Dictionary<string, double> { ["ATGAAA"] = 2, ["AAACCC"] = 4 };
double score = GenomeAnnotator.CalculateCodingPotential("ATGAAACCC", coding, noncoding);
// score = (ln(8/2) + ln(2/4)) / 2 = (1.3862943611 + (-0.6931471806)) / 2 = 0.3465735903
```

**Numerical walk-through:** sequence `ATGAAACCC` (length 9), wordSize 6, stepSize 3. In-frame hexamers: `ATGAAA` (i=0), `AAACCC` (i=3); i=6 yields `CCC` (length 3 < 6) → skipped. `ATGAAA`: ln(8/2)=1.38629436; `AAACCC`: ln(2/4)=−0.69314718. Sum 0.69314718, count 2 ⇒ 0.34657359027997264 [2].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomeAnnotator_CalculateCodingPotential_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Annotation/GenomeAnnotator_CalculateCodingPotential_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [ANNOT-CODING-001-Evidence.md](../../../docs/Evidence/ANNOT-CODING-001-Evidence.md)

## 8. References

1. Wang L, Park HJ, Dasari S, Wang S, Kocher J-P, Li W. 2013. CPAT: Coding-Potential Assessment Tool using an alignment-free logistic regression model. Nucleic Acids Research 41(6):e74. https://doi.org/10.1093/nar/gkt006
2. CPAT / lncScore reference implementation, `cpmodule/FrameKmer.py` (`kmer_ratio`, `word_generator`, `kmer_freq_file`). https://github.com/WGLab/lncScore/blob/master/tools/cpmodule/FrameKmer.py
3. Fickett JW, Tung C-S. 1992. Assessment of protein coding measures. Nucleic Acids Research 20(24):6441–6450. https://doi.org/10.1093/nar/20.24.6441
4. Fickett JW. 1982. Recognition of protein coding regions in DNA sequences. Nucleic Acids Research 10(17):5303–5318. https://doi.org/10.1093/nar/10.17.5303
