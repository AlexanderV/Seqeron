# Pairwise Alignment Statistics and Formatting

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-STATS-001 |
| Related Projects | Seqeron.Genomics.Alignment, Seqeron.Genomics.Infrastructure |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Given a completed pairwise alignment (two equal-length strings over the alphabet plus the gap symbol `-`), this unit computes the standard summary statistics reported by global/local aligners — **Identity**, **Similarity**, and **Gaps** (each as a count and a percentage of the alignment length) — and renders a human-readable three-line alignment with a markup/consensus line. The statistics follow the EMBOSS `needle`/`water` definitions [2] (corroborated by NCBI BLAST [4] and the `pseqsid` reference implementation [5]); the markup follows the EMBOSS `srspair` legend [3]. The computation is exact and deterministic, O(n) in the alignment length.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

After a dynamic-programming aligner produces two gapped sequences of equal length L, biologists summarise the relationship with three percentages over the L aligned columns. The denominator L **includes gap columns**, which is the EMBOSS `needle` convention [2] and the convention corroborated for percent identity generally [secondary corroboration]. These summaries let users judge relatedness independently of the raw score.

### 2.2 Core Model

For an alignment of two rows `a`, `b` of equal length L, each column `i` is classified:

- **gap** if `a[i] = '-'` or `b[i] = '-'`;
- **identical** if not a gap and `a[i] = b[i]`;
- **mismatch** otherwise (not a gap and `a[i] ≠ b[i]`).

A non-gap, non-identical column is additionally **similar** iff its substitution score `s(a[i], b[i])` is **positive** [2][4][5]. Identical columns are always similar (their score is positive). Let `M` = identical count, `X` = mismatch count, `G` = gap count, and `Sim⁺` = number of mismatch columns whose substitution score is positive. Then [2]:

- Identity% = `M / L × 100`
- Similarity% = `(M + Sim⁺) / L × 100`
- Gaps% = `G / L × 100`

with `M + X + G = L`. The published EMBOSS worked example (HBA_HUMAN vs HBB_HUMAN, EBLOSUM62) gives L = 149, Identity 65/149 = 43.6%, Similarity 90/149 = 60.4%, Gaps 9/149 = 6.0% [2], confirming both the formula and the gap-inclusive denominator.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `Matches + Mismatches + Gaps = AlignmentLength` | The three classes partition every column [2][4] |
| INV-02 | `Identity% ≤ Similarity% ≤ 100%` | The similar set ⊇ the identical set [2][5] |
| INV-03 | `Identity% + (Similarity−Identity)% + mismatch-only% + Gap% = 100%` | Length-denominator partition [2] |
| INV-04 | Markup line length equals the displayed sequence-block length; characters ∈ {`|`, `:`, space} | srspair legend [3] |

### 2.5 Comparison with Related Methods (Optional)

| Aspect | EMBOSS-style (this unit) | "Non-gap fraction" (naive) |
|--------|--------------------------|-----------------------------|
| Similarity definition | identical + positive-score columns | identical + all mismatch columns |
| Requires substitution model | yes | no |
| Conforms to EMBOSS/BLAST | yes [2][4] | no |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `alignment` | `AlignmentResult` | required | Completed alignment (two equal-length gapped rows) | non-null |
| `scoring` | `ScoringMatrix?` | `SimpleDna` | Model deciding whether a mismatch is "similar" (positive score) | — |
| `lineWidth` | `int` | 60 | Residues per display block (FormatAlignment only) | > 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `Matches` | int | Identical columns |
| `Mismatches` | int | Non-gap, non-identical columns |
| `Gaps` | int | Columns with a gap on either row |
| `AlignmentLength` | int | L = total columns (includes gaps) |
| `Identity` | double | `Matches/L × 100` |
| `Similarity` | double | `(Matches + positive-score mismatches)/L × 100` |
| `GapPercent` | double | `Gaps/L × 100` |
| (FormatAlignment) | string | Three-line blocks: row1, markup (`|`/`:`/space), row2 |

### 3.3 Preconditions and Validation

`alignment` must be non-null (`ArgumentNullException` otherwise). An empty alignment (`AlignedSequence1` null/empty) returns `AlignmentStatistics.Empty` / `""` respectively. `lineWidth` must be positive (`ArgumentOutOfRangeException` otherwise). The two rows are assumed equal length (an aligner postcondition); the gap symbol is `-`. Comparison is case-sensitive on the characters as stored (the aligners upper-case input before alignment).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate inputs; return Empty/`""` for an empty alignment.
2. Walk the L columns once; classify each as gap, identical, or mismatch.
3. For a mismatch, increment the similar counter iff the substitution score is positive (for the scalar DNA model: `scoring.Mismatch > 0`).
4. Compute the three percentages with denominator L (gap-inclusive).
5. (FormatAlignment) emit `|` for identical, `:` for similar mismatch, space for gap or non-positive mismatch, wrapping every `lineWidth` columns.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures (Optional)

srspair markup legend [3]: `|` = identical, `:` = similar (positive substitution score), space = gap or mismatch. EMBOSS additionally defines `.` for small positive scores; the scalar `ScoringMatrix` used here has no graded positive tier, so positive non-identical columns render `:` and the `.` tier is unreachable (rendering-only; does not affect any counted statistic).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateStatistics | O(L) | O(1) | single pass over L columns |
| FormatAlignment | O(L) | O(L) | builds the output string |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.CalculateStatistics(AlignmentResult, ScoringMatrix?)`: returns `AlignmentStatistics` (Identity/Similarity/Gaps).
- `SequenceAligner.FormatAlignment(AlignmentResult, int, ScoringMatrix?)`: returns the three-line srspair-style display.

### 5.2 Current Behavior

Similarity is parameterised by the `ScoringMatrix`: a mismatch column counts as similar only when its substitution score is positive. For all DNA models exposed by this class (`SimpleDna`, `BlastDna`, `HighIdentityDna`), `Mismatch < 0`, so no mismatch is similar and **Similarity equals Identity**; a caller supplying a model with `Mismatch > 0` will see Similarity exceed Identity. The denominator is the full alignment length including gap columns. This unit performs no search/matching, so the repository suffix tree is **not applicable** (N/A) — it is a single linear pass over an already-computed alignment.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Identity = identical/L × 100 with gap-inclusive L [2].
- Similarity = (identical + positive-score)/L × 100, positive-substitution-score rule [2][4][5].
- Gaps = gap-columns/L × 100 [2].
- srspair markup legend `|`/`:`/space [3].

**Intentionally simplified:**

- srspair `.` (small-positive-score) tier: not emitted; **consequence:** positive non-identical columns are all marked `:` rather than splitting into `:`/`.`, because the scalar `ScoringMatrix` has no graded positive scores. Counted statistics are unaffected.

**Not implemented:**

- Protein substitution matrices (BLOSUM/PAM): not provided by `ScoringMatrix`; **users should rely on:** supplying a model whose `Mismatch` reflects positive substitutability, or an external protein aligner for true BLOSUM-based similarity.

### 5.4 Deviations and Assumptions (Optional)

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Prior Similarity = `(matches+mismatches)/L` (non-gap fraction) | Deviation | Over-reported similarity vs EMBOSS/BLAST | fixed | Corrected to positive-score rule in ALIGN-STATS-001 |
| 2 | srspair `.` tier not emitted | Assumption | Display only | accepted | See 5.3 |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty alignment | `AlignmentStatistics.Empty` / `""` | Denominator undefined; documented empty handling |
| Null alignment | `ArgumentNullException` | Validation contract |
| `lineWidth ≤ 0` | `ArgumentOutOfRangeException` | Validation contract |
| All-gap column run | counted as gaps; Identity 0% if no identical column | Gap classification [2][4] |
| Perfect identity | Identity = Similarity = 100%, Gap 0% | Partition invariant |

### 6.2 Limitations

Only the scalar DNA `Match`/`Mismatch` model is available, so similarity beyond identity requires a caller-supplied `Mismatch > 0` model; there is no built-in BLOSUM/PAM matrix. Statistics assume the two rows are equal length (aligner postcondition) and use `-` as the sole gap symbol. The `.` (small positive score) markup tier is not produced.

## 7. Examples and Related Material (Optional)

### 7.1 Worked Example

**API usage example:**

```csharp
var aln = SequenceAligner.GlobalAlign("ACGTACGT", "ACCTAACGT", SequenceAligner.SimpleDna);
var stats = SequenceAligner.CalculateStatistics(aln); // SimpleDna ⇒ Similarity == Identity
string view = SequenceAligner.FormatAlignment(aln, lineWidth: 60);
```

**Numerical / biological walk-through (optional):**

EMBOSS HBA vs HBB [2]: L = 149, M = 65 ⇒ Identity 65/149 = 43.624…% ≈ 43.6%; similar columns = 90 ⇒ Similarity 90/149 = 60.402…% ≈ 60.4%; G = 9 ⇒ Gaps 9/149 = 6.040…% ≈ 6.0%.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAligner_CalculateStatistics_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_CalculateStatistics_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [ALIGN-STATS-001-Evidence.md](../../../docs/Evidence/ALIGN-STATS-001-Evidence.md)
- Related algorithms: [Global_Alignment_Needleman_Wunsch](../Alignment/Global_Alignment_Needleman_Wunsch.md)

## 8. References

1. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2
2. EMBOSS needle application documentation, release 6.6. EMBOSS. https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
3. EMBOSS Alignment Formats (markup legend). EMBOSS. https://emboss.sourceforge.net/docs/themes/AlignFormats.html
4. NCBI. BLAST QuickStart — Comparative Genomics. NCBI Bookshelf NBK1734. https://www.ncbi.nlm.nih.gov/books/NBK1734/
5. Pérez-Mejías A (amaurypm). pseqsid — pairwise sequence identity/similarity. GitHub. https://github.com/amaurypm/pseqsid
