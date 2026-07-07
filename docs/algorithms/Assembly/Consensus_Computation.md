# Consensus Computation

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-CONSENSUS-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Consensus computation collapses a set of pre-aligned reads into a single sequence by choosing, for each alignment column, the most-supported residue. It is the "Consensus" step of the Overlap-Layout-Consensus assembly paradigm and is also used to summarise any multiple sequence alignment. The rule here is specification-driven, following Biopython's `dumb_consensus`: a residue is committed to the consensus only when it is the unique most frequent residue in the column and its frequency among non-gap residues meets a plurality threshold; otherwise an ambiguity symbol is emitted [1]. The method is exact and deterministic.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A consensus sequence is "the calculated sequence of most frequent residues, either nucleotide or amino acid, found at each position in a sequence alignment" [3]. Positions where no single residue dominates are written with an IUPAC degenerate symbol — for DNA, `N` denotes "any base" [3]. Tools such as EMBOSS `cons` formalise this with a *plurality* cut-off: support below the cut-off means "there is no consensus" at that column [2].

### 2.2 Core Model

For an alignment of `R` reads with column index `n` over `0 ≤ n < L` (where `L` is the longest read), let the per-column tally count only non-gap residues (characters other than `-` and `.`), and let `num_atoms` be the number of such residues. Let `max_size` be the maximum count and `max_atoms` the set of residues achieving it. Biopython's decision rule [1] is, verbatim:

```python
elif (len(max_atoms) == 1) and ((float(max_size) / float(num_atoms)) >= threshold):
    consensus += max_atoms[0]
else:
    consensus += ambiguous
```

So a residue `x` is committed iff it is the **unique** maximum-count residue **and** `max_size / num_atoms ≥ threshold`; otherwise the ambiguity symbol is emitted [1]. Gap characters `-` and `.` are excluded from the tally [1]. The consensus length equals the full alignment length `con_len = get_alignment_length()` [1].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Output length equals the longest input read. | Consensus is built over `con_len = get_alignment_length()` [1] |
| INV-02 | A column whose single most-frequent residue has frequency `≥ threshold` (among non-gap residues) emits that residue; otherwise ambiguous. | Biopython decision rule [1]; EMBOSS plurality cut-off [2] |
| INV-03 | A column with ≥2 residues tied for the maximum count emits the ambiguous symbol. | `len(max_atoms) == 1` guard [1] |
| INV-04 | Gap characters `-`/`.` never appear in the output and never count toward `num_atoms`. | Tally skips `-` and `.` [1] |
| INV-05 | An all-gap (or empty) column emits the ambiguous symbol with no division by zero. | `len(max_atoms)==1` is false, the division short-circuits [1] |

### 2.5 Comparison with Related Methods

| Aspect | This method (`dumb_consensus`-style) | Naive `MaxBy` majority |
|--------|--------------------------------------|------------------------|
| Tie handling | Ambiguous symbol (deterministic) [1] | Arbitrary winner (nondeterministic) |
| Plurality cut-off | Yes (`≥ threshold`) [1][2] | None |
| All-gap column | Ambiguous, no error [1] | Undefined |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `alignedReads` | `IReadOnlyList<string>` | required | Pre-aligned reads (columns must correspond); ragged lengths allowed | Not null |
| `threshold` | `double` | `0.5` | Minimum column frequency (max count / non-gap count) to commit a residue | 0–1; pass 0.7 to reproduce Biopython default |
| `ambiguous` | `char` | `'N'` | Symbol emitted when a column has no committed residue | Any char |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `string` | Consensus, length equal to the longest read; committed residues uppercased, ambiguous columns set to `ambiguous` |

### 3.3 Preconditions and Validation

`alignedReads` null → `ArgumentNullException`. Empty list → empty string. Reads need not be equal length; columns past a shorter read contribute nothing from that read (INV-01). Residues are normalized with `char.ToUpperInvariant`; `-` and `.` are treated as gaps and skipped. The comparison is `≥ threshold` (inclusive), per the source [1].

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; empty list → `""`.
2. `L` = length of the longest read.
3. For each column `0 ≤ n < L`: tally non-gap residues (skip `-`, `.`), counting `num_atoms`.
4. Find `max_size` and whether a unique residue achieves it.
5. Commit that residue iff it is unique **and** `max_size / num_atoms ≥ threshold`; else emit `ambiguous`.
6. Concatenate column results.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Gap set: `{ '-', '.' }` excluded from tally [1].
- Commit predicate: `uniqueMax ∧ numAtoms>0 ∧ (maxSize/numAtoms ≥ threshold)` [1].
- Defaults: `threshold = 0.5` (EMBOSS simple-majority plurality [2]); `ambiguous = 'N'` (IUPAC any-base [3]). Biopython's documented defaults (`0.7`, `'X'`) are reachable via parameters [1].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `ComputeConsensus` | O(L × R) | O(L + a) | L = longest read, R = read count, a = alphabet size per column (small dictionary) |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.ComputeConsensus(IReadOnlyList<string>, double, char)`: column-wise majority/threshold consensus.

### 5.2 Current Behavior

Single linear pass per column tallying into a `Dictionary<char,int>`. The maximum-count residue and tie detection are computed in one loop (tracking `maxSize` and the count of residues sharing it), so a tie sets `maxCount > 1` and forces the ambiguous branch. Residues are uppercased; gaps skipped. No search/matching is performed (this is a column reduction, not occurrence finding), so the repository suffix tree is **not** applicable here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Commit predicate `unique-max ∧ frequency ≥ threshold`, else ambiguous [1].
- Gap characters `-`/`.` excluded from the tally [1].
- Tie for the maximum count → ambiguous symbol [1].
- Consensus length = longest read / alignment length [1].
- All-gap / empty column → ambiguous, no division by zero [1].

**Intentionally simplified:**

- Default `threshold = 0.5` instead of Biopython's documented `0.7`; **consequence:** the default behaves as a simple-majority vote; callers reproduce Biopython exactly with `threshold: 0.7`.
- Default ambiguous symbol `'N'` (DNA IUPAC) instead of Biopython's `'X'`; **consequence:** DNA-appropriate output by default; callers pass `ambiguous: 'X'` for protein alignments.

**Not implemented:**

- Sequence weighting and scoring-matrix-based scoring (EMBOSS `cons` style); **users should rely on:** dedicated MSA scoring tools; this method implements the unweighted plurality/majority rule only.
- `require_multiple` Biopython flag (force ambiguous when only one read covers a column); **users should rely on:** filtering low-coverage columns upstream.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Default threshold 0.5 vs Biopython 0.7 | Assumption | Default vote strictness; all values reachable | accepted | ASM/parameter; M9 verifies 0.7 |
| 2 | Default ambiguous 'N' vs 'X' | Assumption | Emitted symbol for no-consensus columns | accepted | DNA IUPAC; configurable |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty read list | `""` | No columns |
| Null read list | `ArgumentNullException` | Repository contract |
| All-gap column | `ambiguous` | INV-05 [1] |
| Tie for max count | `ambiguous` | INV-03 [1] |
| Sub-threshold majority | `ambiguous` | INV-02 [1] |
| Ragged reads | length = longest read | INV-01 [1] |

### 6.2 Limitations

Requires reads that are already aligned column-for-column; it does not perform alignment. It is an unweighted plurality rule (no per-sequence weights, no scoring matrix). It does not emit blended IUPAC degeneracy codes for ties (a single ambiguity symbol is used), and it does not implement `require_multiple`.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reads = new[] { "ACGT", "ACGT", "TCGT" };
string c = SequenceAssembler.ComputeConsensus(reads, threshold: 0.5); // "ACGT" (col0: A=2/3 ≥ 0.5)
string d = SequenceAssembler.ComputeConsensus(reads, threshold: 0.7); // "NCGT" (col0: A=2/3 < 0.7)
```

**Numerical walk-through:** column 0 = {A,A,T}, `num_atoms=3`, `max_size=2` (A unique). At threshold 0.5: 2/3≈0.667 ≥ 0.5 → `A`. At threshold 0.7: 0.667 < 0.7 → `N` [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_ComputeConsensus_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAssembler_ComputeConsensus_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [ASSEMBLY-CONSENSUS-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-CONSENSUS-001-Evidence.md)
- Related algorithms: [Coverage_Calculation](../Assembly/Coverage_Calculation.md)

## 8. References

1. Cock PJA, Antao T, Chang JT, et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 — `dumb_consensus` source (v1.79): https://raw.githubusercontent.com/biopython/biopython/biopython-179/Bio/Align/AlignInfo.py
2. Rice P, Longden I, Bleasby A. 2000. EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2 — `cons` documentation: https://emboss.sourceforge.net/apps/cvs/emboss/apps/cons.html
3. Wikipedia. Consensus sequence. https://en.wikipedia.org/wiki/Consensus_sequence (definition; IUPAC degenerate notation; cites Schneider & Stephens 1990, https://doi.org/10.1093/nar/18.20.6097).
