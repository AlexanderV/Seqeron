# Sequence Similarity (k-mer Jaccard Index)

| Field | Value |
|-------|-------|
| Algorithm Group | Analysis |
| Test Unit ID | GENOMIC-SIMILARITY-001 |
| Related Projects | Seqeron.Genomics.Analysis, Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-06-14 |

## 1. Overview

`GenomicAnalyzer.CalculateSimilarity` measures the similarity of two DNA sequences using the
**Jaccard index of their k-mer sets**. Each sequence is reduced to the set of its distinct
length-k substrings (k-mers); the Jaccard index is the fraction of k-mers shared between the
two sets out of all distinct k-mers in either set [1][2]. The result is exact (not heuristic or
probabilistic — there is no MinHash sketching here) and is reported as a percentage in [0, 100].
It is an alignment-free comparison suitable for quick resemblance estimates between sequences.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Alignment-free sequence comparison summarizes a sequence by its k-mer composition and compares
those compositions instead of aligning residues. The Jaccard index is the canonical set-resemblance
measure used for this purpose; Mash popularized k-mer Jaccard (via MinHash sketches) for genome and
metagenome distance estimation [2].

### 2.2 Core Model

For two finite sets A and B, the **Jaccard index** is [1]:

```
J(A,B) = |A ∩ B| / |A ∪ B| = |A ∩ B| / (|A| + |B| − |A ∩ B|)
```

with `0 ≤ J(A,B) ≤ 1` [1]. Here A and B are the **sets of distinct k-mers** of the two sequences.
Applied to k-mer sets, J is "the fraction of shared k-mers out of all distinct k-mers in A and B" [2].
This method reports `J × 100` as a percentage.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `0 ≤ result ≤ 100` | `0 ≤ J ≤ 1` [1], scaled ×100 |
| INV-02 | Identical non-empty sequences → 100.0 | A = B ⇒ A∩B = A∪B ⇒ J = 1 [1] |
| INV-03 | Disjoint k-mer sets → 0.0 | A∩B = ∅, A∪B ≠ ∅ ⇒ J = 0 [1] |
| INV-04 | Symmetric: result(a,b,k) = result(b,a,k) | ∩ and ∪ are commutative [1] |
| INV-05 | k-mers compared as sets (within-sequence repeats counted once) | Jaccard over distinct k-mers [2] |

### 2.5 Comparison with Related Methods

| Aspect | k-mer Jaccard (this) | Alignment identity |
|--------|----------------------|--------------------|
| Basis | set resemblance of k-mer composition | residue-by-residue over an alignment |
| Cost | O(n+m) | O(n·m) (dynamic programming) |
| Sensitivity to order | low (composition only) | high (positional) |
| Sketching | none here (exact sets); Mash adds MinHash [2] | n/a |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| sequence1 | DnaSequence | required | first sequence | non-null; alphabet A/C/G/T (enforced by DnaSequence) |
| sequence2 | DnaSequence | required | second sequence | non-null; alphabet A/C/G/T |
| kmerSize | int | 5 | k-mer length | ≥ 1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | double | k-mer Jaccard index ×100, a percentage in [0, 100] |

### 3.3 Preconditions and Validation

- `sequence1` / `sequence2` null → `ArgumentNullException`.
- `kmerSize < 1` → `ArgumentOutOfRangeException`.
- Sequences are uppercase A/C/G/T (normalized and validated by `DnaSequence`); comparison is therefore case-insensitive at the `DnaSequence` boundary.
- Empty sequences (or sequences shorter than k) produce empty k-mer sets; if both are empty the union is empty and the method returns `0.0` (see §5.4).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate arguments (non-null sequences; `kmerSize ≥ 1`).
2. Build the distinct-k-mer set of each sequence (sliding window, `HashSet<string>`).
3. Compute `|A ∩ B|` (count of A's k-mers present in B) and `|A ∪ B| = |A| + |B| − |A ∩ B|`.
4. Return `0` if the union is empty; otherwise `|A∩B| / |A∪B| × 100`.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateSimilarity | O(n + m) | O(n + m) | n,m = sequence lengths; k-mer extraction + hash-set intersection, each k-mer O(k) treated as O(1) for fixed k |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [GenomicAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs)

- `GenomicAnalyzer.CalculateSimilarity(DnaSequence, DnaSequence, int)`: computes the k-mer Jaccard index ×100.
- `GenomicAnalyzer.GetKmers(string, int)` (private): builds the distinct-k-mer `HashSet<string>`.

### 5.2 Current Behavior

- k-mers are stored in a `HashSet<string>`, so duplicate k-mers within a sequence collapse to one element — exactly the distinct-set semantics required by the Jaccard definition [1][2].
- Intersection is computed as `kmers1.Count(kmers2.Contains)` and union as `|A|+|B|−|A∩B|`, avoiding a second pass over the data.
- **Suffix tree evaluated but not used.** This unit computes a set intersection/union of distinct k-mers, not an occurrence-search ("find all positions of X in Y") problem. A `HashSet<string>` gives O(n+m) construction and O(1) membership, which is optimal for set resemblance; the repository `SuffixTree` would add O(n) construction cost without providing the set-overlap query this metric needs. The suffix tree is therefore correctly not applied here.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- `J(A,B) = |A∩B| / |A∪B|` over the distinct-k-mer sets [1].
- k-mer Jaccard as the fraction of shared k-mers out of all distinct k-mers [2].
- Range and extremes: J=1 for identical sets, J=0 for disjoint sets [1].

**Intentionally simplified:**

- Percentage reporting: result is `J × 100` rather than the raw [0,1] coefficient; **consequence:** values are 0–100, not 0–1. Relative ordering is unchanged.

**Not implemented:**

- MinHash sketching / approximate Jaccard estimation (Mash's `j(A_s,B_s)=|A_s∩B_s|/s`); **users should rely on:** this method computes the exact Jaccard over full k-mer sets, so no sketch approximation is needed for the supported sequence sizes [2].
- Jaccard distance `1 − J`; **users should rely on:** derive it as `1 − result/100` if needed.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Empty-union returns 0.0 | Assumption | Affects both-empty / both-shorter-than-k inputs | accepted | Jaccard undefined for empty union [1]; 0 chosen as "no similarity" (ASM-1 in TestSpec) |
| 2 | ×100 scaling | Assumption | Output in [0,100] not [0,1] | accepted | presentation convention only |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical non-empty sequences | 100.0 | J = 1 [1] |
| Disjoint k-mer sets | 0.0 | J = 0 [1] |
| Both sequences empty | 0.0 | empty union; undefined → impl returns 0 [1] |
| Both sequences shorter than k | 0.0 | both k-mer sets empty → empty union |
| One sequence empty, other non-empty | 0.0 | empty intersection over non-empty union |
| Within-sequence repeated k-mers | counted once | distinct-set semantics [2] |
| Null sequence | ArgumentNullException | input validation |
| kmerSize < 1 | ArgumentOutOfRangeException | input validation |

### 6.2 Limitations

- Composition-only: insensitive to k-mer order/position, so two sequences with the same k-mer composition but different arrangement score identically. For positional/edit comparison use alignment.
- Exact-match k-mers only — no mismatch tolerance; a single substitution alters up to k k-mers.
- Choice of k (default 5) controls resolution: larger k is more specific but more sensitive to errors; the formula is unchanged.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var a = new DnaSequence("ACGTACGT");
var b = new DnaSequence("ACGTACGA");
double pct = GenomicAnalyzer.CalculateSimilarity(a, b, kmerSize: 3); // 80.0
```

**Numerical walk-through (k=3, ACGTACGT vs ACGTACGA):**

- A = distinct 3-mers of ACGTACGT = {ACG, CGT, GTA, TAC} (|A| = 4; ACG and CGT repeat but collapse).
- B = distinct 3-mers of ACGTACGA = {ACG, CGT, GTA, TAC, CGA} (|B| = 5).
- A ∩ B = {ACG, CGT, GTA, TAC} → |A∩B| = 4.
- A ∪ B = {ACG, CGT, GTA, TAC, CGA} → |A∪B| = 5.
- J = 4/5 = 0.8 → ×100 = **80.0**.

### 7.3 Related Tests, Evidence, or Documents

- Tests: [GenomicAnalyzer_CalculateSimilarity_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/GenomicAnalyzer_CalculateSimilarity_Tests.cs) — covers INV-01..INV-05
- Evidence: [GENOMIC-SIMILARITY-001-Evidence.md](../../../docs/Evidence/GENOMIC-SIMILARITY-001-Evidence.md)

## 8. References

1. Jaccard, Paul. 1901. Étude comparative de la distribution florale dans une portion des Alpes et des Jura. Bulletin de la Société vaudoise des sciences naturelles 37(142):547–579. https://en.wikipedia.org/wiki/Jaccard_index
2. Ondov BD, Treangen TJ, Melsted P, Mallonee AB, Bergman NH, Koren S, Phillippy AM. 2016. Mash: fast genome and metagenome distance estimation using MinHash. Genome Biology 17:132. https://doi.org/10.1186/s13059-016-0997-x
3. Mash documentation — Distance Estimation. marbl/Mash. https://github.com/marbl/Mash/blob/master/doc/sphinx/distances.rst
