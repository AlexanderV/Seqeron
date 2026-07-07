# Reversal Distance (Breakpoint Lower Bound)

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-REVERSAL-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-14 |

## 1. Overview

Reversal distance is the minimum number of reversals (inversions of contiguous blocks) that transform one gene order into another — a classical measure of genome rearrangement. Computing the *exact* reversal distance is non-trivial (signed permutations require the Hannenhalli–Pevzner cycle/hurdle analysis [1][2]). This implementation returns the **unsigned breakpoint lower bound** ⌈b/2⌉, where b is the number of breakpoints between the two gene orders [1][2]. It is exact for the trivial cases (identity ⇒ 0) and is a guaranteed lower bound otherwise; it is fast (O(n)) and well-defined.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A genome's gene order can be modeled as a permutation of marker ids; an evolutionary reversal flips a contiguous segment. Comparing two genomes by the number of reversals separating them estimates rearrangement distance [1].

### 2.2 Core Model

Given a permutation π, form the **extended permutation** (π_0, π_1, …, π_n, π_{n+1}) = (0, π_1, …, π_n, n+1) [2]. A pair of consecutive elements (π_i, π_{i+1}) is a **breakpoint** iff π_{i+1} ≠ π_i + 1; the identity permutation is the only one with 0 breakpoints [1, §2]. For the **unsigned** model used here, a pair is a breakpoint iff the two values are not consecutive integers, i.e. |π_{i+1} − π_i| ≠ 1 [3].

Let b(π) be the number of breakpoints. A single reversal cuts exactly two adjacencies, so it removes at most two breakpoints: b(α) − b(αρ) ≤ 2 [2]. Summing over a sorting sequence of t reversals, b(α) ≤ 2t, and since the true distance d(α) ≥ t,

> d(α) ≥ b(α) / 2 [2].

The smallest integer satisfying this is ⌈b(α)/2⌉, which is the returned value.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | d(π, π) = 0 | identity is the unique permutation with 0 breakpoints [1, §2] |
| INV-02 | result ≥ 0 | ⌈b/2⌉ with b ≥ 0 |
| INV-03 | d(α, β) = d(β, α) | reversal distance is symmetric [2] |
| INV-04 | result = ⌈b/2⌉ over the extended relative permutation | breakpoint definition [1][3] |
| INV-05 | result is a lower bound on the true reversal distance | b(α) ≤ 2t ⇒ d ≥ b/2 [2] |

### 2.5 Comparison with Related Methods

| Aspect | This method (breakpoint bound) | Hannenhalli–Pevzner |
|--------|--------------------------------|----------------------|
| Output | lower bound ⌈b/2⌉ | exact reversal distance |
| Model | unsigned breakpoints | signed breakpoint graph (cycles + hurdles) |
| Complexity | O(n) | O(n) distance via b − c + h (more involved) [1][2] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| permutation1 | `IReadOnlyList<int>` | required | source gene order (distinct marker ids) | same length as permutation2; ids in permutation1 must appear in permutation2 |
| permutation2 | `IReadOnlyList<int>` | required | target gene order over the same marker set | same length as permutation1 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `int` | lower bound ⌈b/2⌉ ≥ 0 on the reversal distance |

### 3.3 Preconditions and Validation

Inputs are treated as **unsigned** permutations (no strand sign). Marker ids may be any distinct integers; they are relabelled to the relative permutation (target ⇒ identity 0..n−1). If the two lists differ in length, an `ArgumentException` is thrown. Empty or single-element inputs return 0 (no internal adjacency, hence no breakpoint). Indexing is 0-based internally; the result is a count.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate equal lengths; return 0 if n ≤ 1.
2. Build the relative permutation by mapping each marker of permutation1 to its index in permutation2 (target ⇒ identity).
3. Count breakpoints of the extended relative permutation (−1, relative…, n): the left boundary (relative[0] ≠ 0), each internal pair (|Δ| ≠ 1), and the right boundary (relative[n−1] ≠ n−1).
4. Return ⌈b/2⌉ = (b + 1) / 2.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| CalculateReversalDistance | O(n) | O(n) | one hash map for the position index; single pass over n pairs |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.CalculateReversalDistance(permutation1, permutation2)`: returns the unsigned breakpoint lower bound ⌈b/2⌉.

### 5.2 Current Behavior

The two end boundaries are treated as breakpoints against the sentinels (left = −1 ⇒ breakpoint iff relative[0] ≠ 0; right = n ⇒ breakpoint iff relative[n−1] ≠ n−1), matching the extended-permutation definition [2]. No substring search is performed, so the repository suffix tree is **not applicable** (this is an arithmetic adjacency scan, not pattern matching).

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Extended permutation and breakpoint definition (|Δ| ≠ 1, unsigned) [1][3].
- Lower bound d ≥ b/2 from "a reversal removes at most two breakpoints" [2].
- Identity ⇒ 0 [1, §2].

**Intentionally simplified:**

- Uses the **unsigned** breakpoint model and returns a **lower bound** ⌈b/2⌉; **consequence:** the value may be strictly less than the true reversal distance, and gene strand orientation (sign) is ignored.

**Not implemented:**

- Exact signed Hannenhalli–Pevzner distance (cycles c and hurdles h, d = b − c + h + correction); **users should rely on:** a dedicated signed-rearrangement tool (e.g. GRIMM) for exact distances — no in-repo alternative currently.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Integer rounding ⌈b/2⌉ via (b+1)/2 | Assumption | tightest integer guaranteed by d ≥ b/2 | accepted | not an invented value; matches the theorem bound |
| 2 | Unequal-length inputs throw | Assumption | error mode not specified by sources | accepted | distance undefined across different marker sets |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| perm1 = perm2 (identity) | 0 | identity ⇒ 0 breakpoints [1] |
| empty / single element | 0 | no internal adjacency |
| fully reversed [4,3,2,1] → [1,2,3,4] | 1 | b = 2 ⇒ ⌈2/2⌉ = 1 |
| different lengths | `ArgumentException` | distance undefined |

### 6.2 Limitations

The result is a **lower bound**, not the exact reversal distance; it can be strictly smaller than the true number of reversals (the bound "is not very tight" [2]). The unsigned model ignores gene orientation; for signed/oriented genomes the exact distance requires the breakpoint-graph cycle/hurdle analysis [1][2].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
int d = ComparativeGenomics.CalculateReversalDistance(
    new[] { 2, 3, 1, 6, 5, 4 },
    new[] { 1, 2, 3, 4, 5, 6 });
// d == 2  (extended unsigned permutation has b = 4 breakpoints ⇒ ⌈4/2⌉ = 2)
```

**Numerical walk-through:** For perm1 = [2,3,1,6,5,4], target = [1..6], the relative permutation is [1,2,0,5,4,3]. Extended: (−1, 1, 2, 0, 5, 4, 3, 6). Breakpoints: (−1,1) Δ=2 ✔, (1,2) Δ=1, (2,0) Δ=2 ✔, (0,5) Δ=5 ✔, (5,4) Δ=1, (4,3) Δ=1, (3,6) Δ=3 ✔ ⇒ b = 4 ⇒ ⌈4/2⌉ = 2. (The signed worked example in [2] has b = 6; the unsigned specialization here yields b = 4.)

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_CalculateReversalDistance_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Analysis/ComparativeGenomics_CalculateReversalDistance_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [COMPGEN-REVERSAL-001-Evidence.md](../../../docs/Evidence/COMPGEN-REVERSAL-001-Evidence.md)
- Related algorithms: [Genome_Rearrangement_Detection.md](./Genome_Rearrangement_Detection.md)

## 8. References

1. Bafna V, Pevzner PA. 1998. Sorting by Transpositions. *SIAM Journal on Discrete Mathematics* 11(2):224–240. https://www.ic.unicamp.br/~meidanis/courses/mo640/2008s2/textos/Bafna-Pevzner-1998.pdf
2. Hunter College Computational Biology. Lecture 16: Genome rearrangements, sorting by reversals. https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf
3. Hübotter J. 2020. On Sorting by Reversals. https://jonhue.github.io/min-sbr/paper.pdf
4. Bergeron A, Mixtacki J, Stoye J. 2009. The Inversion Distance Problem. https://gi.cebitec.uni-bielefeld.de/_media/teaching/2018winter/cg/inversionbergeron.pdf
