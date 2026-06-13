# Genome Rearrangement Detection (Breakpoints)

| Field | Value |
|-------|-------|
| Algorithm Group | Comparative Genomics |
| Test Unit ID | COMPGEN-REARR-001 |
| Related Projects | Seqeron.Genomics.Analysis |
| Implementation Status | Simplified |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Genome rearrangement detection compares the gene order of two genomes and reports the
positions where their order is disrupted. Given orthologous markers in genome-1 order, the
algorithm builds the corresponding signed permutation, extends it with the standard sentinels
`0` and `n+1`, and reports every consecutive pair that is not an identity adjacency as a
*breakpoint* [1][3]. Each breakpoint marks a rearrangement boundary; the number of breakpoints
is the breakpoint distance `d_BP = n − (common adjacencies)`, a lower bound on the reversal
distance (`d ≥ b/2`) [2][3]. The procedure is exact and deterministic for the breakpoint model;
it classifies each boundary as an inversion or a transposition from its local signed signature.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Genomes evolve by inversions and transpositions, as well as by deletions, insertions, and
duplications of fragments [1]. When the order of conserved markers is more stable than their
sequence, gene-order comparison detects large-scale rearrangements that sequence alignment
misses [1]. A genome's marker order is modelled as a *signed permutation*: each marker is an
integer carrying a `+`/`−` sign for its strand [3].

### 2.2 Core Model

A signed permutation `α = (α(1), …, α(n))` over `{1,…,n}` is extended to
`(α(0), α(1), …, α(n), α(n+1)) = (0, α(1), …, α(n), n+1)` [1][3]. The target genome is
relabelled to the identity `β = (+1, +2, …, +n)`.

A consecutive pair `(x, y)` of the extended permutation is a **breakpoint** of `α` with respect
to `β` iff "neither `(x, y)` nor `(−y, −x)` appear in (extended) `β`" [3]. Because `β`'s
adjacencies are exactly the pairs `(i, i+1)`, this reduces to the single signed test
**`y ≠ x + 1`**: a reversal negates the signs of the block it reverses, so a reversed segment
`+k, +(k+1)` becomes `−(k+1), −k`, for which `y = x + 1` still holds and no internal breakpoint
is created [3]. The breakpoint count is `b(α)`, with `b(β) = 0` [3], and equals the breakpoint
distance `d_BP(π₁, π₂) = n − sim(π₁, π₂)` where `sim` is the number of common adjacencies [2].

Operation signatures used for classification:
- **Inversion (reversal):** `α[i,j] = (…, −α(j), …, −α(i), …)` reverses a block and **negates its
  signs** [3]; a breakpoint flanked by a sign change is an inversion boundary.
- **Transposition:** `ρ(i,j,k)` "inserts an interval `[i, j−1]` … moving genes … to a new
  location" **preserving their orientation** [1]; an orientation-preserving (sign-stable)
  discontinuity is a transposition boundary.

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Identical signed gene order ⇒ 0 breakpoints. | `b(β) = 0` [3]. |
| INV-02 | A pair `(x,y)` is a non-breakpoint iff `y = x + 1` (subsumes the `(x,y)` and `(−y,−x)` clauses). | Reversals negate signs [3]. |
| INV-03 | The permutation is extended with `0` and `n+1`; boundary pairs participate. | Extended-permutation definition [1][3]. |
| INV-04 | Breakpoint count ∈ `[0, n+1]` for `n` markers (n+1 internal pairs of the extended permutation). | n+1 consecutive pairs exist [3]. |
| INV-05 | A boundary with a sign reversal ⇒ Inversion; an orientation-preserving discontinuity ⇒ Transposition. | Reversal negates signs [3]; transposition preserves orientation [1]. |

### 2.5 Comparison with Related Methods

| Aspect | Breakpoint detection (this) | Reversal distance (Hannenhalli–Pevzner) |
|--------|-----------------------------|------------------------------------------|
| Quantity | counts disrupted adjacencies | minimum number of reversals |
| Cost | `O(n)` after relabelling | polynomial (cycle/graph analysis) |
| Relation | `d_BP` is a lower bound | `d ≥ d_BP / 2` [3] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `genome1Genes` | `IReadOnlyList<Gene>` | required | Genome-1 genes in chromosomal order. | Non-null; `Strand` ∈ {`'+'`,`'-'`}. |
| `genome2Genes` | `IReadOnlyList<Gene>` | required | Genome-2 genes in chromosomal order (the reference/identity). | Non-null. |
| `orthologMap` | `IReadOnlyDictionary<string,string>` | required | genome-1 gene id → genome-2 gene id (anchors). | Non-null. |
| `rearrangement` | `RearrangementEvent` | required | A breakpoint event for `ClassifyRearrangement`. | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (return) | `IEnumerable<RearrangementEvent>` | One event per breakpoint, in genome-1 order. `Position` is the genome-1 `Start` of the marker just after the boundary; `Length` is `|x − y|` of the signed pair; `TargetPosition` is `"x->y"`. |
| (return) | `RearrangementType` | `ClassifyRearrangement` returns Inversion or Transposition (per INV-05). |

### 3.3 Preconditions and Validation

`DetectRearrangements` throws `ArgumentNullException` for any null argument (validated eagerly
before iteration). Fewer than 2 mappable orthologs yields no events (no internal adjacency).
Orthologs whose target gene is absent in genome 2 are skipped. Gene coordinates are taken as
0-based `Start`/`End`; strand `'+'`/`'-'` maps to the permutation sign. Markers common to both
genomes are relabelled to their genome-2 order rank (the relative permutation).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate arguments (null checks).
2. Index genome-2 genes by 1-based rank.
3. Read genome-1 orthologs in order; for each, record `sign × rank` (sign from relative strand).
4. If fewer than 2 markers, return no events.
5. Relabel the markers to their genome-2 order positions `1..n` (relative signed permutation).
6. Walk the extended sequence `[0, relabelled…, n+1]`; for each consecutive pair `(x, y)` with
   `y ≠ x + 1`, emit a breakpoint event classified by its signed signature.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- **Breakpoint rule:** pair `(x, y)` is a breakpoint iff `y ≠ x + 1` (INV-02) [3].
- **Sentinels:** `LeftSentinel = 0`, right sentinel `= n + 1` [1].
- **Classification:** sign reversal across the boundary (or a negative value involved) ⇒
  Inversion; orientation-preserving positive discontinuity ⇒ Transposition (INV-05) [1][3].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| `DetectRearrangements` | O(n log n) | O(n) | n = common markers; the `log n` term is the relabelling sort, the walk is O(n). |
| `ClassifyRearrangement` | O(1) | O(1) | parses the stored signed pair. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [ComparativeGenomics.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ComparativeGenomics.cs)

- `ComparativeGenomics.DetectRearrangements(genome1Genes, genome2Genes, orthologMap)`: returns one breakpoint event per disrupted adjacency.
- `ComparativeGenomics.ClassifyRearrangement(rearrangement)`: maps a breakpoint event to its `RearrangementType`.

### 5.2 Current Behavior

The signed pair spanning each breakpoint is stored in `RearrangementEvent.TargetPosition` as
`"x->y"` so `ClassifyRearrangement` can re-derive the type without recomputation. The search is a
single linear scan over the relative permutation; no suffix tree is used because this is not a
substring/occurrence search — it is an adjacency comparison over an already-ordered marker list
(see §9 search-reuse note). Markers are relabelled to their genome-2 order rank so the breakpoint
test is the exact signed-permutation criterion regardless of absolute coordinates.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Extended signed permutation with sentinels `0` and `n+1` [1][3].
- Breakpoint criterion `(x,y)` is a breakpoint iff neither `(x,y)` nor `(−y,−x)` is an identity adjacency, realised as `y ≠ x + 1` [3].
- Inversion = sign-negating reversal; Transposition = orientation-preserving relocation [1][3].

**Intentionally simplified:**

- Classification uses the **local** signed signature of each boundary; **consequence:** a single multi-anchor scenario is reported as several per-boundary events rather than one global operation, so the user sees breakpoint boundaries, not a minimal rearrangement scenario.

**Not implemented:**

- Translocation, Deletion, Insertion, Duplication classification; **users should rely on:** chromosome-aware or gene-set-difference methods (these require chromosome identifiers or marker-count differences not expressible in a single in-order signed permutation — Evidence Assumption 3). Reversal *distance* (minimum number of reversals) is the separate unit COMPGEN-REVERSAL-001 (`CalculateReversalDistance`).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Anchors supplied as ordered ortholog map | Assumption | anchor generation delegated to ortholog/synteny units | accepted | Evidence Assumption 1 |
| 2 | `Gene.Strand` encodes the permutation sign | Assumption | strand drives inversion classification | accepted | Evidence Assumption 2 |
| 3 | Only Inversion/Transposition classified | Deviation | other enum types out of scope for these methods | accepted | Evidence Assumption 3; §5.3 Not implemented |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Identical gene order | 0 events | `b(β) = 0` [3]. |
| Fewer than 2 mappable orthologs | 0 events | no internal adjacency [3]. |
| Ortholog target absent in genome 2 | anchor skipped, no crash | robustness; mirrors synteny unit. |
| Null argument | `ArgumentNullException` | contract §3.3. |
| Empty genomes | 0 events, no exception | contract §3.3. |

### 6.2 Limitations

Operates on a single (unichromosomal) in-order signed permutation of common markers. It does not
distinguish translocations (needs chromosome ids), or indels/duplications (need marker-count
differences); it reports breakpoint boundaries, not a minimum-length rearrangement scenario; and
it assumes anchors are already 1:1 orthologs.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
// genome1 g0..g5 mapped so the signed permutation is (-2,-3,+1,+6,-5,-4) vs identity (Hunter Lecture 16).
var events = ComparativeGenomics.DetectRearrangements(genome1, genome2, orthologMap).ToList();
// events.Count == 6   (b(alpha) = 6)
var type = ComparativeGenomics.ClassifyRearrangement(events[0]); // Inversion or Transposition
```

**Numerical / biological walk-through:**

For `α = (−2, −3, +1, +6, −5, −4)`, extended `(0, −2, −3, +1, +6, −5, −4, +7)`, the consecutive
pairs failing `y = x + 1` are `(0,−2), (−2,−3), (−3,+1), (+1,+6), (+6,−5), (−4,+7)` → **6
breakpoints**; `(−5,−4)` passes (`−4 = −5 + 1`) and is not a breakpoint — exactly the Hunter
Lecture 16 example [3].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [ComparativeGenomics_DetectRearrangements_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/ComparativeGenomics_DetectRearrangements_Tests.cs) — covers `INV-01`..`INV-05`
- Evidence: [COMPGEN-REARR-001-Evidence.md](../../../docs/Evidence/COMPGEN-REARR-001-Evidence.md)
- Related algorithms: [Synteny_Block_Detection](./Synteny_Block_Detection.md), [Ortholog_Identification](./Ortholog_Identification.md)

## 8. References

1. Bafna V, Pevzner PA. 1998. Sorting by Transpositions. *SIAM Journal on Discrete Mathematics* 11(2):224–240. https://doi.org/10.1137/S089548019528280X
2. Tannier E, Zheng C, Sankoff D. 2009. Multichromosomal median and halving problems under different genomic distances. *BMC Bioinformatics* 10:120 — definitions via PMC "On the Complexity of Rearrangement Problems under the Breakpoint Distance": https://pmc.ncbi.nlm.nih.gov/articles/PMC3887456/
3. Hunter College CSCI Computational Biology. Lecture 16: Genome rearrangements, sorting by reversals (exposition of Hannenhalli–Pevzner / Bafna–Pevzner). https://www.cs.hunter.cuny.edu/~saad/courses/compbio/lectures/lecture16.pdf
