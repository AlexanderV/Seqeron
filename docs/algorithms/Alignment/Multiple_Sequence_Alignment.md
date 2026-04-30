# Multiple Sequence Alignment (MSA)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-MULTI-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Simplified |
| Last Reviewed | 2026-04-30 |

## 1. Overview

Multiple sequence alignment aligns three or more sequences in one shared coordinate system by inserting gap characters until every aligned sequence has the same length. Exact optimal MSA is computationally expensive, so practical implementations usually rely on heuristics. In this repository, `SequenceAligner.MultipleAlign(...)` implements a simplified center-star workflow: it selects one center sequence, aligns every other sequence to that center through an anchor-based pairwise aligner, merges the resulting gap patterns, builds a consensus, and reports a sum-of-pairs score on the final alignment. The implementation is restricted to `DnaSequence` inputs, so the documented contract is for validated DNA strings rather than generic protein or RNA alphabets.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Given sequences $S_1, S_2, \ldots, S_k$, an MSA produces aligned sequences $S'_1, S'_2, \ldots, S'_k$ of common length $L$ such that removing gap characters from each $S'_i$ recovers the corresponding original sequence. The usual MSA representation requires a common aligned length and excludes columns that consist entirely of gaps. These are the structural properties cited in the multiple-sequence-alignment reference and enforced by the repository tests.

### 2.2 Core Model

The exact MSA optimization problem is NP-complete; the current document's cited summary states that exact dynamic programming for $n$ sequences of length $L$ requires `O(L^n)` time and space. Practical tools therefore use heuristics such as progressive or center-star alignment. In a center-star approach, one sequence is chosen as the center, each other sequence is aligned to it, and the pairwise alignments are reconciled into one multi-sequence result.

For scoring a completed alignment, the document and the repository both use column-based sum-of-pairs (SP) scoring:

$$
SP = \sum_{c=1}^{L} \sum_{i<j} score\left(S'_i[c], S'_j[c]\right)
$$

where the score function is applied to every unordered pair of characters in each alignment column. (Multiple sequence alignment)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | All aligned sequences have equal length. | That is the defining representation of an MSA and is enforced by the repository after gap reconciliation and padding. |
| INV-02 | Removing gaps from each aligned sequence recovers the original input sequence. | The algorithm only inserts gap characters; it does not rewrite the input sequence content. |
| INV-03 | No output column consists entirely of gaps. | This is a structural MSA validity condition and is asserted by the repository tests. |
| INV-04 | `Consensus.Length` equals the aligned sequence length. | The repository builds consensus one column at a time across the final aligned strings. |
| INV-05 | `TotalScore` is the SP score over the final aligned sequences. | `ComputeSumOfPairsScore` recomputes the public score from the finished alignment columns. |

### 2.5 Comparison with Related Methods

| Aspect | Repository center-star variant | Guide-tree progressive alignment |
|--------|--------------------------------|----------------------------------|
| Global structure | Pick one center sequence, then align all others to it | Build a guide tree from pairwise distances, then align along the tree |
| Sensitivity | Depends strongly on center selection | Depends on the guide tree |
| Pairwise pre-processing | Center selection by 4-mer cosine similarity in the repository | Classical progressive methods compute pairwise distances before tree construction |
| Refinement | No explicit refinement stage in the repository method | Progressive toolchains often include additional tree-driven or iterative refinement |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `sequences` | `IEnumerable<DnaSequence>` | required | Sequences to align. | The collection itself must be non-null. Individual `DnaSequence` instances are validated DNA sequences; the tests also cover empty sequences inside the collection. |
| `scoring` | `ScoringMatrix` | `SequenceAligner.SimpleDna` | Match, mismatch, and gap values used during pairwise gap alignment and final SP scoring. | Optional. |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| `AlignedSequences` | `string[]` | Final aligned sequences, all padded to the same length. |
| `Consensus` | `string` | Majority-voted consensus built column-by-column from the aligned sequences. |
| `TotalScore` | `int` | Public SP score computed on the final aligned strings. |

### 3.3 Preconditions and Validation

`SequenceAligner.MultipleAlign(...)` throws `ArgumentNullException` when the sequence collection itself is null. An empty collection returns `MultipleAlignmentResult.Empty`. A single sequence returns that sequence unchanged as both the sole aligned sequence and the consensus, with `TotalScore = 0`. `DnaSequence` normalizes inputs to uppercase, validates `A`, `C`, `G`, and `T`, and allows empty strings; the MSA tests explicitly cover empty sequences inside the collection.

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate the input collection and handle trivial empty or single-sequence cases.
2. Select a center sequence with the highest total 4-mer cosine similarity to the other sequences.
3. Build one suffix tree on the center sequence.
4. Align each remaining sequence to the center using the anchor-based pairwise aligner.
5. Merge the independent center-vs-other gap patterns into a single coordinate system and pad all sequences to equal length.
6. Build a consensus over the final columns and compute the SP score on the final aligned strings.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

`SelectCenterSequence` hashes DNA 4-mers into a fixed 256-entry profile and uses cosine similarity to score every sequence pair. `MultipleAlign` then builds a suffix tree on the chosen center sequence and calls `AnchorBasedAligner.AlignWithAnchors(...)` for each remaining sequence. During consensus construction, gaps participate in the vote and ties between a nucleotide and `-` are resolved in favor of the nucleotide. `ComputeSumOfPairsScore` evaluates every unordered sequence pair in every column, using `GapExtend` for gap-versus-nucleotide columns and `0` for gap-versus-gap columns.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Exact MSA dynamic programming (theoretical baseline) | `O(L^n)` | `O(L^n)` | Complexity quoted in the current document from Wang and Jiang (1994). |
| Repository center-sequence selection | `O(k^2 * L)` | `O(k * 4^4)` | `SelectCenterSequence` builds one 256-entry 4-mer profile per input sequence. |
| Repository anchor-based pairwise phase | `O(L)` to build the center suffix tree, plus `O(k * L + \sum_i \delta_i^2)` across anchor finding and gap alignment | Not explicitly characterized in current source comments | `SequenceAligner.MultipleAlign` and `AnchorBasedAligner` document this split. |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`: public MSA entry point.
- `SequenceAligner.SelectCenterSequence(...)`: chooses the center sequence using 4-mer cosine similarity.
- `SequenceAligner.ReconcileAlignments(...)`: merges the center-vs-other pairwise results into one final alignment.
- `SequenceAligner.BuildConsensus(...)`: builds the public consensus string from the final columns.
- `SequenceAligner.ComputeSumOfPairsScore(...)`: computes the public `TotalScore` from the final aligned strings.
- `SequenceAligner.MultipleAlignClassic(...)`: internal classic star-alignment helper retained in the same file but not used by the public `MultipleAlign` path.

**Anchor-based pairwise helper:** [AnchorBasedAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AnchorBasedAligner.cs)

- `AnchorBasedAligner.AlignWithAnchors(...)`: finds exact-match anchors on the center sequence and aligns intervening gaps.
- `AnchorBasedAligner.ChainAnchors(...)`: keeps a consistent increasing chain of anchors.
- `AnchorBasedAligner.AlignGap(...)`: either aligns gap segments with `GlobalAlign` or falls back to simpler cases.

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

`MultipleAlign` handles empty and single-sequence inputs before any heuristic work. For two or more sequences, it selects a center by 4-mer cosine similarity rather than by pairwise dynamic-programming distances. Pairwise center alignments are reconciled in parallel when `k >= 4` and sequentially otherwise. If `AnchorBasedAligner` finds no usable anchors for a pair, it falls back to `SequenceAligner.GlobalAlign(...)` for that pairwise step. The public `TotalScore` is always recomputed from the final aligned strings by `ComputeSumOfPairsScore`; it is not the sum of the intermediate pairwise scores accumulated during reconciliation. The final consensus includes gaps in the vote and prefers a nucleotide over a gap on ties.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Common-length aligned outputs with gap removal restoring the original sequences.
- Column-based SP scoring on the final multi-sequence alignment.
- A center-star workflow in which one sequence is chosen as the hub and the remaining sequences are aligned to it.

**Intentionally simplified:**

- The repository does not build a guide tree; it chooses the center by 4-mer cosine similarity; **consequence:** the output depends on that center-selection heuristic rather than on a tree derived from pairwise distances.
- The public path uses anchor-based pairwise alignment against the center instead of exact all-sequence dynamic programming; **consequence:** the result is heuristic rather than globally optimal for the full MSA objective.

**Not implemented:**

- Exact optimal MSA or guide-tree progressive refinement; **users should rely on:** no current alternative in `SequenceAligner.MultipleAlign`.

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null collection | Throws `ArgumentNullException`. | Public API guard. |
| Empty collection | Returns `MultipleAlignmentResult.Empty`. | Explicit early return in `MultipleAlign`. |
| Single sequence | Returns the original sequence as the only aligned sequence and as the consensus, with `TotalScore = 0`. | There are no pairwise columns to score. |
| Different-length inputs | Shorter sequences are padded with gaps in the final alignment. | Gap reconciliation and final padding normalize all outputs to one length. |
| Empty sequence inside the collection | Handled without rejecting the whole alignment. | Covered by `SequenceAligner_MultipleAlign_Tests`. |
| Gap-versus-nucleotide tie in consensus voting | Consensus prefers the nucleotide. | Implemented explicitly in `BuildConsensus`. |

### 6.2 Limitations

The result depends on the chosen center sequence, so center-selection errors propagate into the final alignment. The repository does not build a phylogenetic guide tree and does not perform a later refinement pass. The algorithm is documented and tested for `DnaSequence` inputs, not generic amino-acid or RNA alphabets. Because the public alignment is heuristic, the repository does not claim global optimality for the full multi-sequence objective.

## 8. References

1. [Multiple sequence alignment](https://en.wikipedia.org/wiki/Multiple_sequence_alignment)
2. [Clustal](https://en.wikipedia.org/wiki/Clustal)
3. [Consensus sequence](https://en.wikipedia.org/wiki/Consensus_sequence)
4. Feng, D. F.; Doolittle, R. F. (1987). "Progressive sequence alignment as a prerequisite to correct phylogenetic trees." Journal of Molecular Evolution 25(4): 351-360.
5. Wang, L.; Jiang, T. (1994). "On the complexity of multiple sequence alignment." Journal of Computational Biology 1(4): 337-348.
6. Thompson, J. D.; Higgins, D. G.; Gibson, T. J. (1994). "CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice." Nucleic Acids Research 22(22): 4673-4680.
