# Multiple Sequence Alignment (MSA)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-MULTI-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Production (linear-gap model; declared deviations in §5.3) |
| Last Reviewed | 2026-11-04 (campaign 2026-09, batch B09) |

## 1. Overview

Multiple sequence alignment aligns three or more sequences in one shared coordinate system by inserting gap characters until every aligned sequence has the same length. Exact optimal MSA is computationally expensive, so practical implementations usually rely on heuristics. This repository provides five aligners over `DnaSequence` inputs, all in the **linear gap model** of `GlobalAlign` (every gap position scores `GapExtend`; `GapOpen` is not used):

1. `SequenceAligner.MultipleAlign(...)` — **center-star** alignment (Gusfield 1993) with *exact* pairwise alignments: all pair scores come from the canonical Needleman-Wunsch (`GlobalAlign` core), the center is the sequence with the highest total pair score, every other sequence is aligned to it with the same DP and the pairwise alignments are merged column-wise ("once a gap, always a gap"). For two sequences the result is an optimal pairwise global alignment.
2. `SequenceAligner.MultipleAlignAnchored(...)` — the approximate, faster variant of (1): 4-mer-cosine center (counts from the Core `CountKmersSpan`), suffix-tree exact-match anchors (`AnchorBasedAligner`: MUMmer-style seeds, weighted colinear chain with overlap trimming, Needleman-Wunsch on the gaps, exact Hirschberg for very large gaps). Pairwise alignments are *not* guaranteed optimal.
3. `SequenceAligner.MultipleAlignProgressive(...)` — guide-tree progressive alignment: ClustalW pairwise distances (1 − identities / shorter length) → UPGMA guide tree (Higgins & Sharp 1988 CLUSTAL; Feng & Doolittle 1987 progressive scheme) → profile-profile Needleman-Wunsch, single pass.
4. `SequenceAligner.MultipleAlignIterative(...)` — adds **iterative refinement** of the progressive seed by MUSCLE-style tree-dependent restricted partitioning (Edgar 2004), accepting a re-alignment only when the SP score strictly increases.
5. `SequenceAligner.MultipleAlignConsistency(...)` — **consistency-based** alignment (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of percent-identity-weighted global and local pairwise alignments, the triplet extension evaluated for every residue pair, and a progressive alignment that scores columns with the extended-library weights (zero gap penalties) — a distinct objective from the sum-of-pairs (SP) score the other four optimise.

The implementation is restricted to `DnaSequence` inputs, so the documented contract is for validated DNA strings rather than generic protein or RNA alphabets. MCP tools: `multiple_align`, `multiple_align_anchored`, `multiple_align_progressive`, `multiple_align_iterative`, `multiple_align_consistency`.

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
| INV-06 | `MultipleAlignIterative` SP score is ≥ the `MultipleAlignProgressive` seed SP score (monotonic non-decreasing). | Each refinement re-alignment is accepted only on a strict SP improvement (Edgar 2004, step 3.4); the seed is the lower bound. |
| INV-07 | `MultipleAlignIterative` is deterministic and converges. | Edges are visited in a fixed order (decreasing distance from the root, deterministic tie-break); no RNG; a pass with no accepted change stops the loop, bounded by a positive iteration cap. |
| INV-08 | The T-Coffee extended-library weight of a residue pair equals its primary (direct) weight plus the sum over intermediate sequences of `min(W1, W2)`, and is therefore ≥ the primary weight (extension never lowers a weight). | The triplet consistency transformation only adds non-negative `min`-of-legs contributions to the direct weight (Notredame et al. 2000, p.209; worked example 88 → 165). |
| INV-09 | A consistency-supported residue pair has a strictly greater extended weight than an unsupported (inconsistent) pair. | An unsupported pair is absent from the library (weight 0); a supported pair carries the positive direct weight plus triplet support (Notredame et al. 2000, p.209). |
| INV-10 | `MultipleAlignConsistency` is deterministic and produces a structurally valid MSA (equal-length rows, each degaps to its input, no all-gap column, count preserved). | No RNG; primary alignments, UPGMA tree, extension sum and DP tie-breaks are all deterministic; gaps are inserted only as whole all-gap columns. |

| INV-11 | For two sequences `MultipleAlign`, `MultipleAlignProgressive` and `MultipleAlignIterative` return an optimal pairwise global alignment (SP = Needleman-Wunsch optimum = Biopython `PairwiseAligner` with `open_gap_score = extend_gap_score = GapExtend`). | The star and the two-leaf progressive/iterative merges are exactly one canonical NW / profile NW on two single rows. `MultipleAlignAnchored` returns a score ≤ the optimum; `MultipleAlignConsistency` maximises library weights (zero gap cost), not the NW score. |
| INV-12 | The `MultipleAlign` center is argmax over sequences of the total optimal pair score to the others (ties: lowest index). | Gusfield 1993 center-star definition (the string minimising the total pairwise distance). |

### 2.5 Comparison with Related Methods

| Aspect | Repository center-star variant | Guide-tree progressive alignment |
|--------|--------------------------------|----------------------------------|
| Global structure | Pick one center sequence, then align all others to it | Build a guide tree from pairwise distances, then align along the tree |
| Sensitivity | Depends strongly on center selection | Depends on the guide tree |
| Pairwise pre-processing | All C(k,2) optimal pair scores pick the center (Gusfield); `MultipleAlignAnchored` uses a 4-mer cosine heuristic instead | Classical progressive methods compute pairwise distances before tree construction |
| Refinement | The star method has no refinement; the progressive method is single-pass; `MultipleAlignIterative` adds MUSCLE-style tree-dependent iterative refinement (Edgar 2004) | Progressive toolchains often include additional tree-driven or iterative refinement |

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

`MultipleAlign` (center-star, Gusfield 1993):

1. Validate the input collection and handle trivial empty or single-sequence cases.
2. Compute the optimal Needleman-Wunsch score of every pair (canonical `GlobalAlignCore`; parallel over pairs).
3. Choose the center: the sequence with the highest total pair score (ties: lowest input index).
4. Align every other sequence to the center with the canonical Needleman-Wunsch (center = sequence 1).
5. Merge the independent center-vs-other gap patterns: for every center position the number of inserted columns is the maximum over rows; each row's inserted residues are left-justified in the slot (`MergeCenterAlignments`).
6. Build the majority consensus and the column sum-of-pairs score of the final strings.

`MultipleAlignAnchored` replaces steps 2-4 by the 4-mer cosine center and `AnchorBasedAligner.AlignWithAnchors` (suffix tree of the center, `FindExactMatchAnchors`, weighted chain, NW between anchors), then continues with steps 5-6.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

`SelectCenterSequence` (anchored variant only) counts the 4-mers of every sequence with `SequenceExtensions.CountKmersSpan` (DNA input is upper-case ACGT, so every window is a valid 4-mer; the earlier private 256-slot hash gives identical counts) and picks the sequence with the highest total cosine similarity of the count vectors. During consensus construction gaps participate in the vote; ties between a nucleotide and `-` go to the nucleotide. `ComputeSumOfPairsScore` evaluates every unordered sequence pair in every column: match/mismatch from the matrix, `GapExtend` for gap-versus-nucleotide, `0` for gap-versus-gap.

**Anchor chaining (`AnchorBasedAligner.ChainAnchors`).** The anchors are maximal exact matches, one per matching diagonal run found by suffix-link streaming (not the full MEM set), so anchors from different diagonals may overlap. The chain is the maximum-total-length colinear chain (Abouelhoda & Ohlebusch 2005: maximum-weight chain of colinear non-overlapping fragments, weight = length); an anchor that overlaps its predecessor is trimmed at the front by the smallest amount that makes it non-overlapping (a suffix of an exact match is an exact match). MUMmer's longest increasing subsequence (Delcher et al. 1999) is the unit-weight, overlap-free special case. The regions between chained anchors are aligned with the canonical NW (linear gaps, same scoring as `GlobalAlign`); segments above 25 M cells use the exact linear-space Hirschberg recursion (`GlobalAlignLinearSpace`). The reported `Score` is the linear score of the returned alignment and never exceeds the unconstrained optimum.

**Progressive distance.** `1 − identical aligned residue pairs / min(|a|, |b|)` of the optimal pairwise alignment — the ClustalW distance (`pairalign.c`, clustalw 1.83: `tmat = (100 − 100·count/MIN(len1,len2)) / 100`). Equal averaged distances tie (within 1e-9) and the lowest index pair merges first.

**T-Coffee library.** Primary weight of every aligned residue pair of a pairwise alignment = ⌊100 · identical columns / aligned (non-gap) columns⌋ — `get_seq_sim`, mode sim1, `lib/util_lib/aln_convertion_util.c` of the T-Coffee sources; global (NW) and local (SW) entries of the same residue pair are signal-added. Extended weight of any pair (Si.p, Sj.q) = direct weight + Σ over residues Sk.r of third sequences with both legs in the primary library of min(W1, W2) (`residue_pair_extended_list`, `lib/dp_lib/evaluate.c`); a pair without a direct entry but with triplet support therefore has a positive extended weight.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Exact MSA dynamic programming (theoretical baseline) | `O(L^n)` | `O(L^n)` | Complexity quoted in the current document from Wang and Jiang (1994). |
| `MultipleAlign` (exact center-star) | `O(k^2 * L^2)` pair scores + `O(k * L^2)` center alignments | `O(L^2)` per concurrent pair | Canonical NW per pair, parallel over pairs. |
| `MultipleAlignAnchored` center selection | `O(k^2 * 4^4 + k * L)` | `O(k * min(L, 256))` | `SelectCenterSequence` counts 4-mers per sequence (sparse dictionaries). |
| `MultipleAlignAnchored` pairwise phase | `O(L)` suffix tree, `O(k * L)` anchor streaming, `O(a^2)` chaining over `a` anchors, `O(\sum_i \delta_i^2)` gap DP | `O(L)` + `O(\delta^2)` (Hirschberg above 25 M cells) | `AnchorBasedAligner` documents this split. |
| Iterative refinement (`MultipleAlignIterative`) | Progressive seed `O(k^2 L^2)` + refinement `O(I * E * (k * L^2))` where `E = O(k)` internal edges, `I = maxIterations` passes | `O(k * L)` | Each refinement re-alignment is one profile-profile NW (`O(k L^2)` over the column DP); passes stop early on convergence. Measured baseline: the full ALIGN-MULTI-001 iterative fixture (15 tests incl. a 500-trial random property test, k≤5, L≤8) completes in well under 1 s on the dev machine (2026-06-23). |
| Consistency-based (`MultipleAlignConsistency`) | Primary library + extension `O(k^3 L^2)` worst case; progressive DP `O(k^2 L^2)` | `O(k^2 L^2)` for the library | Library extension is the dominant term (Notredame et al. 2000 quote `O(N^3 L^2)` worst case, closer to `O(N^3 L)` in practice). Measured baseline: the full ALIGN-MULTI-001 consistency fixture (12 tests incl. a 100-trial random property test, k≤5, L≤8) completes in well under 1 s on the dev machine (2026-06-23). |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`: exact center-star (Gusfield 1993); public MSA entry point.
- `SequenceAligner.MultipleAlignAnchored(IEnumerable<DnaSequence>, ScoringMatrix?)`: approximate suffix-tree-anchored center-star.
- `SequenceAligner.SelectCenterByPairwiseScore(...)` / `SelectCenterSequence(...)`: Gusfield pair-score center / 4-mer cosine center (anchored variant; k-mers from the Core `CountKmersSpan`).
- `SequenceAligner.MergeCenterAlignments(...)`: merges the center-vs-other pairwise results into one alignment.
- `SequenceAligner.BuildConsensus(...)`, `ComputeSumOfPairsScore(...)`: consensus and public `TotalScore`.
- `SequenceAligner.MultipleAlignClassic(...)`: internal benchmark baseline (first sequence as center, exact pairwise alignments).
- `SequenceAligner.MultipleAlignProgressive(IEnumerable<DnaSequence>, ScoringMatrix?)`: ClustalW-distance / UPGMA / profile-profile NW progressive MSA.
- `SequenceAligner.MultipleAlignIterative(IEnumerable<DnaSequence>, ScoringMatrix?, int)`: MUSCLE-style tree-dependent iterative refinement of the progressive seed (Edgar 2004, Stage 3). `maxIterations` (default 16) caps the refinement passes.
- `SequenceAligner.MultipleAlignConsistency(IEnumerable<DnaSequence>, ScoringMatrix?)`: T-Coffee consistency MSA. `BuildPrimaryLibrary` (global + local pairwise alignments weighted by `PercentIdentity`, signal-added) feeds the `ExtendedLibrary` (triplet extension evaluated on demand for any residue pair), and `AlignConsistencyProfiles` / `ColumnLibraryScore` run the zero-gap-penalty progressive DP. `scoring` seeds the pairwise library and the reported SP score; the DP uses library weights. `GetLibraryWeights(...)` is an internal accessor for tests.
- `SequenceAligner.RefineByTreePartitioning(...)`, `EnumerateEdgePartitions(...)`, `SplitProfile(...)`: iterative-refinement helpers.
- `SequenceAligner.GlobalAlignLinearSpace(...)`: internal exact Hirschberg (1975) linear-space global alignment used for very large inter-anchor segments.

**Anchor-based pairwise helper:** [AnchorBasedAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/AnchorBasedAligner.cs)

- `AnchorBasedAligner.AlignWithAnchors(...)`: exact-match anchors from the center's suffix tree, chain, gap fill; an empty side yields an all-gap alignment.
- `AnchorBasedAligner.ChainAnchors(...)`: maximum-total-length colinear chain with overlap trimming (§4.2).

**Supporting types:** [AlignmentTypes.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/AlignmentTypes.cs), [DnaSequence.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/DnaSequence.cs)

### 5.2 Current Behavior

All five aligners handle empty and single-sequence inputs before any heuristic work. The public `TotalScore` is always recomputed from the final aligned strings by `ComputeSumOfPairsScore`. The final consensus includes gaps in the vote and prefers a nucleotide over a gap on ties. Validation (independent Python references, B09): the rows of `MultipleAlign`, `MultipleAlignProgressive`, `MultipleAlignIterative` and `MultipleAlignConsistency` equal an independent re-implementation of each method on 300 random sets (k = 2..6, three scoring schemes) with no mismatch; every alignment degaps to its input, has no all-gap column and re-scores to `TotalScore`; for k = 2 the star, progressive and iterative scores equal the Biopython 1.88 `PairwiseAligner` optimum.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Common-length aligned outputs with gap removal restoring the original sequences.
- Column-based SP scoring on the final multi-sequence alignment (match/mismatch from the matrix, residue-gap = `GapExtend`, gap-gap = 0).
- Center-star alignment (Gusfield 1993): optimal pairwise alignments to the pair-score center, merged by once-a-gap-always-a-gap (`MultipleAlign`).
- Progressive alignment along a guide tree (Feng & Doolittle 1987; Higgins & Sharp 1988; Thompson et al. 1994): all pairwise NW alignments → ClustalW distance → UPGMA → profile-profile NW, single pass (`MultipleAlignProgressive`).
- MUSCLE Stage 3 "tree-dependent restricted partitioning" iterative refinement (Edgar 2004): visit guide-tree edges in order of decreasing distance from the root, split into two sub-profiles, re-align with profile-profile NW, keep only on a strictly higher SP score, repeat until convergence or the iteration cap (`MultipleAlignIterative`).
- T-Coffee (Notredame, Higgins & Heringa 2000, `MultipleAlignConsistency`): primary library of global (NW) and local (SW) pairwise alignments weighted by percent identity (aligned columns only, integer division — T-Coffee `get_seq_sim` sim1) and combined by signal addition (p.207); library extension by the triplet consistency transformation for every residue pair, extended weight = direct weight + Σ min(W1, W2) (p.209, worked example 88 → 165; `residue_pair_extended_list`); progressive alignment with the extended-library weights as position-specific substitution scores and zero gap penalties (p.210).

**Intentionally simplified / declared deviations (each with its consequence):**

- **Linear gap model.** `GapOpen` is not used by any MSA method (as in `GlobalAlign`/`LocalAlign`). The sum-of-pairs objective with affine *natural* gap costs (Altschul 1989) charges a gap opening per sequence pair whose projection starts a new run, which depends on the whole column history of both rows; MSA programs (MUSCLE, ClustalW, MAFFT) therefore use position-specific or quasi-natural approximations, not the exact natural cost, and no exact polynomial profile-profile DP for it is known. Consequence: MSAs of sequences that differ by long indels are scored with a per-position gap cost. For two sequences use `GlobalAlignAffine`.
- The star methods build no guide tree; the progressive family uses UPGMA (Higgins & Sharp 1988 CLUSTAL) rather than the neighbor-joining tree of ClustalW/T-Coffee, reusing one tree builder. Consequence: the merge order can differ from ClustalW/T-Coffee; the library and objective of T-Coffee are unchanged.
- The ClustalW distance is computed from the repository's linear NW alignment, not ClustalW's affine Myers-Miller alignment.
- Iterative refinement partitions by guide-tree edge (MUSCLE) rather than by removing one sequence at a time (Barton-Sternberg 1987); acceptance semantics are identical.
- The T-Coffee primary library uses the single best local alignment per pair (not the top-N non-intersecting Lalign alignments), and no sequence-redundancy weights (`1/Σ sim³`).
- `MultipleAlignAnchored` is a heuristic: anchor chaining constrains the path, so pairwise alignments (and the SP score) can fall below `MultipleAlign`'s.

**Not implemented:**

- Exact optimal MSA (NP-complete; Wang & Jiang 1994).

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null collection | Throws `ArgumentNullException`. | Public API guard. |
| Empty collection | Returns `MultipleAlignmentResult.Empty`. | Explicit early return. |
| Single sequence | Returns the original sequence as the only aligned sequence and as the consensus, with `TotalScore = 0`. | There are no pairwise columns to score. |
| Different-length inputs | Shorter sequences are padded with gaps in the final alignment. | Gap reconciliation normalises all outputs to one length. |
| Empty sequence inside the collection | An all-gap row; the other rows are unaffected. | Covered by `SequenceAligner_MultipleAlign_Tests`. |
| Gap-versus-nucleotide tie in consensus voting | Consensus prefers the nucleotide. | Implemented explicitly in `BuildConsensus`. |
| Two equal averaged guide-tree distances | The lowest-index pair merges first (equality within 1e-9). | Deterministic, floating-point-robust UPGMA. |

### 6.2 Limitations

The star result depends on the chosen center and is not SP-optimal (Gusfield's 2 − 2/k bound is for metric distances). The progressive method is single-pass and inherits the "once a gap, always a gap" behaviour; `MultipleAlignIterative` removes that limitation by re-splitting and re-aligning, but it is SP-guided and restricted to guide-tree edge partitions, so it improves on but does not guarantee the global SP optimum. `MultipleAlignConsistency` optimises the T-Coffee consistency objective over a greedy progressive merge order, so it does not guarantee the global optimum of that objective. All aligners are documented and tested for `DnaSequence` inputs only, and the linear-gap deviation of §5.3 applies to all of them.

## 7. Examples and Related Material

### 7.3 Related Tests, Evidence, or Documents

- Tests (consistency / T-Coffee): [SequenceAligner_MultipleAlignConsistency_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlignConsistency_Tests.cs) — covers `INV-08`, `INV-09`, `INV-10` and the primary→extended weight relation (200 → 400).
- Tests (star / anchored / chaining): [SequenceAligner_MultipleAlignStar_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlignStar_Tests.cs), [AnchorBasedAligner_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/AnchorBasedAligner_Tests.cs)
- Tests (iterative refinement): [SequenceAligner_MultipleAlignIterative_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlignIterative_Tests.cs) — covers `INV-06`, `INV-07` and the headline gap-relocation correction.
- Tests (progressive): [SequenceAligner_MultipleAlignProgressive_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlignProgressive_Tests.cs)
- Tests (star): [SequenceAligner_MultipleAlign_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAligner_MultipleAlign_Tests.cs)
- Evidence: [ALIGN-MULTI-001-Evidence.md](../../../docs/Evidence/ALIGN-MULTI-001-Evidence.md)
- TestSpec: [ALIGN-MULTI-001.md](../../../tests/TestSpecs/ALIGN-MULTI-001.md)

## 8. References

1. [Multiple sequence alignment](https://en.wikipedia.org/wiki/Multiple_sequence_alignment)
2. [Clustal](https://en.wikipedia.org/wiki/Clustal)
3. [Consensus sequence](https://en.wikipedia.org/wiki/Consensus_sequence)
4. Feng, D. F.; Doolittle, R. F. (1987). "Progressive sequence alignment as a prerequisite to correct phylogenetic trees." Journal of Molecular Evolution 25(4): 351-360.
5. Wang, L.; Jiang, T. (1994). "On the complexity of multiple sequence alignment." Journal of Computational Biology 1(4): 337-348.
6. Thompson, J. D.; Higgins, D. G.; Gibson, T. J. (1994). "CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice." Nucleic Acids Research 22(22): 4673-4680.
7. Edgar, R. C. (2004). "MUSCLE: multiple sequence alignment with high accuracy and high throughput." Nucleic Acids Research 32(5): 1792-1797. https://academic.oup.com/nar/article/32/5/1792/2380623
8. Barton, G. J.; Sternberg, M. J. (1987). "A strategy for the rapid multiple alignment of protein sequences. Confidence levels from tertiary structure comparisons." Journal of Molecular Biology 198(2): 327-337. https://pubmed.ncbi.nlm.nih.gov/3430611/
9. Wallace, I. M.; O'Sullivan, O.; Higgins, D. G. (2005). "Evaluation of iterative alignment algorithms for multiple alignment." Bioinformatics 21(8): 1408-1414. https://academic.oup.com/bioinformatics/article/21/8/1408/249176
10. Notredame, C.; Higgins, D. G.; Heringa, J. (2000). "T-Coffee: A novel method for fast and accurate multiple sequence alignment." Journal of Molecular Biology 302(1): 205-217. https://doi.org/10.1006/jmbi.2000.4042 (full text https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf)
11. T-Coffee Technical Documentation. https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html
12. Gotoh, O. (1982). "An improved algorithm for matching biological sequences." Journal of Molecular Biology 162(3): 705-708. https://doi.org/10.1016/0022-2836(82)90398-9
13. Gusfield, D. (1993). "Efficient methods for multiple sequence alignment with guaranteed error bounds." Bulletin of Mathematical Biology 55(1): 141-154. https://doi.org/10.1007/BF02460299 (center-star method, bound 2 − 2/k)
14. Delcher, A. L.; Kasif, S.; Fleischmann, R. D.; Peterson, J.; White, O.; Salzberg, S. L. (1999). "Alignment of whole genomes." Nucleic Acids Research 27(11): 2369-2376 (MUMmer: suffix-tree matches, longest colinear subset, DP between them).
15. Abouelhoda, M. I.; Ohlebusch, E. (2005). "Chaining algorithms for multiple genome comparison." Journal of Discrete Algorithms 3(2-4): 321-341 (maximum-weight chain of colinear non-overlapping fragments; PDF opened).
16. Hirschberg, D. S. (1975). "A linear space algorithm for computing maximal common subsequences." Communications of the ACM 18(6): 341-343.
17. Higgins, D. G.; Sharp, P. M. (1988). "CLUSTAL: a package for performing multiple sequence alignment on a microcomputer." Gene 73(1): 237-244 (UPGMA-based progressive alignment).
18. ClustalW 1.83 source, `pairalign.c` (pairwise distance, https://raw.githubusercontent.com/epruesse/ARB/master/GDE/CLUSTALW/pairalign.c, opened); T-Coffee sources `lib/dp_lib/evaluate.c` (`residue_pair_extended_list`) and `lib/util_lib/aln_convertion_util.c` (`get_seq_sim`), https://github.com/cbcrg/tcoffee (raw files opened); T-Coffee paper PDF https://tcoffee.org/Publications/Ps_pdf/tcoffee.pdf (opened).
