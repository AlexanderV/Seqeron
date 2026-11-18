# Multiple Sequence Alignment (MSA)

| Field | Value |
|-------|-------|
| Algorithm Group | Alignment |
| Test Unit ID | ALIGN-MULTI-001 |
| Related Projects | Seqeron.Genomics |
| Implementation Status | Production (linear-gap model by default, additive affine-gap model `MsaGapModel.Affine`; declared deviations in §5.3) |
| Last Reviewed | 2026-11-04 (campaign 2026-09, batch B09) |

## 1. Overview

Multiple sequence alignment aligns three or more sequences in one shared coordinate system by inserting gap characters until every aligned sequence has the same length. Exact optimal MSA is computationally expensive, so practical implementations usually rely on heuristics. This repository provides five aligners over `DnaSequence` inputs, all in the **linear gap model** of `GlobalAlign` (every gap position scores `GapExtend`; `GapOpen` is not used):

1. `SequenceAligner.MultipleAlign(...)` — **center-star** alignment (Gusfield 1993) with *exact* pairwise alignments: all pair scores come from the canonical Needleman-Wunsch (`GlobalAlign` core), the center is the sequence with the highest total pair score, every other sequence is aligned to it with the same DP and the pairwise alignments are merged column-wise ("once a gap, always a gap"). For two sequences the result is an optimal pairwise global alignment.
2. `SequenceAligner.MultipleAlignAnchored(...)` — the approximate, faster variant of (1): 4-mer-cosine center (counts from the Core `CountKmersSpan`), suffix-tree exact-match anchors (`AnchorBasedAligner`: MUMmer-style seeds, weighted colinear chain with overlap trimming, Needleman-Wunsch on the gaps, exact Hirschberg for very large gaps). Pairwise alignments are *not* guaranteed optimal.
3. `SequenceAligner.MultipleAlignProgressive(...)` — guide-tree progressive alignment: ClustalW pairwise distances (1 − identities / shorter length) → guide tree: UPGMA by default (Higgins & Sharp 1988 CLUSTAL; Feng & Doolittle 1987 progressive scheme) or ClustalW's neighbour-joining tree (`MsaOptions.GuideTree = NeighborJoining`, §4.5) → profile-profile Needleman-Wunsch, single pass.
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

**Anchor-chain algorithm (audit A1-16).** `ChainAnchors` evaluates best[i] = Lᵢ + max(0, maxⱼ<ᵢ best[j] − trim(j,i)), trim(j,i) = max(0, endRef(j) − refStart(i), endQuery(j) − queryStart(i)), over the predecessors that leave a positive remainder (trim < Lᵢ ⇔ endRef(j) < endRef(i) ∧ endQuery(j) < endQuery(i)), with range-maximum queries instead of the all-pairs scan (the sparse dynamic programming of Abouelhoda & Ohlebusch 2005). With the diagonal δ = query − ref, the predecessors split into δⱼ ≤ δᵢ (reference overlap is the larger one, trim = (endRef(j) − refStart(i))⁺) and δⱼ > δᵢ (trim = (endQuery(j) − queryStart(i))⁺); inside each half trim is a function of ONE end coordinate (0 when it is ≤ the new anchor's start), so each half is two orthogonal range-maximum queries — a Fenwick tree over diagonal ranks whose nodes carry max segment trees over the end coordinate. Time O(a log² a), memory O(a log a); 10⁵ anchors chain in well under a second. The chain, its score and every tie decision (earliest predecessor, earliest anchor) are identical to the former O(a²) scan, which is kept as the oracle of `AnchorBasedAligner_ChainEquivalence_Tests`.

**Progressive distance.** `1 − identical aligned residue pairs / min(|a|, |b|)` of the optimal pairwise alignment — the ClustalW distance (`pairalign.c`, clustalw 1.83: `tmat = (100 − 100·count/MIN(len1,len2)) / 100`). Equal averaged distances tie (within 1e-9) and the lowest index pair merges first. `MsaOptions.Distance = MsaDistanceMode.ClustalWAffine` (MCP `distance=clustalw`) takes the identity count from ClustalW 2.1's own pairwise stage instead (`SequenceAligner.ClustalWPairwiseAlign`, §4.4); `Linear` (default) is the repository's linear NW alignment.

**T-Coffee library.** Primary weight of every aligned residue pair of a pairwise alignment = ⌊100 · identical columns / aligned (non-gap) columns⌋ — `get_seq_sim`, mode sim1, `lib/util_lib/aln_convertion_util.c` of the T-Coffee sources; global (NW) and local (SW) entries of the same residue pair are signal-added. Extended weight of any pair (Si.p, Sj.q) = direct weight + Σ over residues Sk.r of third sequences with both legs in the primary library of min(W1, W2) (`residue_pair_extended_list`, `lib/dp_lib/evaluate.c`); a pair without a direct entry but with triplet support therefore has a positive extended weight.

**ClustalW pairwise stage (§4.4 reference, `ClustalWPairwiseAlign`).** Port of `FullPairwiseAlign::pairwiseAlign` (ClustalW 2.1; the full-DP mode of `clustalw -ALIGN`/`-TREE`): an affine **local** Smith-Waterman-Gotoh alignment (`forwardPass`: end cell = first row-major maximum, scores floored at 0; `reversePass`: start cell; `diff`: Myers-Miller divide-and-conquer with ClustalW's midpoint tie rule) with the DNA matrix IUB × 100 (`matrices.h` `swgapdnamt`: identical codes +1000, non-compatible −900; the 2.1 transition-weight lines address residue codes A-E, C-M, C-N, so A/C/G/T pairs are unchanged) and gap cost 3000 + 666·k (`PWGAPOPEN 15` → `2·15·100`, `PWGAPEXT 6.66` → `666`). Percent identity = 100 · identical aligned code pairs / min(len1, len2) in single precision; distance = (100 − that)/100; an empty sequence gives distance 1; `-ALIGN` prints `(int)` of the percent identity as "Score:". Because the Myers-Miller choice among co-optimal paths decides the identity count, the alignment is reproduced step by step (not `LocalAlignAffine`). Verified on 3 600 pairs against the installed `clustalw` 2.1.

**Affine-gap MSA (`MsaOptions.GapModel = MsaGapModel.Affine`, MCP `gapModel=affine`).** *Objective (exact):* the sum of pairs with **natural** affine gap costs (Altschul 1989, J Theor Biol 138:297): for every pair of rows delete the columns that are gaps in both rows, score residue pairs with Match/Mismatch and charge every maximal run of gap columns of the projection `GapOpen + length·GapExtend` (a gap in row i directly followed by a gap in row j = two runs). `SequenceAligner.SumOfPairsScore(rows, scoring, MsaGapModel.Affine)` computes it in O(k²·L) and it is the `TotalScore` of every method under the affine model; for k = 2 it is the Gotoh score of `GlobalAlignAffine`, so a pairwise optimum is an MSA optimum. *Optimisation:* (a) star methods (`MultipleAlign`, anchored): every sequence is aligned to the center with the exact affine pairwise aligner (`GlobalAlignAffineLinearSpace`, identical to `GlobalAlignAffine`; anchored: `GlobalAlignAffineLinearSpace` between anchors), the center maximises the total optimal affine pair score; (b) progressive / iterative: a Gotoh three-state profile–profile DP (`AlignProfilesAffine`): column pair M(i, j) = Σ over cross pairs (Match/Mismatch for residue pairs, GapExtend for a residue against an existing gap, 0 gap against gap); a run of new all-gap columns against column i costs `r_i·|B|·(GapOpen + GapExtend)` at its first column and `r_i·|B|·GapExtend` afterwards (r_i = residues of column i; ClustalW/MUSCLE-style occupancy scaling); exact integer arithmetic, same Flouri border and M > own gap state > other gap state tie rule as `GlobalAlignAffine`, so for two sequences the rows equal `GlobalAlignAffine`. The iterative method accepts a re-alignment only if the *exact* natural objective strictly improves; (c) consistency: the T-Coffee primary library is built from Gotoh global and Smith–Waterman–Gotoh local alignments, the progressive DP keeps T-Coffee's zero gap penalty; guide-tree identities come from `GlobalAlignAffine` when the distance mode is `Linear`. *Optimal vs heuristic:* exact only for k = 2 (and the exact objective value for any rows). For k ≥ 3 the natural cost of a pair depends on the whole column history of both rows, which a three-state profile DP cannot track, so (a)–(c) are heuristics; on 227 random 3-sequence sets of length ≤ 6 (exact optimum by dynamic programming over (i, j, l, run state of each pair), validated against exhaustive enumeration) the affine star reaches the optimum in 170, progressive and iterative in 189 sets and none exceeds it (`SequenceAligner_MultipleAlignAffine_Tests`).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Exact MSA dynamic programming (theoretical baseline) | `O(L^n)` | `O(L^n)` | Complexity quoted in the current document from Wang and Jiang (1994). |
| `MultipleAlign` (exact center-star) | `O(k^2 * L^2)` pair scores + `O(k * L^2)` center alignments | `O(L^2)` per concurrent pair | Canonical NW per pair, parallel over pairs. |
| `MultipleAlignAnchored` center selection | `O(k^2 * 4^4 + k * L)` | `O(k * min(L, 256))` | `SelectCenterSequence` counts 4-mers per sequence (sparse dictionaries). |
| `MultipleAlignAnchored` pairwise phase | `O(L)` suffix tree, `O(k * L)` anchor streaming, `O(a \log^2 a)` chaining over `a` anchors (A1-16; was `O(a^2)`), `O(\sum_i \delta_i^2)` gap DP | `O(L)` + `O(\delta^2)` (Hirschberg above 25 M cells) | `AnchorBasedAligner` documents this split. |
| Iterative refinement (`MultipleAlignIterative`) | Progressive seed `O(k^2 L^2)` + refinement `O(I * E * (k * L^2))` where `E = O(k)` internal edges, `I = maxIterations` passes | `O(k * L)` | Each refinement re-alignment is one profile-profile NW (`O(k L^2)` over the column DP); passes stop early on convergence. Measured baseline: the full ALIGN-MULTI-001 iterative fixture (15 tests incl. a 500-trial random property test, k≤5, L≤8) completes in well under 1 s on the dev machine (2026-06-23). |
| Consistency-based (`MultipleAlignConsistency`) | Primary library + extension `O(k^3 L^2)` worst case; progressive DP `O(k^2 L^2)` | `O(k^2 L^2)` for the library | Library extension is the dominant term (Notredame et al. 2000 quote `O(N^3 L^2)` worst case, closer to `O(N^3 L)` in practice). Measured baseline: the full ALIGN-MULTI-001 consistency fixture (12 tests incl. a 100-trial random property test, k≤5, L≤8) completes in well under 1 s on the dev machine (2026-06-23). |

### 4.5 Guide trees (audit A1-12)

`MsaOptions.GuideTree` selects the tree that orders the merges of `MultipleAlignProgressive`, `MultipleAlignIterative` (also the edges of its refinement) and `MultipleAlignConsistency`:

- `Upgma` (default): `PhylogeneticAnalyzer.BuildTreeFromMatrix(…, TreeMethod.UPGMA)` with `TreeBuildOptions { TieTolerance = 1e-9 }` — proportional (cluster-size weighted) averaging, Sokal & Michener 1958; identical to the former private copy on 900 random and tie-heavy matrices.
- `NeighborJoining`: Saitou & Nei 1987 as ClustalW 2.1 uses it (`-CLUSTERING=NJ`, default; T-Coffee/Clustal-style pipelines). NJ is unrooted (trifurcating centre), so ClustalW (1) writes the unrooted tree (`.dnd`: 5 decimals, |b| < 1e-4 → 0, children ordered by lowest leaf index), (2) re-reads it with single-precision lengths and (3) roots it with `Tree::reRoot` (the branch giving the shallowest tree whose root-to-leaf mean difference is smallest) before `Tree::createSets` lists the merges in post-order, left subtree first. The canonical NJ builder (with `TieOrder = SecondIndexMajor` = ClustalW's `NJTree::generateTree` loop nest) supplies the unrooted tree; `SequenceAligner.GuideTree.cs` ports steps (1)-(3) line by line (`createTree`/`insertNode`/`reRoot`/`calcMean`/`insertRoot`/`groupSeqs`, clustalw 2.1+lgpl-7 `src/tree/Tree.cpp`, `ClusterTreeOutput.cpp`). With `MsaDistanceMode.ClustalWAffine` the distances are ClustalW's single-precision `(100f − 100f·count/min)/100f`.

Verification: 420 DNA sets (families, two clades, duplicate sequences, two-letter alphabets, tiny sequences; 2-11 sequences) — the unrooted tree equals the `.dnd` of the installed `clustalw` 2.1 (`-ALIGN -MAXDIV=0 -NEWTREE`; edge for edge, branch lengths ≤ 1.5·10⁻⁵) and the sizes of the progressive merges equal the "Group n: Sequences: s" lines (0 mismatches; the ordinary first-index tie order fails 21 sets). The canonical NJ equals Biopython 1.88 `DistanceTreeConstructor().nj` and the canonical UPGMA equals scipy `linkage('average')` on 400 random matrices each (Biopython's own `upgma()` averages without cluster sizes — WPGMA — and is not the reference).

### 4.6 ClustalW 2.1 multiple-alignment stage (audit A1-10b / A1-15b)

`SequenceAligner.MultipleAlignClustalW(sequences, ClustalWOptions?, ScoringMatrix?)` is a step-by-step port of `clustalw -TYPE=DNA -ALIGN` (ClustalW 2.1+lgpl-7 sources: `Clustal.cpp`, `tree/Tree.cpp`, `tree/NJTree.cpp`, `tree/TreeInterface.cpp`, `multipleAlign/{MSA,MyersMillerProfileAlign,ProfileBase,ProfileStandard,ProfileWithSub}.cpp`, `substitutionMatrix/SubMatrix.cpp`, `pairwise/FastPairwiseAlign.cpp`; Thompson, Higgins & Gibson 1994; Larkin et al. 2007). It is a separate method (not an `MsaOptions` switch of the Seqeron profile aligners) because the result is defined by ClustalW's own arithmetic, not by a scoring model:

1. **Pairwise distances** — the full DP stage (§4.4, `ClustalWOptions.PairwiseMode = Full`, `-PWGAPOPEN/-PWGAPEXT`, `-TRANSWEIGHT` also moves the C-M / C-N scores of IUB codes) or the **fast k-tuple stage** `ClustalWQuickPairwiseAlign` (`PairwiseMode = Quick`, A1-15b, below), stored as in ClustalW's `DistMatrix` (double; the full stage computes `(100 - 100·count/MIN(len1,len2))/100` in single precision).
2. **Guide tree** — ClustalW's neighbour-joining arithmetic (`ClustalWNeighborJoining`, a port of `NJTree::generateTree`: height-corrected matrix, update `(d_ik + d_jk)/2`, pair = `fnseqs2·d_ij + 2·sumd − diq − djq` in a jj-outer / ii-inner scan with strict `<`), 5-decimal print, single-precision re-read and `Tree::reRoot` (A1-12). The exact arithmetic matters for tied distances (fast-mode scores and duplicate sequences tie often): the canonical `PhylogeneticAnalyzer` NJ differs from ClustalW on such ties by rounding noise and gave a different topology on 2 of 700 oracle sets once the fast mode was included, so the ClustalW guide tree uses the ported arithmetic (the canonical builder stays the oracle-checked implementation of Saitou-Nei NJ for everything else).
3. **Sequence weights** (`Tree::calcSeqWeights`, Thompson 1994): w_i = ⌊100·Σ_branches on the path leaf→root (branch length / number of leaves below it)⌋ (single precision, root edge excluded), normalised to `INT_SCALE_FACTOR` = 1000 (Unix build) with a minimum of 1; `-NOWEIGHTS` = equal weights 1000/k. The pairwise identities used from here on come from the TREE (`Tree::calcSimilarities`: 100 − 100·(path length), clamped to [0.01, 1]), not from the pairwise stage.
4. **Progressive order and delayed divergent sequences** (`MSA::multiSeqAlign`): the internal nodes of the rooted tree in post-order (left subtree first) are the alignment steps. A sequence whose best (tree-derived) identity to any other sequence is ≤ `-MAXDIV` (30 %) is skipped in this first pass; afterwards the remaining sequences are added one at a time, the closest to the aligned set first, with weights = tree weight × identity to the new sequence, renormalised.
5. **Profile alignment** (`MyersMillerProfileAlign::profileAlign`, `ProfileBase::calcGapCoeff`, `ProfileWithSub`, `ProfileStandard`): the larger group is profile 1; profile 1 holds, for every column, `⌊(Σ_aa weight_aa·M[aa][res]/Σ weights)·(non-gap fraction)⌋` for every residue, profile 2 holds `⌊10·weight fraction⌋`; the column score is Σ p1·p2 / 10. **Position-specific gap penalties**: gap open = ⌊100·GAPOPEN⌋ = 1500, extension = ⌊100.0·GAPEXT⌋ with GAPEXT a float (6.66f → 665, not 666); a column that already holds gaps (end gaps do not count) gets ⌊(non-gap fraction)·open⌋ and ⌊0.5·extension⌋ (so alignments open gaps where gaps already are); every penalty is at least 1; the penalties at both ends of a profile are 0 (end gaps are free). **Matrix**: IUB × 100 shifted to be non-negative (−(lowest score) = +900 on every residue pair unless `-NEGATIVE`), then A-G, C-T, C-U := ⌊`TRANSWEIGHT`·M[A][A]⌋ (0.5 → 950). The two profiles are aligned by Myers & Miller (1988) with ClustalW's midpoint rule (`progDiff`).
6. For **DNA** ClustalW's residue-specific (Pascarella & Argos), hydrophilic-run and variable-region gap penalties and the gap-separation distance (`-GAPDIST`, `-NOPGAP`, `-NOHGAP`, `-NOVGAP`) are guarded by `if (!_DNAFlag)` in `calcGapCoeff`; they do not exist for DNA and are therefore not options of `ClustalWOptions` (they belong to the protein extension, A1-14).

**Fast k-tuple pairwise stage (`ClustalWQuickPairwiseAlign`, `MsaDistanceMode.ClustalWQuick`; A1-15b).** Port of `FastPairwiseAlign::pairAlign` (Wilbur & Lipman 1983): every k-mer (k = `-KTUPLE`, DNA 1-4, default 2; A=0, C=1, G=2, T=U=3, any IUB code breaks the k-mer) of sequence 1 is matched against the position list of sequence 2; the k-mer hits are counted per diagonal; the `-TOPDIAGS` (4) diagonals with the highest counts and the `-WINDOW` (4) diagonals around each are kept (the quicksort is ClustalW's own, unstable and leaves the last diagonal outside its range, so it is ported verbatim); a chain of exactly matching fragments on those diagonals is scored with gap penalty `-PAIRGAP` (5; the command line sets `ktuple + 4` whenever `-KTUPLE` is given, so `-KTUPLE=2` alone means 6); the best chain score S gives the percent 100·S/MIN(len1, len2) (printed as `Score:` with `%lg`, `-SCORE=ABSOLUTE` = S) and the distance (100 − percent)/100 in double precision. When the fragment list reaches 2·maxAlnLength (= 4 × the longest sequence of the whole input) ClustalW stops with "(Partial alignment)" and keeps the best chain so far; the port reproduces this, which is why the method takes the input's longest length.

**Verification.** `clustalw_quick_oracle.tsv`: 8 257 pairs of 605 sets (`clustalw -QUICKTREE -ALIGN`, k = 1-4, window / topdiags / pairgap / score variants, IUB-N sets, repeat-rich sets that reach the partial-alignment exit) — every printed `Score:` equal. `clustalw_msa_full_oracle.tsv`: 1 000 sets (2-25 sequences, ≤ 600 nt; options -GAPOPEN -GAPEXT -TRANSWEIGHT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE -KTUPLE -WINDOW -TOPDIAGS -PAIRGAP -PWGAPOPEN -PWGAPEXT; IUB-ambiguity, duplicate and tiny sets) — every output row of `clustalw -TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` identical; a development panel of 2 500 further sets (another seed) was identical after the fixes listed in F56. Because the binary runs its own pairwise stage, tree and weights, this one comparison covers the whole chain; `-USETREE` with the same tree cannot add information.

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAligner.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAligner.cs)

- `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`: exact center-star (Gusfield 1993); public MSA entry point.
- `SequenceAligner.MultipleAlignAnchored(IEnumerable<DnaSequence>, ScoringMatrix?)`: approximate suffix-tree-anchored center-star.
- `SequenceAligner.SelectCenterByPairwiseScore(...)` / `SelectCenterSequence(...)`: Gusfield pair-score center / 4-mer cosine center (anchored variant; k-mers from the Core `CountKmersSpan`).
- `SequenceAligner.MergeCenterAlignments(...)`: merges the center-vs-other pairwise results into one alignment.
- `SequenceAligner.BuildConsensus(...)`, `ComputeSumOfPairsScore(...)`: consensus and public `TotalScore`.
- `SequenceAligner.MultipleAlignProgressive(IEnumerable<DnaSequence>, ScoringMatrix?)`: ClustalW-distance / UPGMA (or NJ, `MsaOptions.GuideTree`) / profile-profile NW progressive MSA.
- `SequenceAligner.BuildProgressiveGuideTree` / `BuildClustalWNeighborJoiningGuideTree` / `ClustalWNeighborJoining` (`SequenceAligner.GuideTree.cs`): the canonical `PhylogeneticAnalyzer.BuildTreeFromMatrix` UPGMA tree, or ClustalW's own NJ arithmetic (port of `NJTree::generateTree`, F56), adapted to the binary rooted guide tree (§4.5).
- `SequenceAligner.MultipleAlign/MultipleAlignAnchored/MultipleAlignProgressive/Iterative/Consistency(IEnumerable<DnaSequence>, ScoringMatrix?, MsaOptions?)`: the same with additive options (`MsaOptions.GapModel` linear/affine, `MsaOptions.Distance` guide-tree distance mode for the progressive family).
- `SequenceAligner.MultipleAlignClustalW(sequences, ClustalWOptions?, ScoringMatrix?)` (`SequenceAligner.ClustalWMultiple.cs`): ClustalW 2.1 multiple alignment of DNA, row-for-row equal to `clustalw -ALIGN` (§4.6, audit A1-10b); `SequenceAligner.ClustalWQuickPairwiseAlign(a, b, ClustalWQuickOptions?, maxSequenceLength?)` (`SequenceAligner.ClustalWQuick.cs`): the fast k-tuple pairwise stage (`-QUICKTREE`), also `MsaDistanceMode.ClustalWQuick` / MCP `distance=clustalw-quick` (audit A1-15b).
- `SequenceAligner.SumOfPairsScore(rows, scoring, MsaGapModel)`: the sum-of-pairs objective of any alignment under either gap model (§4.2).
- `SequenceAligner.ClustalWPairwiseAlign(string | DnaSequence, …, gapOpen = 15, gapExtend = 6.66)`: ClustalW 2.1 pairwise stage with its percent identity (§4.4).
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
- Note on the refinement variant: `MultipleAlignIterative` is the MUSCLE Stage 3 scheme as published (Edgar 2004), which partitions by guide-tree edge; the earlier Barton-Sternberg (1987) scheme removes one sequence at a time and re-aligns it to the profile of the rest. They are two published refinement schemes, so this is a choice of published method, not a deviation of the MUSCLE implementation.
- T-Coffee (Notredame, Higgins & Heringa 2000, `MultipleAlignConsistency`): primary library of global (NW) and local (SW) pairwise alignments weighted by percent identity (aligned columns only, integer division — T-Coffee `get_seq_sim` sim1) and combined by signal addition (p.207); library extension by the triplet consistency transformation for every residue pair, extended weight = direct weight + Σ min(W1, W2) (p.209, worked example 88 → 165; `residue_pair_extended_list`); progressive alignment with the extended-library weights as position-specific substitution scores and zero gap penalties (p.210).

**Intentionally simplified / declared deviations (each with its consequence):**

- **Gap model.** The default stays the linear model (every gap position costs `GapExtend`, `GapOpen` unused, column sum of pairs). The affine model is the additive option `MsaGapModel.Affine` (§4.2): exact natural affine objective, exact for two sequences, heuristic optimisation for k ≥ 3 (position-specific profile DP; the natural cost cannot be tracked by a three-state DP). Not implemented: ClustalW's own multiple-alignment stage (`prfalign.c`: sequence weights, residue-specific and hydrophilic gap penalties, gap-separation distance, delayed divergent sequences) and MUSCLE's profile scoring (audit A1-10b); the exact optimum for three or four sequences (audit A1-11).
- The star methods build no guide tree; the progressive family uses UPGMA by default (Higgins & Sharp 1988 CLUSTAL; the historical behaviour) and neighbour-joining as ClustalW/T-Coffee on request (`MsaOptions.GuideTree = NeighborJoining`, A1-12, §4.5). The UPGMA tree comes from the canonical `PhylogeneticAnalyzer` builder; since A1-10b the NJ tree uses ClustalW's own NJ arithmetic (`ClustalWNeighborJoining`, same topology as Saitou-Nei, ties resolved as clustalw does) instead of the canonical NJ plus tie order (the private UPGMA copy and the "circular reference" it was justified by are gone). With `NeighborJoining` + `ClustalWAffine` the unrooted tree equals the `.dnd` file of `clustalw -ALIGN` and the progressive merge order equals ClustalW's (420 DNA sets); the final rows are not claimed to equal ClustalW's, because its profile stage (sequence weights, position-specific gap penalties, delayed divergent sequences) is audit item A1-10b.
- The default guide-tree distance is computed from the repository's linear NW alignment; `MsaDistanceMode.ClustalWAffine` reproduces ClustalW's own pairwise stage (full-DP mode, exact; §4.4) and `MsaDistanceMode.ClustalWQuick` its fast k-tuple mode (`-QUICKTREE`, exact; §4.6, audit A1-15b). ClustalW's protein pairwise stage is not implemented (audit A1-15b part (b), needs the protein MSA of A1-14). The Seqeron profile aligners (`MultipleAlignProgressive/Iterative`) remain Seqeron's own occupancy-scaled models; ClustalW's own profile stage is the separate, exact `MultipleAlignClustalW` (§4.6).
- The T-Coffee primary library uses the single best local alignment per pair and no sequence-redundancy weights (`1/Σ sim³`). The non-intersecting top-N local alignments of Lalign are available as `SequenceAligner.LocalAlignTopN` (Waterman-Eggert declumping, [Local_Alignment_Smith_Waterman.md §5.5](Local_Alignment_Smith_Waterman.md), audit A1-7) but are not yet fed into the library (audit A1-14).
- `MultipleAlignAnchored` is a heuristic: anchor chaining constrains the path, so pairwise alignments (and the SP score) can fall below `MultipleAlign`'s.

**Not implemented:**

- Exact optimal MSA (NP-complete; Wang & Jiang 1994).
- MUSCLE's profile score (Edgar 2004 `profile.cpp`; audit A1-10c) and ClustalW's protein-only gap penalties / matrices / `-CLUSTERING=UPGMA` tree / iteration (A1-14, A1-10c): `MultipleAlignClustalW` is DNA only (`-DNAMATRIX=IUB`).

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
| Tied neighbour-joining pairs (e.g. identical sequences) | The pair with the smallest second index, then the smallest first index, joins first (ClustalW's `NJTree` loop order, equality within 1e-9). | `TreeBuildOptions { TieOrder = SecondIndexMajor }`. |

### 6.2 Limitations

The star result depends on the chosen center and is not SP-optimal (Gusfield's 2 − 2/k bound is for metric distances). The progressive method is single-pass and inherits the "once a gap, always a gap" behaviour; `MultipleAlignIterative` removes that limitation by re-splitting and re-aligning, but it is SP-guided and restricted to guide-tree edge partitions, so it improves on but does not guarantee the global SP optimum. `MultipleAlignConsistency` optimises the T-Coffee consistency objective over a greedy progressive merge order, so it does not guarantee the global optimum of that objective. Memory: the exact center-star (`MultipleAlign`) aligns every sequence to the center with the exact linear-space pairwise program (identical rows to the full-matrix Needleman-Wunsch, [Linear_Space_Alignment.md](Linear_Space_Alignment.md)) and has no quadratic matrix; the progressive / iterative / consistency profile alignments keep an `O(m·n)` double matrix of profile columns and throw `ArgumentException` above `SequenceAligner.ProfileMatrixMaxCells` = 50 000 000 cells (400 MB) instead of an `OutOfMemoryException`. All aligners are documented and tested for `DnaSequence` inputs only, and the linear-gap deviation of §5.3 applies to all of them.

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
