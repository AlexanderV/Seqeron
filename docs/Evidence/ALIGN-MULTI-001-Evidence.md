# ALIGN-MULTI-001: Multiple Sequence Alignment - Evidence Document

## Test Unit
**ID:** ALIGN-MULTI-001
**Area:** Alignment
**Methods:** `SequenceAligner.MultipleAlign(IEnumerable<DnaSequence>, ScoringMatrix?)`

---

## Authoritative Sources

### Source 1: Wikipedia - Multiple Sequence Alignment
**URL:** https://en.wikipedia.org/wiki/Multiple_sequence_alignment

**Key Information Extracted:**
1. **Definition:** MSA is the process of aligning three or more biological sequences (protein, DNA, or RNA)
2. **Objective:** Insert gaps into sequences until all sequences conform to length L ≥ max{n_i}
3. **Constraint:** No column may consist only of gaps
4. **Complexity:** Finding optimal MSA is NP-complete (Wang & Jiang, 1994; Just, 2001; Elias, 2006)
5. **Progressive Alignment:** Most widely used heuristic approach developed by Feng & Doolittle (1987)
   - Stage 1: Create guide tree representing sequence relationships
   - Stage 2: Build MSA by adding sequences according to guide tree
6. **Sum-of-Pairs Score:** "sum of all of the pairs of characters at each position in the alignment (the so-called sum of pair score)" — common objective function for MSA optimization. Sums scores across all C(k,2) sequence pairs.
7. **Note:** Consensus sequence derivation is NOT described on the MSA page. See Source 3 (Wikipedia Consensus sequence).
8. **Error Propagation:** Progressive methods propagate early alignment errors to final result
9. **Reversibility:** "To return from each particular sequence S'_i to S_i, remove all gaps"
10. **Benchmark:** BAliBase evaluation found at least 24% of aligned amino acid pairs were incorrect (Nuin et al., 2006)

### Source 2: Wikipedia - Clustal
**URL:** https://en.wikipedia.org/wiki/Clustal

**Key Information Extracted:**
1. **Progressive Method:** Aligns sequences in most-to-least similarity order
2. **Guide Tree:** Can be generated via UPGMA or neighbor-joining
3. **Pairwise Alignment:** First step computes distance matrix between all pairs
4. **Gap Penalties:** "The gap opening penalty and gap extension penalty parameters can be adjusted by the user."
5. **Output Symbols:**
   - `*` asterisk: fully conserved position
   - `:` colon: conservation between strongly similar groups
   - `.` period: conservation between weakly similar groups
   - ` ` blank: non-conserved
6. **k-tuple similarity:** "the number of k-tuple matches between two sequences, accounting for a set penalty for gaps" (Clustal/ClustalV)

### Source 3: Wikipedia - Consensus sequence
**URL:** https://en.wikipedia.org/wiki/Consensus_sequence

**Key Information Extracted:**
1. **Consensus definition:** "the calculated sequence of most frequent residues, either nucleotide or amino acid, found at each position in a sequence alignment"
2. **Majority voting:** Consensus = most frequent residue at each column position
3. **Note:** Does not specify whether gaps participate in the vote or how ties are resolved

---

## Implementation Details

The `SequenceAligner.MultipleAlign()` method implements an **anchor-based star alignment**:

1. **Center selection** via 4-mer cosine similarity — selects the sequence with highest total similarity to all others (O(k²·L)). *Note: ClustalV uses k-tuple match counts (Wikipedia Clustal); cosine similarity on k-mer frequency vectors is an implementation design choice.*
2. **Suffix tree construction** on center sequence (O(L)) — *implementation-specific; not from external sources*
3. **Anchor-based pairwise alignment** of each other sequence to center via suffix tree exact-match anchors, with Needleman-Wunsch for gaps between anchors — *implementation-specific*
4. **Gap reconciliation** — merges gap columns from independent pairwise alignments into a single MSA coordinate space — *implementation-specific*
5. **Pad sequences** to equal length
6. **Consensus generation** via majority voting at each position (Wikipedia Consensus sequence: "most frequent residues... at each position"). *Implementation choices: gaps participate in the vote; on tie between gap and nucleotide, nucleotide is preferred.*
7. **Sum-of-pairs scoring** — column-based SP score across all C(k,2) pairs (Wikipedia MSA). *Implementation choice: gap-gap pairs scored as 0 (standard bioinformatics convention, not stated in Wikipedia).*

### Complexity
- **Time:** O(k² × m) where k = number of sequences, m = average sequence length
- **Space:** O(k × L) for aligned sequences, where L = alignment length

### Scoring
- **TotalScore** is the true sum-of-pairs (SP) score per Wikipedia MSA: "sum of all of the pairs of characters at each position in the alignment"
- Column-based: for each position, for each pair (i,j): match → +Match, mismatch → +Mismatch, gap-nucleotide → +GapExtend, gap-gap → 0
- Sums all C(k,2) pair scores
- *Note: gap-gap = 0 is standard bioinformatics convention, not explicitly stated in Wikipedia MSA*

### Consensus
- **Majority voting** at each column, per Wikipedia Consensus sequence: "the calculated sequence of most frequent residues... found at each position in a sequence alignment"
- All characters (including gaps) participate in the vote — *implementation design choice, not specified in Wikipedia*
- On tie between gap and nucleotide, nucleotide is preferred — *implementation design choice, not specified in Wikipedia*

---

## Test Datasets from Sources

### From Wikipedia MSA Definition:
1. **Minimum sequences:** MSA requires 3+ sequences (pairwise alignment handles 2)
2. **Empty set:** Edge case returning empty result
3. **Single sequence:** Trivial case returning the sequence unchanged
4. **Identical sequences:** Should produce perfect alignment
5. **Variable length sequences:** Requires gap insertion

### From ClustalW/Progressive Alignment:
1. **Highly similar sequences:** Expected to align well
2. **Divergent sequences:** Progressive methods may produce suboptimal results
3. **Gap-heavy sequences:** Tests gap padding behavior

---

## Edge Cases and Corner Cases

| Case | Expected Behavior | Source |
|------|-------------------|--------|
| Empty input | Return `MultipleAlignmentResult.Empty` | Wikipedia MSA (edge case) |
| Single sequence | Return alignment containing just that sequence | Wikipedia MSA (trivial case) |
| Null input | Throw `ArgumentNullException` | .NET conventions |
| Two identical sequences | Perfect alignment, consensus equals input | Wikipedia MSA |
| All identical sequences | All outputs equal input, consensus equals input | Wikipedia MSA |
| Sequences of different lengths | All aligned sequences equal length (padded) | Wikipedia MSA: "conform to length L ≥ max{n_i}" |
| All-gap column | Must not occur per MSA invariant | Wikipedia MSA: "no column consists of only gaps" |

---

## Invariants (from sources)

1. **Length invariant:** All aligned sequences have equal length (Wikipedia MSA)
2. **Gap-only column prohibition:** No column consists only of gaps (Wikipedia MSA)
3. **Reversibility:** Removing gaps from each aligned sequence recovers the original (Wikipedia MSA)
4. **Count preservation:** Number of aligned sequences equals number of input sequences (Wikipedia MSA)
5. **Consensus validity:** Consensus characters come from {A, C, G, T, -} (Wikipedia MSA)
6. **Sum-of-pairs correctness:** TotalScore = Σ over all pairs (i<j) of column-based score (Wikipedia MSA)

---

## Known Limitations (from sources)

1. **Error propagation:** Early alignment errors propagate through progressive alignment (Wikipedia MSA)
2. **Center dependence:** Result depends on which sequence is chosen as center (star alignment property)
3. **Distantly related sequences:** Performance degrades with evolutionary distance (Wikipedia MSA)
4. **Suboptimality:** Heuristic methods do not guarantee optimal solution (Wikipedia MSA: NP-complete)
5. **No guide tree:** Simplified star alignment lacks the phylogenetic guide tree used by ClustalW (Wikipedia Clustal)

---

## References

1. Feng DF, Doolittle RF (1987). "Progressive sequence alignment as a prerequisite to correct phylogenetic trees." J Mol Evol. 25(4):351-360.
2. Wang L, Jiang T (1994). "On the complexity of multiple sequence alignment." J Comput Biol. 1(4):337-348.
3. Nuin PA, Wang Z, Tillier ER (2006). "The accuracy of several multiple sequence alignment programs for proteins." BMC Bioinformatics. 7:471.
4. Thompson JD, Higgins DG, Gibson TJ (1994). "CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment." Nucleic Acids Res. 22(22):4673-80.
5. Wikipedia. "Consensus sequence." https://en.wikipedia.org/wiki/Consensus_sequence

---

# Addendum: Iterative Refinement (MultipleAlignIterative) — collected 2026-06-23

This addendum covers the **iterative refinement** aligner `SequenceAligner.MultipleAlignIterative`,
which removes the single-pass "once a gap, always a gap" limitation of the progressive method by
re-splitting the alignment along guide-tree edges and re-aligning the two sub-profiles, keeping the
result only on a non-decreasing sum-of-pairs (SP) score.

## Online Sources (retrieved this session)

### Source A — Edgar RC (2004), MUSCLE (PRIMARY, authority rank 1)

**URL:** https://academic.oup.com/nar/article/32/5/1792/2380623
**Retrieved by:** WebSearch "MUSCLE Edgar 2004 tree-dependent restricted partitioning iterative
refinement sum of pairs SP score profile realignment" → WebFetch of the article page.
**Accessed:** 2026-06-23   **Authority rank:** 1 (peer-reviewed, Nucleic Acids Research)

**Key Extracted Points (verbatim from the fetched Stage 3 description):**

1. **Edge selection (3.1):** "An edge is chosen from TREE2 (edges are visited in order of
   decreasing distance from the root). TREE2 is divided into two subtrees by deleting the edge."
2. **Profile + realign (3.2/3.3):** "The profile of the multiple alignment in each subtree is
   computed. A new multiple alignment is produced by re-aligning the two profiles."
3. **Accept on SP improvement (3.4):** "If the SP score is improved, the new alignment is kept,
   otherwise it is discarded."
4. **Convergence / cap:** "Steps 3.1–3.4 are repeated until convergence or until a user-defined
   limit is reached."
5. **Name:** this is MUSCLE's Stage 3, "Refinement", a variant of *tree-dependent restricted
   partitioning*.

### Source B — Barton GJ & Sternberg MJ (1987) (PRIMARY, authority rank 1)

**URL:** https://pubmed.ncbi.nlm.nih.gov/3430611/
**Retrieved by:** WebSearch "Barton Sternberg 1987 ... iterative refinement remove sequence
realign J Mol Biol 198 327" → WebFetch of the PubMed abstract.
**Accessed:** 2026-06-23   **Authority rank:** 1 (peer-reviewed, J Mol Biol)

**Key Extracted Points:**

1. **Progressive seed:** "Initially, two sequences are aligned, then the third sequence is aligned
   against the alignment of both sequences one and two ... repeated until all sequences have been
   aligned."
2. **Iteration:** "Iteration is then performed to yield a final alignment." (The abstract confirms
   iterative refinement of an existing alignment; the exact remove-one mechanism is in the body.)
3. **Citation:** Barton, G.J. and Sternberg, M.J. "A strategy for the rapid multiple alignment of
   protein sequences. Confidence levels from tertiary structure comparisons." J Mol Biol. 1987 Nov
   20;198(2):327-37.

### Source C — Wallace, O'Sullivan & Higgins (2005), "Evaluation of iterative alignment algorithms" (rank 1)

**URL:** https://academic.oup.com/bioinformatics/article/21/8/1408/249176
**Retrieved by:** WebSearch (MUSCLE / iterative refinement) → WebFetch.
**Accessed:** 2026-06-23   **Authority rank:** 1 (peer-reviewed, Bioinformatics)

**Key Extracted Points:**

1. **Barton-Sternberg "Remove First" scheme:** "In each iteration step a sequence is removed from
   the alignment and realigned to the remaining alignment."
2. **Iteration bound:** the procedure "continues until the alignment score converges or reaches a
   computational limit," confirming a convergence-or-cap stopping rule.

### Source D — Wikipedia, "Multiple sequence alignment" (rank 4; used for the SP-score definition)

**URL:** https://en.wikipedia.org/wiki/Multiple_sequence_alignment
**Retrieved by:** WebFetch.   **Accessed:** 2026-06-23   **Authority rank:** 4

**Key Extracted Points:**

1. **SP score definition:** the program "optimizes the sum of all of the pairs of characters at
   each position in the alignment (the so-called *sum of pair* score)."
2. **Iterative methods overcome single-pass error:** "iterative methods can return to previously
   calculated pairwise alignments or sub-MSAs incorporating subsets of the query sequence as a
   means of optimizing a general objective function," addressing that in progressive methods "once
   a sequence has been aligned into the MSA, its alignment is not considered further."

## Documented Corner Cases and Failure Modes (iterative refinement)

### From Source A (Edgar 2004)

1. **No improvement found:** if no edge improves the SP score, the seed is returned unchanged (the
   refinement is a strict-improvement filter). Implemented: convergence stops the loop.
2. **User-defined limit:** refinement passes are capped; convergence usually halts earlier.

## Test Datasets (iterative refinement)

### Dataset: Hand-derived gap-relocation correction

**Source:** Derived from the SP definition (Source D) under SimpleDna (match +1, mismatch −1,
gap −1); independently recomputed with a stand-alone SP function during test authoring.

| Parameter | Value |
|-----------|-------|
| Input | CGA, GAGAT, CGC, GAC |
| Progressive seed | `-CG-A` / `GAGAT` / `-CG-C` / `--GAC`, SP = −8 |
| Refined (iterative) | `-CGA-` / `GAGAT` / `-CGC-` / `--GAC`, SP = −6 |
| Improvement | +2 (gap relocated from an internal to a terminal column) |

## Assumptions (iterative refinement)

1. **ASSUMPTION: Edge-partition refinement (not single-sequence removal).** The implementation uses
   MUSCLE's guide-tree edge partitioning (Source A) rather than Barton-Sternberg's remove-one-
   sequence loop (Sources B/C). Both are accept-on-SP-improvement iterative refinement of the same
   seed; the edge scheme reuses the existing UPGMA guide tree and profile–profile NW. This is an
   API/structure choice, not a correctness-affecting parameter: the accept rule (non-decreasing SP)
   is identical and is the property the tests verify.
2. **ASSUMPTION: SP scoring conventions (gap-gap = 0, residue-gap = GapExtend).** Inherited from the
   existing `ComputeSumOfPairsScore`; the gap-gap = 0 convention is standard but not stated verbatim
   in Source D (already documented for the star/progressive aligners above).

## Recommendations for Test Coverage (iterative refinement)

1. **MUST Test:** refined SP ≥ progressive seed SP for every input (monotonicity) — Evidence: Source A step 3.4.
2. **MUST Test:** a constructed case where progressive misplaces a gap and refinement corrects it, with hand-derived SP — Evidence: Source D (SP) + hand derivation.
3. **MUST Test:** all rows equal length and degap to input; no all-gap column — Evidence: Wikipedia MSA invariants (main doc).
4. **MUST Test:** idempotence on an already-optimal alignment; convergence within the cap — Evidence: Source A convergence clause.
5. **MUST Test:** determinism across repeated runs (no RNG) — Evidence: implementation contract.

## References (iterative refinement addendum)

6. Edgar RC (2004). "MUSCLE: multiple sequence alignment with high accuracy and high throughput." Nucleic Acids Res. 32(5):1792-1797. https://academic.oup.com/nar/article/32/5/1792/2380623
7. Barton GJ, Sternberg MJ (1987). "A strategy for the rapid multiple alignment of protein sequences. Confidence levels from tertiary structure comparisons." J Mol Biol. 198(2):327-337. https://pubmed.ncbi.nlm.nih.gov/3430611/
8. Wallace IM, O'Sullivan O, Higgins DG (2005). "Evaluation of iterative alignment algorithms for multiple alignment." Bioinformatics. 21(8):1408-1414. https://academic.oup.com/bioinformatics/article/21/8/1408/249176
9. Wikipedia. "Multiple sequence alignment." https://en.wikipedia.org/wiki/Multiple_sequence_alignment
