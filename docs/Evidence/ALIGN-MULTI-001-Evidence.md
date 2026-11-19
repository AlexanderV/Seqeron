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

The `SequenceAligner.MultipleAlign()` method implements the **center-star alignment** (Gusfield 1993; corrected 2026-11, B09 F4 — the earlier anchor-based pairwise stage is now `MultipleAlignAnchored`):

1. **Pairwise scores** — the optimal Needleman-Wunsch score of every pair (canonical `GlobalAlign` core).
2. **Center selection** — the sequence with the highest total pair score (Gusfield: the string minimising the total pairwise distance). Ties: lowest index.
3. **Pairwise alignment** of every other sequence to the center with the canonical Needleman-Wunsch.
4. **Gap reconciliation** — "once a gap, always a gap" merge: per center position the maximum number of inserted columns over all rows; own insertions left-justified in the slot (Gusfield 1993; Feng & Doolittle 1987).
5. **Consensus generation** via majority voting at each position (Wikipedia Consensus sequence: "most frequent residues... at each position"). *Implementation choices: gaps participate in the vote; on tie between gap and nucleotide, nucleotide is preferred.*
6. **Sum-of-pairs scoring** — column-based SP score across all C(k,2) pairs (Wikipedia MSA). *Implementation choice: gap-gap pairs scored as 0 (standard bioinformatics convention, not stated in Wikipedia).*

`MultipleAlignAnchored` keeps the former implementation idea (4-mer cosine center, suffix-tree anchors per MUMmer/LAGAN, NW between anchors) as an explicitly approximate variant.

### Complexity
- **Time:** O(k² × m²) for the exact center-star (all pair scores); `MultipleAlignAnchored` O(k²·4⁴ + k·m + Σδ²)
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
4. **Suboptimality:** Heuristic methods do not guarantee optimal solution (Wikipedia MSA: NP-complete). Since 2026-11 (F60) the exact sum-of-pairs optimum is available for small inputs (`MultipleAlignExact`: k = 2, k = 3, k >= 4 by Carrillo-Lipman up to 8 linear / 6 affine sequences); beyond that the NP-completeness (Wang & Jiang 1994) stays.
5. **Guide tree (2026-06 text: 'no guide tree: simplified star' - superseded 2026-11):** the star methods build no tree by design (Gusfield center-star, exact for k = 2); the progressive / iterative / consistency family uses the canonical `PhylogeneticAnalyzer` UPGMA or neighbour-joining tree (`MsaOptions.GuideTree`, F53), and `MultipleAlignClustalW` / `MultipleAlignTCoffeeProgram` / `MultipleAlignMuscle` reproduce the guide trees of ClustalW 2.1, T-Coffee 13.41 and MUSCLE 3.8.1551 (F56, F66, F65).

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

---

# Addendum: Consistency-based MSA (T-Coffee) — `MultipleAlignConsistency`, collected 2026-06-23

This addendum covers the **consistency-based** aligner `SequenceAligner.MultipleAlignConsistency`,
which optimises the T-Coffee **consistency objective** (a distinct objective class from the
sum-of-pairs score the other methods optimise), via a primary library + library extension +
progressive alignment on the extended library.

## Online Sources (retrieved this session 2026-06-23)

### Source T1 — Notredame, Higgins & Heringa (2000), "T-Coffee" (PRIMARY, rank 1)

**URL:** https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf (full text of J Mol Biol 302:205–217; DOI 10.1006/jmbi.2000.4042; PubMed 10964570)
**Retrieved by:** WebSearch "Notredame Higgins Heringa 2000 T-Coffee … primary library extension" → WebFetch of the Stanford PDF → Read of PDF pages 3–8 (Figures 1–2, pp.207–211).
**Accessed:** 2026-06-23   **Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points (verbatim / close paraphrase from the read pages):**

1. **Pipeline (Figure 1, p.207):** ClustalW (global) and Lalign (local) primary libraries →
   *Weighting / Signal Addition* → PRIMARY LIBRARY → EXTENSION → EXTENDED LIBRARY →
   PROGRESSIVE ALIGNMENT.
2. **Primary-library weight (p.207):** "Each constraint receives a weight equal to percent
   identity within the pair-wise alignment it comes from." In the GARFIELD worked example the
   SeqA/SeqB pairwise alignment has **Prim Weight = 88**.
3. **Library combination (p.207):** "If any pair is duplicated between the two libraries, it is
   merged into a single entry that has a weight equal to the **sum** of the two weights." Pairs not
   present "will be considered to have a weight of zero."
4. **Extension — triplet/consistency (pp.208–209):** "a triplet approach is used … checking the
   alignment of the two residues with residues from the remaining sequences." Worked: A(G),B(G)
   direct weight 88; through C, W1 = W(A(G),C(G)) = 77, W2 = W(C(G),B(G)) = 100, "we associate that
   alignment with a weight equal to the **minimum** of W1 and W2 … set to 77. In the extended
   library, this new value is **added** to the previous one to give a total weight of **165**
   (i.e. 77 + 88)."
5. **Extended weight = sum of triplet supports (p.209):** "the weight associated with a pair of
   residues will be the **sum of all the weights** gathered through the examination of all the
   triplets involving that pair. The more intermediate sequences supporting the alignment of that
   pair, the higher its weight." Uninformative triplets (e.g. through D) contribute 0.
6. **Progressive (pp.209–210):** guide tree from pairwise distances (NJ in the paper); "The closest
   two sequences … are aligned first using normal dynamic programming. This alignment uses the
   weights in the extended library." Once-a-gap-always-a-gap: "any gaps … cannot be shifted later."
   Group-vs-group: "the average library scores in each column of existing alignment are taken."
7. **No gap penalties in the progressive DP (p.210):** "a dynamic-programming algorithm (Gotoh,
   1982) with **gap-opening penalties and gap-extension penalties set to zero**."
8. **Objective class (p.209):** the substitution scores are the extended-library consistency
   weights, optimising a consistency-based objective — distinct from a fixed-matrix SP score.

### Source T2 — T-Coffee Technical Documentation (canonical project docs, rank 3)

**URL:** https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html
**Retrieved by:** WebFetch.   **Accessed:** 2026-06-23   **Authority rank:** 3

**Key Extracted Points:**

1. Confirms the **min() rule** for the two legs of a triplet through an intermediate sequence.
2. Triplet extension modes propagate constraints; "the job of T-Coffee is to satisfy as many
   constraints as possible."

## Documented Corner Cases (T-Coffee)

1. **Never-aligned residue pairs → weight 0** (T1 p.207, p.209): such pairs are not in the library;
   the DP scores them 0.
2. **Uninformative triplet → contributes 0** (T1 p.209): extension never *lowers* a weight.
3. **Identical sequences:** all pairwise alignments are gapless 100% identity → trivial exact MSA.
4. **k = 2:** no intermediate sequences → extended library = primary library → consistency alignment
   = the single pairwise global alignment.

## Test Datasets (T-Coffee)

### Dataset: GARFIELD worked example (T1, Figure 2)

| Quantity | Value | Source |
|----------|-------|--------|
| Primary weight A–B | 88 | Figure 2(b) |
| W(A(G),C(G)) | 77 | p.209 |
| W(C(G),B(G)) | 100 | p.209 |
| Triplet through C for A(G),B(G) | min(77,100)=77 | p.209 |
| Triplet through D | 0 (uninformative) | p.209 |
| **Extended weight A(G),B(G)** | **165 = 88 + 77** | p.209 |

This proves the MUST property "extended weight of a consistency-supported pair > its primary weight"
with exact integers. (Re-expressed over a DNA alphabet for the DNA-only API; the integer relation
extended = primary + Σ min-triplets is alphabet-independent.)

### Dataset: 3-sequence DNA consistency case (derived from T1 formulas)

**Source:** derived from primary = %identity and extension = Σ min-triplet (T1 pp.207–209);
independently recomputed in test code with a stand-alone library builder, not from implementation
output.

## Assumptions (T-Coffee)

1. **ASSUMPTION: Primary weight = round(100 × identical-aligned-columns / alignment-length)** of the
   pairwise *global* alignment (integer percent identity), matching the integer GARFIELD numbers
   (88, 77, 100). Local-library entries are added from the existing Smith–Waterman local alignment
   of each pair, weighted by the local segment's percent identity; combination is signal addition
   (sum of weights for duplicated residue pairs) exactly per T1 p.207. Local segments that are a
   subset of the global alignment simply reinforce existing pairs (sum), which is the intended
   "stacking" behaviour.
2. **ASSUMPTION: Guide tree = the repository's existing UPGMA tree** (`BuildProgressiveGuideTree`).
   T1 uses NJ; UPGMA vs NJ changes only merge order, not the library or the consistency objective.
   Reused for determinism and parity with `MultipleAlignProgressive`.
3. **ASSUMPTION: Progressive DP gap score = 0** (T1 p.210): a gap column contributes 0 library score;
   group-vs-group columns use the average extended-library score over all cross-group residue pairs
   (T1 p.210). No fixed substitution matrix is consulted in the progressive phase.

## Recommendations for Test Coverage (T-Coffee)

1. **MUST:** extended-library weight of a consistency-supported pair > its primary weight and equals
   primary + Σ min-triplets (GARFIELD 88 → 165 relation) — Evidence: T1 p.209.
2. **MUST:** a residue pair supported by an intermediate has strictly greater extended weight than an
   inconsistent (unsupported) pair — Evidence: T1 p.209.
3. **MUST:** identical inputs → trivial exact alignment (rows == inputs, no gaps) — Evidence: corner case.
4. **MUST:** valid MSA — all rows equal length; degap recovers each input exactly; no all-gap column — Evidence: MSA definition.
5. **MUST:** consistency alignment improves (or equals) the T-Coffee consistency objective vs the plain progressive seed on a case engineered to mislead one pairwise alignment — Evidence: T1 central claim.
6. **SHOULD:** k = 2 → consistency alignment = single pairwise global alignment — Evidence: extended = primary, no triplets.
7. **SHOULD:** null / empty / single-sequence parity with sibling MSA methods; sibling aligners unchanged — Evidence: API contract / additivity.
8. **COULD:** determinism — same input twice → byte-identical output — Evidence: no RNG.

## References (T-Coffee addendum)

10. Notredame C, Higgins DG, Heringa J (2000). T-Coffee: A novel method for fast and accurate multiple sequence alignment. *J Mol Biol* 302(1):205–217. https://doi.org/10.1006/jmbi.2000.4042 (full text https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf ; PubMed https://pubmed.ncbi.nlm.nih.gov/10964570/)
11. T-Coffee Technical Documentation. https://tcoffee.readthedocs.io/en/latest/tcoffee_technical_documentation.html (accessed 2026-06-23)
12. Gotoh O (1982). An improved algorithm for matching biological sequences. *J Mol Biol* 162(3):705–708. https://doi.org/10.1016/0022-2836(82)90398-9


---

# Addendum: B09 re-validation 2026-11-04 (campaign 2026-09)

Sources actually opened in this review (neutral `Seqeron-review/1.0` User-Agent, no personal data sent):

| Source | How | What it confirms |
|---|---|---|
| T-Coffee paper PDF, https://tcoffee.org/Publications/Ps_pdf/tcoffee.pdf (pp.207-210, `pdftotext`) | curl | Constraint weight = "percent identity within the pair-wise alignment it comes from"; signal addition; "Extension will be carried out on each pair of residues of A and B"; triplet weight = minimum of W1, W2 (88 + 77 = 165); group columns use average library scores; DP with gap penalties set to zero; NJ guide tree |
| T-Coffee sources, `lib/dp_lib/evaluate.c` `residue_pair_extended_list` (raw.githubusercontent.com/cbcrg/tcoffee) | curl | Extension is evaluated for *any* residue pair (s1,r1),(s2,r2) the DP asks for: direct weight + Σ MIN(w(s1.r1,t), w(s2.r2,t)) over common partners t |
| T-Coffee sources, `lib/util_lib/aln_convertion_util.c` `get_seq_sim` | curl | Default `-weight=sim` = `(sim*MAXID)/pos0`: identical columns × 100 / aligned non-gap columns, integer division |
| T-Coffee technical documentation (readthedocs) | WebFetch | `-weight=sim` default; `-extend_mode=very_fast_triplet` default |
| ClustalW 1.83 `pairalign.c` (raw.githubusercontent.com/epruesse/ARB/master/GDE/CLUSTALW) | curl | Pairwise distance `tmat = (100 − mm_score)/100`, `mm_score = 100·count/MIN(len1,len2)` (identities over the shorter sequence length) |
| Abouelhoda & Ohlebusch 2005, J Discrete Algorithms 3:321-341 (uni-ulm.de PDF) | curl, `pdftotext` | Chaining problem = maximum-weight chain of colinear non-overlapping fragments; MUMmer cited as a chaining client |
| Gusfield 1993, Bull Math Biol 55:141-154 (WebSearch snippets: Springer/MaRDI) | WebSearch | Center-star method, bound 2 − 2/k |
| Higgins & Sharp 1988, Gene 73:237 (PubMed snippet) | WebSearch | UPGMA-based guide tree, simplified Feng-Doolittle |

## Oracle numbers

- Independent Python implementations (own NW/SW with the documented tie order, Gusfield center-star, ClustalW distance + UPGMA with exact `Fraction` arithmetic, profile-profile NW, MUSCLE edge refinement, T-Coffee primary library + on-demand triplet extension + zero-gap progressive DP) compared with 300 random sets (k = 2..6, three scoring schemes): star 0/300, progressive 0/300, iterative 0/300, consistency 0/300 row mismatches. Before the fixes: consistency 170/300 rows differed from the T-Coffee reference semantics; k = 2 star MSAs were suboptimal in 4 of 2000 related pairs (anchored pairwise stage; e.g. 9 vs Biopython 24).
- Biopython 1.88 `PairwiseAligner(mode='global', match_score=1, mismatch_score=-1, open_gap_score=-1, extend_gap_score=-1)`: GGAAGCCAATAATCGTGCACGAAGGTGACTGCAGGAGAGTACGAAGGT / AGACAGCCAAAATCTTCAGAGTTCGACTGGAGGTGACTACGAAGAT = 24; TAGACGCATCAAACAGCCTAGGCCGTACTGCGATTCAGCTGTCAACAGGTAAAACAGA / TAGACGCCTAGTCAACAGCCTACGGCCCTATCAGCTGTATACATGTTAGTCCACCGTAACAGT = 27; CACTTGGGGACAGTTAGCTCCGCCCCCGTTAGTTTGATGCGGGTCC / ACGGCGCGTAGGAAGTATAGTAATTCTTAGATGGGGACACGATACTA = −3.
- Star rows/SP (python, hand-checked): ACACG/AGAGC/AATGGC → `-A-CACG-`, `-A-GA-GC`, `AATG--GC`, SP −3; TTTGA/AGGTCA/AGTGCGA/ACTGGCA → SP −5; CAT/TAT/CGAT/CTAT/CTA → SP 10.
- Chaining: A(0,0,10)+A(6,20,10) overlap 4 → chain (0,0,10),(10,24,6) = 16 bases.
- T-Coffee: ACGTTACGT/ACGACGT/ACGTACGA primary(S0.0,S1.0) = 200 (global 7/7 + local 7/7), extended 385 = 200 + min(87+100, 85+100) = 200 + 185; CTT/ACT/CAA pair (S1.1,S2.0): primary 0, extended 133 = min(200, 133) via S0.
- Iterative: TCA/TCT/GCAC/TAAC progressive `T--CA T--CT GCAC- TAAC-` SP −7, refined `TCA- TCT- GCAC TAAC` SP −3.

- T-Coffee binary (Debian `t-coffee` 13.45, apt): same explicit primary library + same UPGMA tree + `-extend_mode slow_triplet`. Rows equal in 15/120 (k = 3) and 13/150 (k = 4-5) random sets; the rest differ because the binary rescales each extended pair score by its residue's partner-weight total (`evaluate.c` L5245-5250), which is not in the paper. Under the paper's objective (Σ extended weights over aligned residue pairs) the binary's alignment was never higher than ours for k = 3 (0/120) and higher in 7/150 sets for k = 4-5. `-lib_only -out_lib` does not extend, so the extended weights were not read out of the binary; the extension rule is confirmed by the paper text (p.208-209) and `residue_pair_extended_list`. **Superseded by A1-13 (F58, below): the difference does not come from `-extend_mode` (no effect with the default DP) but from the program's relaxation + linked pair-wise DP.**

## Assumptions (B09)

1. The linear-gap model applies to every MSA method (documented deviation; `GlobalAlignAffine` for pairs).
2. UPGMA guide tree for all progressive-family methods (CLUSTAL 1988) by default; NJ available as `MsaOptions.GuideTree` (A1-12, addendum at the end).
3. T-Coffee local library = best single SW alignment per pair (Lalign supplies up to N non-intersecting alignments in T-Coffee).


---

# Addendum: ClustalW pairwise stage as guide-tree distance (audit A1-15, B09 F51, 2026-11-18)

## Sources actually opened
- ClustalW 2.1 (Debian `clustalw 2.1+lgpl-7`) source via `sources.debian.org/data/main/c/clustalw/2.1+lgpl-7/src/`: `pairwise/FullPairwiseAlign.cpp` (`pairwiseAlign`, `forwardPass`, `reversePass`, `diff`, `tracePath`), `substitutionMatrix/SubMatrix.cpp` (`getPairwiseMatrix`, `getMatrix`, `setUpCrossReferences`), `substitutionMatrix/matrices.h` (`swgapdnamt`), `general/UserParameters.cpp` (`DNAPWGapOpen = 15.0`, `DNAPWGapExtend = 6.66`, `transitionWeight = 0.5`, `aminoAcidCodes`), `Clustal.cpp` (full vs fast pairwise selection, `-PIM`/`countid` path), `pairwise/FastPairwiseAlign.cpp` (opened, not ported: A1-15b). Neutral User-Agent `Seqeron-review/1.0`.
- Installed binary `clustalw` 2.1 (`-TYPE=DNA -ALIGN`): prints `Sequences (i:j) Aligned. Score: N`, N = (int)(100·count/MIN(len1,len2)) in float. `-TREE -PIM` writes a *column-wise* identity matrix of the input rows (e.g. 54 for the pair whose pairwise Score is 91), so it cannot serve as the oracle for the pairwise stage.

## Findings from the source
1. The pairwise stage is an affine **local** alignment (Smith-Waterman-Gotoh, end = first row-major maximum; reverse pass for the start; Myers-Miller `diff` on the region), not a global one.
2. DNA pairwise matrix = IUB ×100: identical +1000, ACGT mismatch −900 (including transitions: the 2.1 transition-weight statements write to code pairs (0,4), (2,11), (2,12) = A-E, C-M, C-N). Gap open `(int)(2·15·100) = 3000`, extend `(int)(6.66f·100) = 666`; k-gap costs 3000 + 666·k.
3. identity count uses residue codes (case-insensitive, N–N identical); distance = (100 − 100·count/MIN(len1, len2))/100; either sequence empty → 1.0.

## Oracle numbers
3 600 pairs equal (printed Score, identities, DP score); 140 MSA sets × {progressive, iterative} equal the independent Python; 65/140 progressive results differ from the linear-distance guide tree.

## Assumptions / not covered
- (State of the first F55 pass; superseded) Full-DP ("slow") mode only. The `-QUICKTREE` fast k-tuple mode and the protein pairwise stages were ported in audit A1-15b (`MsaDistanceMode.ClustalWQuick`, `ClustalWProteinPairwiseOptions`; Multiple_Sequence_Alignment.md §4.2 / §4.6).
- The guide tree itself stays UPGMA (ClustalW uses NJ: audit A1-12), so a `clustalw -USETREE`-driven comparison of final rows is not claimed.


---

# Addendum: affine-gap multiple alignment (audit A1-10, B09 F52, 2026-11-18)

## Sources actually opened / used
- Altschul SF (1989) "Gap costs for multiple sequence alignment", J Theor Biol 138:297-309 (WebSearch record: "if SP substitution scores are used, gap scores should also be defined as the sum of gap scores for the implied pairwise alignments" = natural gap costs; the quasi-natural definition could not be retrieved, so it is not implemented or claimed).
- ClustalW 2.1 / MUSCLE (Edgar 2004) / Gotoh 1982 / Flouri et al. 2015 as already cited by F1; ClustalW `FullPairwiseAlign.cpp` for the gap convention (open + k·extend).
- Installed tools used for an informational comparison only: `clustalw` 2.1, `muscle` 5.2, `mafft` 7.505, `kalign` 3.3.5, `clustalo` 1.2.4, `t_coffee` 13.41.

## Decisions
1. Objective = natural affine sum of pairs (exactly computable). The exact optimisation of that objective for k ≥ 3 is not a three-state DP (the open charge of a pair depends on whether the previous non-gap-gap column of THAT pair was a gap of the same row), so the profile methods use the standard occupancy-scaled profile approximation and are scored with the exact objective; the exact optimum for three sequences is the 7-move DP over (i, j, l, run state per pair) used as the test oracle (audit A1-11 will expose it as an API).
2. k = 2 is exact: the three-state profile DP reduces to the `GlobalAlignAffine` recurrence (same border, same tie order), rows equal.

## Oracle numbers
- 453 random MSAs: independent re-score = C# `SumOfPairsScore` (linear + affine); affine SP ≤ Σ Biopython pairwise optima in all rows.
- 260 tiny sets (k = 2/3, length ≤ 7 (k = 3: ≤ 6)): exact optimum by DP, DP = exhaustive enumeration whenever Σ length ≤ 9; k = 2 optimum = Biopython; C# methods never exceed it; k = 2: star/progressive/iterative = optimum; 227 three-sequence sets: optimum reached by star 170, progressive 189, iterative 189 (iterative never lowers the exact objective).
- 150 sets (k 3-6, ≤ 24 nt, 7 scoring schemes): C# affine star, progressive, iterative and consistency rows = independent Python rows (0 differences); progressive rows differ from the linear model on > 20 sets.
- Informational comparison under the same objective (natural affine SP with match 2, mismatch −3, open −5, extend −2) on 60 sets (k 4-7, 30-80 nt, 8-25 % divergence): mean SP of Seqeron affine iterative −94.7, consistency −84.4, progressive −102.6, star −144.4, anchored −190.1, Seqeron linear progressive −278.3 / linear star −300.1; external tools (each optimises its own objective and parameters, so this favours the Seqeron methods): muscle −240.9, mafft −299.1, clustalw −372.2 (−310.0 with open 12.5 / extend 5), kalign −640.1, clustalo −636.0, t_coffee −148.1. An external tool scores higher than Seqeron's affine iterative result on 11 sets (t_coffee), 6 (kalign), 5 (mafft), 2 (muscle), 1 (clustalw) and 0 (clustalo) of the 60; the Seqeron consistency method (zero-gap-penalty library DP, affine library) is higher than iterative on 25 and lower on 29.

## Assumptions / not covered
- The profile DP is not ClustalW's `prfalign.c` stage (sequence weights, residue-specific / hydrophilic penalties, gap-separation distance) nor MUSCLE's profile score: audit A1-10b.
- Consistency keeps T-Coffee's default gap penalty 0 in its progressive DP; the gap model affects its pairwise library and the reported score only.


---

# Addendum: neighbour-joining guide tree, canonical tree builders, O(a log² a) anchor chaining (audit A1-12 / A1-16, B09 F53 / F54, 2026-11-18)

## Sources actually opened
- ClustalW 2.1 (Debian `clustalw 2.1+lgpl-7`) source via `sources.debian.org/data/main/c/clustalw/2.1+lgpl-7/src/`: `tree/NJTree.cpp` (`generateTree`: loop nest jj-outer / ii-inner with strict `<`, update `(d_ik + d_jk)/2` with `av[]` height correction, |branch| < 1e-4 → 0), `tree/ClusterTreeOutput.cpp` (`printPhylipTree`, `twoWaySplit`: children listed by lowest leaf index, `fixed setprecision(5)`), `tree/Tree.cpp` (`readTree`, `createTree`, `insertNode`, `reRoot`, `calcMean`, `insertRoot`, `createSets`, `groupSeqs`, `calcSeqWeights`), `tree/TreeInterface.cpp`, `tree/UnRootedClusterTree.cpp`, `pairwise/FullPairwiseAlign.cpp` (distance in single precision), `general/clustalw.h` (`NODE`/`LEAF`/`LEFT`/`RIGHT`), `tree/UPGMA/*` (ClustalW's own `-CLUSTERING=UPGMA`, read, not used as oracle).
- Installed `clustalw` 2.1: `-TYPE=DNA -ALIGN -MAXDIV=0 -NEWTREE=x.dnd` writes the unrooted NJ guide tree and prints "Group n: Sequences: s" for every merge of the rooted tree (`-MAXDIV=0` switches off delayed divergent sequences, which would print "Delayed" instead).
- Biopython 1.88 `Bio.Phylo.TreeConstruction` (`nj`; its `upgma` source read: `dm[min_j, k] = (dm[min_i, k] + dm[min_j, k]) / 2` — no cluster sizes, i.e. WPGMA), scipy `cluster.hierarchy.linkage(method='average')` (UPGMA).
- Saitou & Nei 1987 (Mol Biol Evol 4:406), Sokal & Michener 1958, Thompson, Higgins & Gibson 1994 (NAR 22:4673), Abouelhoda & Ohlebusch 2005 (J Discrete Algorithms 3:321; opened in F5) for the chaining recurrence.

## Findings
1. `Seqeron.Genomics.Phylogenetics` referenced `Seqeron.Genomics.Alignment` without using a single Alignment type (`using` directive only); removing the project reference and the global `using` leaves `dotnet build Seqeron.sln` at 0 errors and the Architecture tests (`Modules_ShouldBeFreeOfCycles`) green. `Alignment` now references `Phylogenetics`. The former private UPGMA of `SequenceAligner` (justified by that reference) is deleted; the guide tree is the canonical `PhylogeneticAnalyzer` tree.
2. ClustalW's neighbour-joining update `(d_ik + d_jk)/2` (instead of Saitou-Nei's `(d_ik + d_jk − d_ij)/2`) shifts every Q by a constant per step, so topology and branch lengths equal standard NJ; only exact ties depend on the loop order (jj-outer), which is reproduced by `PairScanOrder.SecondIndexMajor`.
3. The progressive order is not the NJ tree itself: ClustalW roots it with `reRoot` (single precision, lengths re-read from the 5-decimal `.dnd`) and lists the merges left subtree first.
4. Anchor chaining: with δ = query − ref the trim of a predecessor is a function of one end coordinate on each side of δᵢ, so the all-pairs scan is two orthogonal range-maximum queries per side.

## Oracle numbers
420 sets: unrooted tree = `.dnd` (0 mismatches), merge-size sequence = "Group n" lines (0 mismatches; first-index tie order: 21 mismatching sets, tolerance 0 instead of 1e-9 changes nothing on this panel); 400 matrices: NJ = Biopython (≤ 1e-6), UPGMA = scipy average linkage (≤ 1e-6); 900 matrices: UPGMA guide tree = the removed private implementation; chaining: 1 500 random + 400 MEM-like + 12 adversarial + 2 × 3 000-anchor sets equal the O(a²) oracle; 10⁵ anchors in < 1 s.

## Assumptions / not covered
- (State of the F55 pass; superseded) Final MSA rows were not compared with `clustalw -USETREE` then: ClustalW's profile stage (sequence weights from the tree, residue-specific and hydrophilic gap penalties, gap-separation distance, end gaps, delayed divergent sequences) was ported in audit A1-10b (B09 F56; `MultipleAlignClustalW`, rows equal on 1 000 DNA sets and later the protein stage). ClustalW's own UPGMA (`-CLUSTERING=UPGMA`) is a separate implementation with its own tie rules: audit A2-14 (see B09.md).


---

# Addendum: ClustalW fast pairwise stage and multiple-alignment stage (audit A1-15b / A1-10b, B09 F55 / F56, 2026-11-18)

## Sources actually opened
- ClustalW 2.1 (Debian `clustalw 2.1+lgpl-7`) source via `sources.debian.org/data/main/c/clustalw/2.1+lgpl-7/src/`: `pairwise/FastPairwiseAlign.cpp` (`pairwiseAlign`, `pairAlign`, `makeNPtrs`, `putFrag`, `fragRelPos`, `desQuickSort`), `alignment/Alignment.cpp` (`calculateMaxLengths`: `maxAlignmentLength = 2 · longest`), `general/UserParameters.cpp` (DNA defaults: `DNAKtup` 2, `DNAWindow` 4, `DNASignif` 4, `DNAWindowGap` 5, `DNAGapOpen` 15, `DNAGapExtend` 6.66, `divergenceCutoff` 30, `INT_SCALE_FACTOR`), `interface/CommandLineParser.cpp` (`-KTUPLE` sets the window gap to ktup + 4; `-WINDOW/-TOPDIAGS/-PAIRGAP` are accepted only above the k-tuple; `-ENDGAPS` sets the default false), `Clustal.cpp` (`align`), `tree/Tree.cpp` (`calcSeqWeights`, `calcWeight`, `calcSimilarities`, `createSets`, `groupSeqs`), `tree/NJTree.cpp` (`generateTree`), `tree/TreeInterface.cpp`, `tree/UnRootedClusterTree.cpp` (two-sequence tree file), `multipleAlign/MSA.cpp` (`multiSeqAlign`), `multipleAlign/MyersMillerProfileAlign.cpp` (`profileAlign`, `progDiff`, penalty accessors, `addGGaps`), `multipleAlign/ProfileBase.cpp` (`calcGapCoeff`), `ProfileStandard.cpp`, `ProfileWithSub.cpp`, `substitutionMatrix/SubMatrix.cpp` (`getProfileAlignMatrix`, `getMatrix`, `setUpCrossReferences`), Thompson, Higgins & Gibson 1994, Wilbur & Lipman 1983, Myers & Miller 1988.
- Installed binary `clustalw` 2.1: `-TYPE=DNA -QUICKTREE -ALIGN` prints "Sequences (i:j) Aligned. Score: %lg" / "Not Aligned"; `-ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` writes the final rows.

## Findings from the source
1. The fast stage is a Wilbur-Lipman diagonal / fragment-chain score, not an identity; the unstable quicksort and the partial-alignment exit (fragment list ≥ 2·maxAlnLength, maxAlnLength from the longest sequence of the whole input) are part of the result.
2. `-KTUPLE=k` also sets the pair gap to k + 4: `-KTUPLE=2` alone is pair gap 6, the default without `-KTUPLE` is 5.
3. ClustalW's multiple-alignment stage replaces the pairwise identities by tree-derived similarities and takes the progressive order, sequence weights and delayed sequences from the tree; the profile matrix is IUB ×100 shifted by +900, transitions (A-G, C-T, C-U) = ⌊TRANSWEIGHT·1900⌋; `-GAPEXT` is a float (100·6.66f = 665.99998 → 665); the gap-separation / hydrophilic / residue-specific penalties are guarded by `!_DNAFlag` (no effect for DNA); `-ENDGAPS` is a no-op on the command line.
4. `-TRANSWEIGHT` also reaches the pairwise matrix: it writes ⌊w·1000⌋ to the CODE pairs C-M and C-N (ClustalW 2.1 uses residue codes there), which only matters for IUB-ambiguity input.
5. NJ ties: ClustalW's neighbour joining (height-corrected matrix, running sums) decides exact ties by rounding of its own arithmetic.

## Oracle numbers
- Fast stage: 8 257 pairs of 605 sets, all `Score:` equal (see TestSpec CQ1); development panel of 10 733 pairs equal.
- Multiple alignment: 1 000 committed sets, every row equal; development panels of 700 and 2 500 sets equal after the fixes of F56 (first run: 148 of 150 sets equal, then 698 of 700, then 2 496 of 2 500).

## Assumptions / not covered
- (F56 state, superseded: the `CLUSTALW` DNA matrix and user matrices are ported in F68 below, the protein stage in F63, user protein matrices in F64, the end-gap options in F67.) Not ported then: secondary-structure masks, `-ITERATION`, `-CLUSTERING=UPGMA` (the last and the `-TYPE` guess since audit A2-14, F92: below); `-ENDGAPS` on the command line has no effect (F67).


## T-Coffee 13.41 progressive stage (audit A1-13, B09 F58)
- Program: Debian `t-coffee` 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1 (`t_coffee`, "Version_13.41.0.28bdc39"); source of exactly that version read at sources.debian.org (files in B09.md F58).
- Plain default: `t_coffee seqs.fa` on DNA runs `Mproba_pair` (`t_coffee.c` L4371-4377, "new default"); `-mode=regular` gives the same rows; the paper's `lalign_id_pair` + `slow_pair` libraries must be requested (`-method slow_pair,lalign_id_pair`).
- Progressive stage (what `MultipleAlignTCoffee` ports; weights of a library file are used as given, the program's own methods write 10 x percent identity): (1) `-relax_lib 1`: for every directed edge of residues r1 < len: `trunc((w_back + sum over common neighbours min(w1, w2)) / min(deg1, deg2))`, double accumulation; last residues keep their weight; edges with weight <= 10 dropped. (2) `dp_mode linked_pair_wise`: `cl2pair_list_ext` (float: min(w1, w2) / 1000 summed per column pair over the library neighbours outside profile 1, divided by their count, x 1000, kept when > 1 % of the column's mass, not in the last row/column, `%f` round trip), `cl2diag_cap`, `cl2list_borders`, `list2nodup_list`, `list2linked_pair_wise` (three states, gap open = gap extend = 0 (`MANUAL PENALTIES: gapopen=0 gapext=0`), ties M >= I,J then I >= J, smaller of (sequences x columns) first); the tree is processed recursively, left subtree is profile 1.
- `-extend_mode` (`slow_triplet`, `triplet`, `fast_triplet`) does not change the rows with this DP: 0 of 60 sets differ.
- Oracle: `scripts/alignment/generate_tcoffee_oracle.py` -> `tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/tcoffee_oracle.tsv` (400 sets x relax 1/0 = 800 rows; libraries "rand" = NW + SW of this repository, weight 10 x pid, and "binary" = the program's `-method slow_pair,lalign_id_pair -lib_only -out_lib`; random binary trees; DNA, 3-6 sequences). **800 / 800 row sets equal** (`SequenceAligner_TCoffee_Tests`). Python port development panels: 0 / 100 mismatches at relax 0 for k = 3, 4, 5, 6 (first try), relax 1: 0 / 100 for k = 4, 5, 6 after the filter fix (strict `> 10`).
- Not reproduced (A1-13b): (a) primary library - the program's `slow_pair` is Myers-Miller (`util_dp_mm_nw.c`) global alignment with `dna_idmat` (+10/-1, scaled by SCORE_K = 10), gap penalties `gop = get_avg_matrix_mm(M, AA_ALPHABET)*10` (an average over amino-acid letter indices of a DNA matrix), `gep = -1`, TG_MODE 1 (terminal gaps pay `gep` only), `lalign_id_pair` is SIM top-10 (`-lalign_n_top 10`) local alignments, weights 10 x percent identity, and the plain default is `Mproba_pair` (pair-HMM); the gap model fitted on 60 random pairs (M = 10, X = -1, gep = -1, terminal gep only) makes the binary's alignment optimal for 54/60 pairs at gop = -11 (43/60 at -10, 39/60 at -12), so it is not yet an exact model of the Myers-Miller midpoint rule; this repository's NW (+1/-1/-2) + single SW library equals the program's `slow_pair,lalign_id_pair` library in 0/150 sets (pair sets 0/150); (b) the program's default guide tree (`-tree_mode nj` on k-tuple distances, trifurcating root converted to a binary order) is not ported. End to end on 150 related-sequence sets, `MsaConsistencyEngine.TCoffeeLinked` (own library, NJ tree) equals `t_coffee -method slow_pair,lalign_id_pair` in 21 sets (the paper-formula engine in 19) and a plain `t_coffee` in 20 (paper 21): that is the size of the remaining library/tree gap, split off as A1-13b.


## Audit r1 addendum (A1-11, B09 F60): exact sum-of-pairs multiple alignment

- **Objective** = `SumOfPairsScore` of the heuristic methods (linear column sum; Altschul 1989 natural affine). The optimum of a set is a best lattice path (2^k - 1 neighbours, Altschul & Lipman 1989); with natural affine gaps the per-pair run state is a function of the order in which the sequences last moved (weak orders: 13 / 75 / 541 states for k = 3 / 4 / 5).
- **Carrillo-Lipman bound**: Shamir's lecture notes (`cs.tau.ac.il/~rshamir/algmb/98/scribe/html/lec03/node11.html`) opened - pairwise budget `B(u,v) = c' - sum_{(i,j) != (u,v)} D(x_i, x_j)` and pruning of cells whose through-cell pairwise bound exceeds it; the implemented test (the sum over all pairs of the through-cell optimum cannot reach the lower bound `L`) implies it; plus the A*-style state test with the pairwise suffix tables. The original papers (Carrillo & Lipman 1988; Lipman, Altschul & Kececioglu 1989) were found only as abstracts / records.
- **Oracle** `scripts/alignment/generate_msa_exact_oracle.py` -> `msa_exact_oracle.tsv` (804 sets): optimum by two independent Python programs (lattice DP; DP over positions x run state of every pair) which agree with each other, with exhaustive enumeration of every alignment (tiny sets) and (k = 2) with Biopython 1.88; C# equal on 804 / 804 with and without the bound.
- **Heuristic optimality measured** (k = 3..4, total length <= 30, 309 sets per objective): linear - star 218, progressive 276, iterative 276, consistency 186, anchored 187, best of five 285; affine - 212 / 227 / 228 / 164 / 178, best of five 265. Never above the optimum.
- **Speed** (Debug build): k = 3, 300 nt (28 M lattice cells, 10 % divergence) 0.18 s / 1 300 cells visited; k = 4, 100 nt (110 M cells) 30 ms; k = 5, 50 nt (331 M cells) 25-66 ms; k = 6, 40 nt affine 0.6 s; k = 4, 60 nt at 30 % divergence affine 0.8 s / 25 660 cells.


## Protein / substitution-matrix MSA (audit A1-14, B09 F61)
- Model: the pairwise matrix model of F41 (EMBOSS needle / water numbers, Biopython `PairwiseAligner(substitution_matrix=..., open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend)`) lifted to the multiple aligners: sum of pairs = Σ matrix(a, b) over residue pairs of a column, residue / gap = `GapExtend` (linear) or the natural affine run cost `GapOpen + k·GapExtend` of every pair projection (Altschul 1989); gap-gap columns neutral.
- Oracles: `scripts/alignment/generate_msa_matrix_oracle.py` (independent Python, exact Fractions, scaled integers; nothing of the C# code reused). `msa_matrix_oracle.tsv` 238 sets x 5 methods (k = 2..6; BLOSUM62 / PAM250 / EDNAFULL; 15 gap schemes incl. 9.5/0.5, 4.9/0.1, 0.25/0.75; linear and affine): **all rows equal**; k = 2 asserted equal to the Biopython 1.88 optimum inside the generator. `msa_matrix_exact_oracle.tsv` 308 sets: lattice DP = Biopython (k = 2) = exhaustive enumeration (tiny sets); C# equal with and without the Carrillo-Lipman bound. `msa_matrix_sp_oracle.tsv` 420 random MSAs: SP linear / affine and the pairwise-optimum bound. `tcoffee_protein_oracle.tsv`: `t_coffee` 13.41, 160 random protein sets x relax 1 / 0 = 320 rows of the explicit-library progressive stage: **320 / 320 equal** (`generate_tcoffee_oracle.py OUT 160 17 ACDEFGHIKLMNPQRSTVWY`; the DNA file is regenerated byte-identically by the default alphabet).
- Informational (`scripts/alignment/compare_msa_protein_tools.py`, 10 random families of 4-6 sequences, 45-70 residues, 10-40 % divergence; SP under BLOSUM62 / gapopen 10 / gapextend 0.5 natural affine): mean SP mafft 1072.5, seqeron-consistency 1067.4, muscle 1054.1, kalign 1024.4, seqeron-star 1021.5, seqeron-iterative 1013.3, t_coffee 1009.9, seqeron-progressive 995.8, clustalw 910.7, clustalo 753.4; upper bound (sum of the pairwise optima) 1257.2. The programs optimise their own objectives (profile scores, consistency, HMM); lower SP is not an error and nothing is asserted.
- Anchored aligner: the oracle reproduces the suffix-tree anchor semantics (first peak of the matching statistic per run, one anchor per run; sets with a repeated anchor text are skipped because the tree reports an arbitrary occurrence) + the colinear chain with trimming + Needleman-Wunsch / Gotoh between anchors; 276 of 794 center-vs-sequence pairs have at least one anchor.


## ClustalW PROTEIN pairwise stages (audit A1-15b part b, B09 F62)
- Sources actually opened (sources.debian.org, clustalw 2.1+lgpl-7, `src/`): `pairwise/FullPairwiseAlign.cpp` (all of it: the protein gap-penalty branch `_matAvgScore <= 0 ? 2 * (int)((PWGAPOPEN + ln MIN) * intScale) : (int)(2 * avg * (PWGAPOPEN + ln MIN) * gapOpenScale)`, `_gapExtend = (int)(PWGAPEXT * intScale)`), `pairwise/FastPairwiseAlign.cpp` (`makePPtrs`, protein `pot[]` = (maxAA + 1)^(j - 1)), `substitutionMatrix/SubMatrix.cpp` (`getMatrix`, `getPairwiseMatrix`: matrix choice, intScale 100 / GONNET 10, `matrixAvgScore`), `substitutionMatrix/matrices.h` (tables `gon250mt`, `blosum30mt`, `pam350mt`, `idmat`), `general/UserParameters.cpp` (protein defaults: PWGAPOPEN 10, PWGAPEXT 0.1, KTUP 1, WINDOW 5, SIGNIF 5, WINDOWGAP 3; `aminoAcidCodes`), `interface/CommandLineParser.cpp` (`-PWMATRIX`, `-KTUPLE <= 2` sets the pair gap to k + 3, window / topdiags / pairgap accepted only above k).
- Oracle: `scripts/alignment/generate_clustalw_protein_pairwise_oracle.py OUT matrices.h EBLOSUM30` -> `clustalw_protein_pairwise_oracle.tsv` (4 204 pairs / 366 sets: matrices GONNET / BLOSUM / PAM / ID, gap variants (10, 0.1) (5, 0.5) (20, 0.05) (0, 0.1) (10, 1.0) (7.5, 0.25) (15, 0.2) (2, 2), 20-letter and 23-letter (B, Z, X) alphabets, two-letter and tandem-repeat sequences, 5-400 residues) - every `Score:` of the installed binary is equal to the line-by-line Python port (CW class of `generate_clustalw_pairwise_oracle.py` + the matrices.h tables + the formulas above in single precision) before the row is written; `clustalw_protein_quick_oracle.tsv` (3 774 pairs / 335 sets, `-QUICKTREE`, 357 "Not Aligned"); both equal the C# ports on every row.
- Findings of the oracle: `-PWMATRIX` really changes the scores (the three-globin hand run happens to print 40 / 38 / 31 for all four matrices); the gap open of the identity matrix uses the `avg <= 0` branch (matrixAvgScore 0); matrixAvgScore = 82 (GONNET), 113 (BLOSUM30), 142 (PAM350), 0 (ID); `J` and `*` are silently dropped by the program (rejected here), `O` / `U` are residues without table row (score 0, accepted).
- Matrix tables: BLOSUM30 (`blosum30mt`) = EMBOSS `EBLOSUM30` on all 529 cells; GONNET250 (`gon250mt`) = round(10 x Biopython 1.88 `GONNET1992`) on 400 cells; `pam250mt` = Biopython `PAM250` on the 20 amino acids (parser check); `blosum62mt2` = 2 x BLOSUM62; PAM350 only from ClustalW's own table (no other local reference) - its numbers are exercised by the PAM rows of the binary oracle.
- MSA: `msa_matrix_clustalw_oracle.tsv` (60 sets x progressive / iterative / consistency): rows of an independent Python implementation whose UPGMA distances are `1 - identities / MIN(len)` of the port; C# equal on all 180 rows. NJ (single-precision) and quick distances through the MSA are covered by the invariants only (the port of the NJ arithmetic was verified on DNA in F56).

## ClustalW PROTEIN profile stage (audit A1-10c part b, B09 F63)
- Sources actually opened (sources.debian.org, clustalw 2.1+lgpl-7, `src/`): `multipleAlign/ProfileBase.cpp` + `ProfileBase.h` (all of it: `calcGapCoeff`, `calcVPenalties`, `calcResidueSpecificPen` with the Pascarella & Argos table `ACDEFGHKILMNPQRSTVYW` / `87 87 104 69 80 139 100 104 68 79 71 137 126 93 128 124 111 75 100 77`, `calcHydrophilicPen`, `localPenalty`, `percentId`; the member `gdist` is declared but never initialised), `multipleAlign/MyersMillerProfileAlign.cpp` (`profileAlign`: the percent identity of the two groups, `len1` / `len2`, `logmin` / `logdiff`, the protein gap coefficients, `gaps` sizing), `multipleAlign/{ProfileWithSub,ProfileStandard,MSA}.cpp`, `substitutionMatrix/SubMatrix.cpp` (`getProfileAlignMatrix` protein branch, `getMatrix`, `setUpCrossReferences`), `substitutionMatrix/matrices.h`, `general/UserParameters.cpp` (defaults `noVarPenalties = true`, `gapDist = 4`, `hydResidues = GPSNDQEKR`, `setHydResidues` with `MAXHYDRESIDUES = 9` from `general/clustalw.h`, `PrfScaleValues` with `float scale / intScale`), `interface/CommandLineParser.cpp` (`-MATRIX`, `-NOPGAP`, `-NOHGAP`, `-NOVGAP` = `setNoVarPenalties(false)`, `-HGAPRESIDUES`, `-GAPDIST`, `-ENDGAPS` = `setUseEndGaps(false)`, `-NEGATIVE`, `-MAXDIV`). Thompson, Higgins & Gibson 1994 (NAR 22:4673; the gap-penalty heuristics) as the description of the same steps; Pascarella & Argos 1992 (Protein Eng 5:121) for the residue table. The whole source tree of the Debian package was also fetched and compiled (g++) for a planned ENDGAPS verification, but the sandbox refused to run the self-built binary - see the hand-off package A1-10d.
- Oracle: `scripts/alignment/generate_clustalw_protein_msa_oracle.py OUT 1000` -> `clustalw_protein_msa_oracle.tsv` (1 000 sets; family / hydrophilic-block / low-complexity / repeat / tiny / duplicate sets, B Z X O U residues; 333 default runs; `-MATRIX` 216, `-PWMATRIX` 92, `-GAPOPEN` 197, `-GAPEXT` 187, `-MAXDIV` 128, `-NOWEIGHTS` 74, `-NEGATIVE` 103, `-QUICKTREE` 149, `-NOPGAP` 135, `-NOHGAP` 138, `-NOVGAP` 189, `-GAPDIST` 180, `-HGAPRESIDUES` 107, `-ENDGAPS` 54, `-PWGAPOPEN` 85, `-PWGAPEXT` 88) = rows of the installed clustalw 2.1; the intermediate state (opening / extension penalty of every column and every profile column of every profile alignment) of the binary was compared with `-DEBUG=5` to localise each divergence during development.
- Findings: (1) the series matrix, the length / matrix-average scaling and the three per-column penalty families reproduce every row of the binary; (2) `-NOVGAP` is `setNoVarPenalties(false)`, i.e. it enables the variable-region penalties (the option name is inverted w.r.t. its effect); (3) `PrfScaleValues.scale` is a `float`, so `scale * matAvgScore * (float)(gapOpen / (logdiff * logmin))` is single precision (a case of the oracle gives 175, double precision 174); (4) `MAXHYDRESIDUES` = 9 cuts a longer `-HGAPRESIDUES`; (5) the lowest cell of the matrix shift is searched over `j = 1..i` only (column A is excluded); (6) `ProfileBase::gdist` is uninitialised on the `-NOVGAP` two-sequence > 60 % path: the binary's values were 0, 4, 5, 16, 20 or 30 depending on the heap, reproduced by forcing (see the algorithm doc §4.6); (7) PAM350 (`pam350mt`) equals the EMBOSS data file EPAM350 on all 529 cells (this supersedes the F62 remark that no local reference exists).
- Matrix series tables (`clustalw_protein_series_reference.tsv`): `blosum40mt` / `blosum45mt` / `blosum80mt` / `pam20mt` / `pam60mt` / `pam120mt` / `pam350mt` = EMBOSS EBLOSUM40/45/80, EPAM20/60/120/350 on all 23 x 23 cells; `blosum62mt2` = 2 x EBLOSUM62 (= 2 x Biopython BLOSUM62 on the shared alphabet); `gon40mt` / `gon80mt` / `gon120mt` / `gon160mt` / `gon350mt` (Gonnet, Cohen & Benner 1992 at other PAM distances) have no independent local reference (matrices.h only; Biopython carries GONNET1992 = 250 PAM only) - their use is verified by the oracle rows (the GONNET series is the default).


## ClustalW user-defined matrix files (audit A1-10d file part, B09 F64)
- Sources actually opened (sources.debian.org, clustalw 2.1+lgpl-7, `src/substitutionMatrix/SubMatrix.cpp`, `interface/CommandLineParser.cpp`, `general/UserParameters.cpp`, `general/clustalw.h`, `general/Utility.cpp`): `SubMatrix::readUserMatrix` (header = alphabetic characters of the first non-comment line, `Xref` over `aminoAcidCodes = "ABCDEFGHIKLMNOPQRSTUVWXYZ-"`, rows: optional alphabetic label token, `numargs < maxRes` check, per-row `.` scale detection, `ix != k + 1` row count, `#` / `!` terminator), `readMatrixSeries` (`CLUSTAL_SERIES` prefix test, `MATRIX %d %d %s`, limits 0-100, `ulimit > llimit`, `MAXMAT` = 10), `getMatrix` (`Xref` scatter, `matrixAvgScore`, the shift loop over `Xref`), `getProfileAlignMatrix` (the `userSeries` and `userMat` branches, scale `0.5 + (pcid - llimit) / ((ulimit - llimit) * 2.0)` into the float `PrfScaleValues.scale`), `getPairwiseMatrix` (`pwUserMat`, intScale 100, never shifted), `SubMatrix::commentline`, `Utility::lineType`, the `-MATRIX` / `-PWMATRIX` handlers of `CommandLineParser.cpp` (file name = matrix name, `exit(1)` on a reading error).
- Oracle: `scripts/alignment/generate_clustalw_user_matrix_oracle.py OUT 240` -> `clustalw_user_matrix_oracle.tsv` (240 sets with the matrix files embedded as JSON; matrices from Biopython 1.88 BLOSUM45 / 50 / 62 / 80 / 90, PAM30 / 70 / 250, GONNET1992, JONES, DAYHOFF, MDM78, RISLER, JOHNSON and random matrices; layouts NCBI / EMBOSS, unlabelled, integer / `%.1f` / `%g` rows, tab separated, header subsets, shuffled residue order, junk above the diagonal, trailers; series of 2-4 entries; quirk matrices with J / lower-case letters / 2-3 residues): every output row and 2 000+ printed pairwise `Score:` lines of the installed `clustalw -TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT -MATRIX=file -PWMATRIX=file` equal; development run (seed 31337, 800 sets): 2 differences, both `-NOVGAP` (the documented uninitialised `gdist`).
- Findings: (1) ClustalW requires exactly *letters + 1* rows, i.e. a matrix file without a non-letter header column (the `*` of the NCBI files) and its row is rejected "wrong format in matrix" (verified with the binary); (2) only the lower triangle is read; (3) the x10 scale of a floating-point matrix is decided per row by the first value token (`%g`-formatted whole numbers lose it); (4) a header letter ClustalW does not know (J, lower case) is only a warning and misreads the stored values of the following rows (the profile matrix of such a file was reproduced cell by cell: the profile rows of a `-DEBUG=5` run of a 23-residue + J matrix with unique cell ids decode to the same permutation as the port); (5) a two-residue matrix divides by zero in `matrixAvgScore` (INT_MIN on x86) and still runs; (6) a series with no matching range uses the first entry with a scale outside [0.5, 1].

## ClustalW nucleotide matrices (audit A1-10d part 3, B09 F68)
- Sources actually opened (sources.debian.org, clustalw 2.1+lgpl-7, `src/`): `substitutionMatrix/SubMatrix.cpp` (`setUpCrossReferences`: `defaultDNAXref` over `nucleicAcidOrder = "ABCDGHKMNRSTUVWXY"`; `getPairwiseMatrix` DNA branch: `swgapdnamt` / `clustalvdnamt` / `pwUserDNAMat`, `gapOpenScale` 0.6667 and `gapExtendScale` 0.751 for CLUSTALW, `getMatrix(..., negFlag = true, intScale)`, the transition writes on the residue codes (0,4), (2,11), (2,12); `getProfileAlignMatrix` DNA branch: scale 1.0 / 0.66 (`PrfScaleValues.scale` float), `getMatrix(..., _negMatrix, 100)`, the transition writes through `(*_matXref)[0 / 4 / 2 / 11 / 12]`; `getMatrix`; `getUserMatFromFile` (DNA: single matrix via `readUserMatrix`, `getUserMatAddress` / `getUserXrefAddress` per (residue type, stage)); `readUserMatrix`; `setCurrentNameAndNum`), `substitutionMatrix/matrices.h` (`swgapdnamt`, `clustalvdnamt`), `interface/CommandLineParser.cpp` (`-DNAMATRIX` / `-PWDNAMATRIX`: `iub`, `clustalw`, else `getUserMatFromFile(file, DNA, ...)` then `setCurrentNameAndNum(file, 3, Protein, ...)`), `interface/InteractiveMenu.cpp` (`readMatrix`: built-in items call `setCurrentNameAndNum(menu.opt[i - 1].string, ...)` after the loop, user item `setCurrentNameAndNum(userFile, option, alignResidueType, alignType)`), `pairwise/FullPairwiseAlign.cpp` (gap formulas, `forwardPass` / `reversePass` / `diff`), `multipleAlign/MyersMillerProfileAlign.cpp` (`gapcoef = (int)(100.0 * GAPOPEN * scale)`).
- Oracle: `scripts/alignment/generate_clustalw_dna_matrix_oracle.py OUT 360` -> `clustalw_dna_matrix_oracle.tsv` (installed `clustalw` 2.1; 360 sets: 125 command-line runs `-TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT -DNAMATRIX / -PWDNAMATRIX` = IUB / CLUSTALW (35) / file (75) or a malformed file (15: exit 1, no output), 235 interactive-menu runs `clustalw -INTERACTIVE -TYPE=DNA` through stdin (user matrix for the profile stage, the pairwise stage or both; gap penalties, transition weight, delay-divergent, negative-matrix and fast / slow toggles through the menu items); rows, `Sequences (i:j) Aligned. Score:` and `Group g: Sequences: n Score: s` / `Delayed` lines recorded). Result: **every row, 3 485 pairwise scores and 1 393 group lines equal** (the 29 non-crashing short-header runs: 28 equal, 1 differs - undefined behaviour below).
- Findings (all measured on the binary): (1) `clustalw -DNAMATRIX=file -PWDNAMATRIX=file` has no effect on a DNA alignment (rows and every printed score equal the run without the options; the matrix of the file is read and validated - a malformed file / a series file ends the program with exit status 1 and no output); (2) the interactive menu applies the file (e.g. a 5-letter ACGTN matrix with +1 / -3 changes the printed pairwise scores from 90 / 94 / 92 / 78 / 82 / 86 to 33 / 70 / 92 / 32 / 31 / 62 and the group scores from 865 / 826 / 728 to 188 / 151 / 113); (3) `-DNAMATRIX=CLUSTALW` changes the group scores (865 -> 450 on the same set) through the 0.66 gap scale and the identity matrix; (4) the menu's sequence reader classifies a set with >= 15 % non-ACGTUN letters as PROTEIN unless the type is given: the first oracle run recorded protein alignments for such sets (3 pairwise scores of set 116 differed in the port); the generator now starts the menu as `clustalw -INTERACTIVE -TYPE=DNA`; (5) undefined behaviour: the transition writes through `Xref[0 / 2 / 4 / 11 / 12]` with a header of fewer than 13 residues write outside `matrix[NUMRES][NUMRES]` (Xref -1): 8 of 37 such menu runs crashed (SIGSEGV; not in the oracle), 28 of the 29 others equal the port (which drops the writes); a matrix with fewer than three residues divides by zero in `getMatrix` (`av1 /= (maxRes * maxRes) / 2`) - not modelled.
- Contract: `UserDnaMatrixUse.Applied` (default) = the menu behaviour, `CommandLine` = what the command-line option really does; built-in matrices as `DnaMatrix` / `PairwiseDnaMatrix`.

## ClustalW end-gap options and the program's abort (audit A1-10d part 2, B09 F67)
- Sources actually opened (sources.debian.org, clustalw 2.1+lgpl-7, `src/`): `multipleAlign/ProfileBase.cpp` (`calcGapCoeff`: the `useEndGaps || endGapPenalties` test around the end-gap skipping loops; `profile[0][GAPCOL] = gapCoef` / `LENCOL = lencoef` when `endGapPenalties`, else the profile ends get 0), `multipleAlign/MyersMillerProfileAlign.cpp` (`openPenalty1/2`, `extPenalty1/2`, `gapPenalty1/2`: `if (!endGapPenalties && (i == 0 || i == prfLength1)) return 0`; `progDiff` returns the optimum, `profileAlign` returns `score / 100` with `int` division, `maxRes == 0 || maxRes == -1` returns -1), `multipleAlign/MSA.cpp` (`multiSeqAlign`: first pass `if (score < 0) return -1`, second pass ignores the score, `return _numSeqs`), `Clustal.cpp` (`align`: `if (count <= 0) return;` - no output written), `interface/CommandLineParser.cpp` (`-ENDGAPS`: `setUseEndGaps(false)`), `general/UserParameters.cpp` (`useEndGaps`, `endGapPenalties` default false), `interface/InteractiveMenu.cpp` ("Toggle End Gap Separation").
- Oracle: Docker package `docs/Validation/review-2026-09/handoff/A1-10d` (run by the user, x86-64, Debian 12 g++ -O2; the returned archive `returned/A1-10d-out.zip` is committed): `clustalw` 2.1+lgpl-7 built from the Debian source with ONE edit (the `-ENDGAPS` handler reads `CW_USEENDGAPS` / `CW_ENDGAPPEN`), 340 sets (120 DNA of `clustalw_msa_full_oracle.tsv`, 220 protein of `clustalw_protein_msa_oracle.tsv`, no `-NOVGAP`) x 4 modes (control, useendgaps, endgappen, both); `control_check.txt`: control rows equal the installed binary on 340 / 340 sets; `scripts/alignment/build_clustalw_end_gaps_oracle.py` merges sets and returned rows into `clustalw_end_gaps_oracle.tsv`.
- Findings: (1) `useendgaps` and `endgappen` change the alignment of many sets (the oracle is discriminating: see the test); (2) in the modes `endgappen` / `both` 34 sets (31 DNA + 3 protein; ids in the TSV) produced NO output file (rc 0): the diagnosis (a trace of the port with a score hook, one example re-verified here) shows a negative profile score in the FIRST pass, `MSA::multiSeqAlign` `score < 0` -> `return -1`, `Clustal::align` `count <= 0 -> return`; a negative optimum is impossible with free end gaps (the all-end-gap alignment costs 0), so the abort is specific to `EndGapPenalties`; the score is `optimum / 100` with truncating `int` division, so optima in (-100, 0) do not abort; (3) the second pass (delayed sequences) never aborts; (4) the five sets (D113, P116, P221, P61, P68, all `-QUICKTREE`) that differed in the diagnostician's scratch harness equal the oracle in the repo's tests in all four modes - the differences came from the scratch harness's option mapping (-KTUPLE sets the pair gap to ktuple + 3 / + 4, as in the committed oracle tests), not from the port.
- Contract: the abort is `ClustalWAbortedException : InvalidOperationException` (`Group` = 1-based number of the first-pass profile alignment, `Score` the negative value); the library never returns an alignment that ClustalW itself does not produce.

## MUSCLE 3.8.1551 oracle and the remaining options (audit A1-10c a + a2, B09 F69 / F70)
- **Oracle (F69).** The user ran the Docker package `handoff/A1-10c-muscle` (`debian:12.7-slim`, linux/amd64; the unmodified Debian `muscle3` 3.8.1551-2 source built with the upstream flags): `ENV.txt` MUSCLE v3.8.1551 by Robert C. Edgar, x86_64, binary sha256 `ab1a833e93c972e627c9e1802d1d62cb345d746caff26ef763b49b0b05e08b67`, `run.log` "420 sets run; failed: 0". `scripts/alignment/generate_muscle3_oracle.py` writes `TestData/Alignment/muscle3_oracle.tsv`; `MultipleAlignMuscle` equals it on **420 / 420 sets, rows and output order** (`SequenceAligner_MuscleOracle_Tests`, 841 cases). No difference, so no stage was changed. Platform caveat: 1-ulp differences of the single-precision `logf` between libm versions are possible in principle; none occurs here.
- **Options without an oracle yet (F70).** Source files read (text only, never executed): see the F70 entry of `B09.md` (`options.cpp`, `params.cpp`, `muscle.cpp`, `domuscle.cpp`, the refinement, objective score, weight, distance, clustering, tree, profile and mode files). The source decides the following, all reproduced: `VER_3_52 = 0` (diagonals compiled out), `HYDRO` undefined, `DOUBLE_AFFINE 0`; `-maxiters` default 8; leaf-index vs sequence-id weights after the tree refinement and for `-usetree`; `ReadMx` reads `HeadingCount - 1` columns; `DistKbit20_3` skips every second residue and counts the union of the bit vectors; `DistKmer20_3` returns a similarity; `RootByMinAvgLeafDist` swaps the two averages of its comments; `SetPPScore` replaces a protein score by SPN for nucleotides and SPN by LE for proteins; `-ppscore/-sw/-refinew` initialise the score matrix only for nucleotides and `-profdb` not the alphabet (the program crashes for proteins).
- **The second oracle round (F76).** `handoff/A1-10c-muscle-r2` (638 deterministic sets: 449 align, 36 profile, 35 refine, 22 spscore, 23 maketree, 27 ppscore, 18 refinew, 23 sw, 5 profdb) was run by the user (`ENV.txt`: MUSCLE v3.8.1551 by Robert C. Edgar, x86_64, Debian 12, binary sha256 `ab1a833e93c972e627c9e1802d1d62cb345d746caff26ef763b49b0b05e08b67`; `run.log` "638 sets run; non-zero or no output: 31"); the archive `returned/A1-10c-muscle-r2-out.zip` is committed alone (10ac6c8b) and was parsed as untrusted data (`python3 -I`), never executed. `scripts/alignment/generate_muscle3_r2_oracle.py` writes `TestData/Alignment/muscle3_r2_oracle.tsv` (630 rows; 8 `-matrix blosum62nostar` rows are `RUNNER_ERROR` in the archive: the harness could not decode the program's 0xff console bytes, so they are re-run in round 3). First comparison: 602 / 639 test cases equal; the 37 failures are the 8 unrecorded sets and 29 sets with 6 causes, all defects of the port, fixed stage by stage: (1) `DistKbit20_3` / `CommonBitCount`: `b = *p1 | (*p2 << 8)` with `b & 0x101` counts the bits set in either byte for the shifts 0-7 and the bits of the SECOND vector again for the shifts 8-15 (3 sets); (2) `-brenner` after the tree refinement: `RealignDiffsE` makes the root MSA with `MakeRootMSABrenner` again, rows in input order (10 sets); (3) `ProgressiveAlignE` reads `Weights[uId]` of an array filled by leaf INDEX: the port used the leaf node (differs for `-usetree`; 5 sets); (4) `printf("%g", -0.0)` = `-0`: the `-log(1.0)` distances of identical sequences give negative-zero edge lengths in the tree / distance files (6 sets); (5) `-maketree` never calls `SetPPScore`, so `-distance2` is never parsed (always `pctid_kimura`; 3 sets); (6) a one-sequence profile re-rooted with `-root1 midlongestspan | minavgleafdist` kills the program (`UnrootByDeletingRoot` writes through `GetLeft(root) = 0xFFFFFFFF`; 2 sets, modelled as `NotSupportedException`). After the fixes: 630 / 630 equal. Exceptions recorded as program behaviour: `-stable` ("not supported in this version of muscle"), `-distance2 edit` with neighbour joining ("MSADist::ComputeDist, invalid DISTANCE_11"), `-sw` on proteins with multi-sequence profiles ("DistCalcMSA: Invalid DISTANCE_1": `-distance1` is never parsed there), one-sequence neighbour joining ("Clust::SetLeafCount(1)"), `-weight2 henikoffpb` of `-profile` ("g_MuscleWeights = 0"), protein `-ppscore` / `-sw` / `-refinew` and every `-profdb` (segmentation fault).
- **The third round (audit A2-1, F89).** `apt-get install muscle3` works in the sandbox (Debian `muscle3` 3.8.1551-3 from the Debian archive, `muscle3 -version` = MUSCLE v3.8.1551, x86_64), so `handoff/A1-10c-muscle-r3` (60 sets: the 8 `-matrix blosum62nostar` sets, `-matrix` with DNA / RNA x8, one- / two-sequence corners x44 incl. the 10 refusals, `-le -sp -sv -spn` x8) was run locally with the package's own harness `run_sets.py` (console bytes decoded with `backslashreplace`, so nothing is lost) by `scripts/alignment/generate_muscle3_r3_oracle.py` -> `TestData/Alignment/muscle3_r3_oracle.tsv` (60 rows; `returned/ENV-local.txt`). First comparison 52 / 60 equal (+ 1 count of my own test): the 8 failures were `-matrix` on nucleotides (the port refused it) and `-maketree` of one sequence (the leaf is the root, so it is written as a group `(\ns0)\n;\n`, the port wrote `s0;`). Fixed from the source (`domuscle.cpp` ReadMx under the nucleotide tables, then `g_Alpha = ALPHA_Amino`, `g_PPScore = PPSCORE_SP`): `MuscleEngine._gNuc` separates the alphabet of the character tables from `g_Alpha`; 60 / 60 equal. **Fresh re-verification of the F69 / F76 claims**: `scripts/alignment/muscle3_fresh_sample.py --seed N` (seeds 7, 11-15; 300 random `align` sets each = 1 800 sets: protein / DNA / RNA, 1-9 sequences, random mixtures of iterations, distance1, cluster1/2, root1/2, weights, objective scores, term gaps, anchors, `-brenner`, score flags, gap parameters, `-matrix` blosum62 / pam250 / blosum62nostar on every alphabet; 1 602 runs, 70 segmentation faults, 128 + 6 program errors): 1 800 / 1 800 equal to the port (rows, written order, files, refusals), after two classification additions in the test mapping (the program's own `DistUnaligned, unsupported distance method` for the pctid_* distances of unaligned sequences and `DistKmer4_6 requires nucleo alphabet` are `ArgumentException` refusals). No defect of F69 / F76 was found. **Undefined behaviour found**: with `-matrix` on DNA / RNA `scorepp.cpp` reads `m_uSortOrder[4]`, which `SortCounts` never writes (it orders `g_AlphaSize = 4` letters) whenever a profile column holds all four nucleotides; Debian's binary segfaulted on 8 / 60 sampled runs and a `-O3` build of the same source (instrumented to print the value) read 0, 22 or 46 there and finished; under ASan the garbage is 0xBEBEBEBE. The port uses 0 (equal to the program on every run that finished, among them runs whose profiles did contain such columns, e.g. R644 and F7_11, F7_160, F7_282) and reports `MuscleRunResult.ReferenceBehaviourUndefined`; the test asserts that all 8 + 8 crashed runs are flagged.

## MUSCLE 3.8.1551 progressive alignment (audit A1-10c part a, B09 F65)
- Sources actually opened (sources.debian.org, muscle3 3.8.1551-2, flat tree of 183 files; downloaded into an isolated directory, never executed: the sandbox refused to compile it): `domuscle.cpp` (the command-line pipeline: alphabet, `SetPPScore`, `MHackStart`, `TreeFromSeqVect`, `ProgressiveAlignE`, `RefineTreeE`, `RefineVert(..., maxiters - 2)`), `muscle.cpp` (library `MUSCLE()`; same stages), `params.cpp` (defaults: `SetDefaultsLE` VTML_LA / -2.9 / -0.52 / normalised counts, `SetDefaultsSP` PAM200 / -1439, `SetDefaultsSV` VTML_SP / -300, `SetDefaultsSPN_DNA` NUC_SP / -400, `SetDefaultsSPN_RNA` -420 / center -300; `g_Distance1` Kmer6_6 / Kmer4_6, `g_dSUEFF` 0.1, `g_uMaxIters` 8, `g_bLow`), `alpha.cpp` / `alpha.h` (letter enums AX / NX, `Res` / `Wild`, DNA falls through into RNA), `fasta2.cpp`, `seq.cpp` (`FixAlpha`), `seqvect.cpp` (`GuessAlpha`), `mhack.cpp`, `fastdist.cpp`, `fastdistmafft.cpp` (`DistKmer6_6`, the `ResidueGroup` table whose entries are listed in the order B, Z, X while the enum is X, B, Z), `fastdistnuc.cpp` (`DistKmer4_6`), `fastclust.cpp`, `upgma2.cpp` (all of it incl. the nearest-neighbour fix-up and the `Test()` matrix), `phy.cpp` (`Tree::Create`, depth-first traversal, `LeafIndexToNodeIndex`), `phy3.cpp` (`FixRoot`: pseudo root = keep), `clwwt.cpp` (all of it incl. the worked example), `intmath.cpp` (`Normalize`), `msa.cpp` / `msa2.cpp` (`GetFractionalWeightedCounts`, `SetClustalWWeightsMuscle`, `NormalizeWeights`, `MSAFromColRange`), `cons.cpp` (`GetPctIdentityPair`), `profilefrommsa.cpp`, `aligngivenpath.cpp` (`AlignTwoProfsGivenPath`, `ProfScoresFromFreqs`, all `SetGaps*`, `SetFreqs`), `aligntwoprofs.cpp`, `glbalign.cpp` (`GlobalAlign` = `NWSmall`), `nwsmall.cpp` (all of it), `bittraceback.cpp`, `scorepp.cpp`, `termgaps.cpp`, `hydro.cpp` (`HYDRO` is never defined: the hydrophobicity heuristic is off in the profiles), `progressivealign.cpp` / `progalign.cpp`, `makerootmsa.cpp`, `refinetreee.cpp`, `difftreese.cpp`, `realigndiffse.cpp`, `msadistkimura.cpp`, `distcalc.cpp`, `refinevert.cpp` / `refinehoriz.cpp` (zero iterations leave the alignment untouched: `-maxiters 2` = draft + tree refinement), `vtml2.cpp`, `mpam200.cpp`, `nucmx.cpp` (tables transcribed by `scripts/alignment/transcribe_muscle3_matrices.py`; VTML_LA / VTML_SP / PAM200 / NUC_SP symmetric, 20 x 20 / 4 x 4 cells), `enums.h`, `muscle.h`, `types.h`, `intmath.h`.
- Edgar, R. C. (2004) NAR 32:1792 (the log-expectation score, k-mer distance, UPGMB, tree-dependent refinement) is consulted through the source; the formulae listed above are those of the code.
- Why no oracle run here: the installed `muscle` is the 5.x program (Debian package 1:5.1.0-1; `muscle -version` prints 5.2.linux64): a different pipeline, and the sandbox classifier refused to compile the downloaded 3.8.1551 source (`make` on the fetched tree), so a Docker package (`handoff/A1-10c-muscle`, user's Mac) builds the unmodified source and runs 420 deterministic sets (270 protein / 110 DNA / 40 RNA, `-maxiters 1|2`, `-le/-sp/-sv`, `-gapopen/-gapextend/-center`, `-cluster1/-cluster2`, `-sueff`, `-distance1`, `-maxtrees 0`). The result `out/muscle3_oracle.tsv` becomes the test oracle (README of the package); until then F65 is "implemented and checked against definitions", not "equal to MUSCLE".
- Source quirks reproduced on purpose (each is a deviation from what the paper / comments suggest, found by reading the code): (1) `ResidueGroup[]` order (X -> group 2, B -> 2, Z -> 0); (2) nucleotide wildcards in `GetFractionalWeightedCounts` compare the NX letter with `AX_R` (14) / `AX_Y` (19): X counts as G + A, N / R / Y ... as w/20 each of A, C, G, T; (3) the profile of a merged node is built from frequencies, the gap statistics follow the path (`SetGaps*`), the penalties are recomputed in double precision; (4) `AlignTwoProfs` passes `uLengthB` for A's length to `AlignTwoProfsGivenPath` (unused); (5) the single-column special cases of `NWSmall` (`uLengthA == 1` adds the close penalty of column 0); (6) `SetTermGaps` multiplies the zeroed terminal penalties by -1 (signed zero) and mutates the profile in place; (7) iteration 2 keeps the old weights of unchanged nodes and builds the new tree over the rows in draft order.

## T-Coffee 13.41 program: primary libraries, guide tree, whole pipeline (audit A1-13b, B09 F66)
- Program and source: as A1-13 (Debian `t-coffee` 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1; sources.debian.org `t_coffee_source/`: `t_coffee_lib/t_coffee.c`, `util_lib/util_constraints_list.c` (`method2pw_cl`, method files, `aln2constraint_list`, `seqpair2weight`, `add_entry2list`, `fork_read_n_constraint_list`, `fork_subset_produce_list`, `freeze_constraint_list`, `safe_dump_constraint_list`, `cl2distance_matrix`), `dp_lib/util_dp_drivers.c`, `util_dp_mm_nw.c`, `util_dp_sim.c`, `util_dp_suboptimal_nw.c` (`proba_pair_wise`, `ProbaMatrix2CL`), `util_dp_fasta_nw.c` (`ktup_dist_mat`), `evaluate.c` (`get_avg_matrix_mm`, `cw_profile_get_dp_cost`, `slow_get_dp_cost`), `evaluate_dirichlet.c` (`aln2prf`, dirichlet code), `util_lib/aln_convertion_util.c` (`get_seq_sim`, `make_group_aa`, `fix_seq_seq`), `util_lib/util.c` (`LOG_ADD`, `EXP`, `sort_int_inv`), `util_lib/util_make_tree.c` (`slow_nj_tree`, `fast_nj_tree`, `guide_tree`), `io_lib/tree_util.c` (`read_tree`, `reroot`, `calc_mean`, `insert_root`), `io_lib/io_func.c` (`read_matrice`, `get_string_type`), `matrices.h`, `coffee_defines.h`) - downloaded as text into an isolated directory, read, never compiled.
- Findings that change F58's description: (1) `read_matrice("dna_idmat")` returns the array `idmat` (+10 / **0**), not `dna_idmat` (+10 / -1); `slow_pair` / `lalign_id_pair` take GOP -10, GEP -1 from the method file for DNA and RNA (protein `slow_pair`: GOP `(int)(10 x get_avg_matrix_mm(blosum62mt))`, GEP -1; protein `lalign_id_pair`: BLOSUM50, -10, -4); F58's fitted "gop -11" is `g + h k` with g = -10, h = -1; (2) the Myers-Miller routine of the program is not optimal for its own model: the sub-problem with one residue left charges the open penalty on both flanks while the first row / column start with `tb = te = 0` (TTTCAACTC vs TTTTCAACTC -> `T-TTCAACTC`); the oracle values are the program's; (3) the residues of a SIM alignment are mapped to the full sequence by a second Myers-Miller alignment (`fix_seq_seq`; identity matrix, -4 / -1, also for proteins because the type of `aln2seq`'s sequences is never set), so the SIM library holds pairs outside the local alignments; (4) two methods are not merged as a union: `freeze_constraint_list` / `safe_dump_constraint_list` return list positions, so entries inserted before older entries are lost and older entries are returned again (counts are doubled); the in-memory list is asymmetric and the relaxation / DP run on it; the result does not depend on `-n_core`; (5) `get_string_type` decides DNA at > 98 % letters of `AGCTUNRYMKSWHBVD`; (6) the default tree is `ktup_dist_mat` + `slow_nj_tree` (`fast_nj_tree` from 100 sequences) + `reroot`, not the paper's; (7) input names of which one is a prefix of another (s1 / s10) change the tree for 11 or more sequences - the oracle uses q000, q001, ...
- Oracles (`scripts/alignment/generate_tcoffee_program_oracle.py`, `-n_core 1 -no_warning`, seed 7): `tcoffee_program_library_oracle.tsv` - 820 sets of 2-4 sequences (DNA: slow_pair 90, lalign_id_pair 90, proba_pair 90, slow_pair,lalign_id_pair 120, lalign_id_pair,slow_pair 60, proba_pair,slow_pair,lalign_id_pair 60; protein: slow_pair 70, lalign_id_pair 70, proba_pair 70, slow_pair,lalign_id_pair 100), the whole `-lib_only -out_lib` (pairs, weights, counts): **820 / 820 equal**; `tcoffee_program_msa_oracle.tsv` - 1 480 sets of 2-8 sequences (DNA plain 260, protein plain 260, DNA slow_pair,lalign_id_pair 260, protein slow_pair,lalign_id_pair 260, DNA slow_pair 80, lalign_id_pair 80, lalign_id_pair,slow_pair 80, proba_pair,slow_pair,lalign_id_pair 80, protein slow_pair 60, lalign_id_pair 60): guide tree parsed from the program's "Group n" log and rows of `-output fasta_aln`: **1 480 / 1 480 trees and row sets equal**; `tcoffee_program_big_oracle.tsv` - 100 and 101 sequences (`fast_nj_tree`): equal. Development runs besides (a few thousand further random sets, all equal): libraries of every method and merge order, alignments incl. 1-6 residue sets, near-identical sequences, two-letter alphabets, `-n_core` 2 / 4 / 8 (the result does not depend on it), 99-110 sequences.
- Residual: RNA `proba_pair` (a U; separate model and `get_tot_prob`), the other `-method`s, `-tree_mode upgma`, `-quicktree`, `-seq_weight`, `-relax_lib` > 1 (registered A1-13c / A1-13d); the name-prefix quirk; a rooting tie that crashes the program.


## T-Coffee 13.41 program: RNA `proba_pair` and the command-line options (audit A1-13c / A1-13d, B09 F71-F75)
- Program and sources: as A1-13b (Debian `t-coffee` 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `t_coffee` 13.41 installed; sources.debian.org `t_coffee_source/`, text only, never compiled): `t_coffee_lib/t_coffee.c` (option table L1280-3400; order of the library operations L5000-5180: type, penalties, `choose_extension_mode`, `filter_lib`, `relax_lib`, `shrink_lib`, `evaluate_constraint_list_reference`, `weight_constraint_list`; `-inorder` sort L4706), `util_lib/util_constraints_list.c` (`seq2distance_matrix`, `cl2distance_matrix`, `weight_constraint_list`, `compute_t_coffee_weight`, `re_weight_constraint_list`, `evaluate_constraint_list_reference`, `filter_constraint_list`, `shrink_constraint_list`, method files L5380-7100, `@EP@` parameters L1096-1125), `dp_lib/util_dp_gotoh_nw.c` (`list2linked_pair_wise`, `cl2list_borders`, `idscore_pairseq`), `util_dp_suboptimal_nw.c` (`proba_pair_wise`, `get_tot_prob`, `get_emitPairs`), `util_dp_drivers.c` (`make_tree`, `method2pw_cl`), `util_lib/util_make_tree.c`, `util_lib/aln_convertion_util.c`, `io_lib/io_func.c` (`read_matrice`, type decision), `io_lib/reformat.c` (`read_seq_weight`), `matrices.h`.
- Findings that change F66 / F58: (1) the "name-prefix quirk" is the program's alphabetical sort of the sequences (`-inorder aligned`); (2) the sparse DP is not gap-free: `-gapopen` / `-gapext` act on it and the border cells carry `gep x 10 x position`; F58's "gap open = gap extend = 0" is the default, not a property; (3) `idscore` distances are not multiplied by 10; (4) UPGMA children are delivered in reverse order by `read_tree`.
- Oracles (`scripts/alignment/generate_tcoffee_options_oracle.py`, seed 7, `-n_core 1 -no_warning`): 4 715 sets in `tcoffee_program_options_oracle.tsv`: rna 430, relax / filter 500, gap model 560, tree / distance 730, weights 470, name order 560, matrix 675, method parameters 390, clustalw_pair 400. All rows equal (`SequenceAligner_TCoffeeOptions_Tests`).
- Measured as without effect on the rows (0 of 24 sets, DNA and protein): `-extend_mode slow_triplet`, `-extend`, `-evaluate_mode`, `-shrink_lib`, `-tg_mode`, `-ktuple`, `-weight no` (without `-seq_weight`).
- Not ported (registered DOABLE, see B09.md): the fasta heuristics (`fast_pair`, `-quicktree`, `-distance_matrix_mode very_fast`), `clean_slow_pair`, `hash_pair`, `biphasic_pair`, `exon*_pair`, `promo_pair`, `blastr_pair`, `lalign_rs_s_pair`, `subop*_pair`, `-dp_mode` other than `linked_pair_wise`, `-tree_mode cwph / kmeans`, the `slow` / `cscore` / `geometric` / `aln` / `random` distance modes. Blocked: methods that need BLAST / PDB / structure programs or a network database (`blastp_msa`, `blastp_o2a`, `tblastx_msa`, `align_pdbpair`, `lalign_pdbpair`, `sap_pair`, `TMalign_pair`, `mustang_pair`, `dali_pair`, `fugue_pair`, `thread_pair`, `-mode 3dcoffee / expresso`); `mafft_pair` / `muscle_pair` fail inside the installed `t_coffee` (rc 1: its fixed command line does not suit the installed `mafft` / `muscle` 5.x).
- Residual differences of the committed file (3 of 4 715 sets, each in `KnownDifferences` of `SequenceAligner_TCoffeeOptions_Tests`): protein `-nomatch 50 -gapext -1` (7 sequences); protein `-seq_weight t_coffee -distance_matrix_mode idscore -tree_mode upgma` (8 sequences, minimal 5; uniform weights 0.04 / 0.05 reproduce it while 0.03 / 0.045 / 0.055 do not; the weights equal the program's); `clustalw_pair` on the protein pair YAGR / NSMYC (installed clustalw `---YAGR / NSMYC--`, port `-YAGR / NSMYC-`).
- Program core dump (reproduced as `TCoffeeAbortedException`): `lalign_id_pair@EP@GOP@-5` on 3 of 25 committed sets; SIM's script array `S[min(M, N) x 5 / 4 + 1]` (`util_dp_sim.c` L326) overflows its glibc chunk (usable bytes ((4 x capacity + 8 + 15) / 16 x 16) - 8); the rule predicts the 3 crashing and the 22 surviving sets and 1 of 1 crash in 240 random pairs.
- Blocked proofs (installed program, this sandbox): `t_coffee -method blastp_msa` stops with "If legacy_blast.pl is NOT installed obtain a copy" (no BLAST database, NCBI unreachable); `sap_pair`: "SAP is not installed"; `fugue_pair`: "fugueali is NOT installed"; `TMalign_pair`, `mustang_pair`, `dali_pair`: "Method TMalign cannot be applied to [q000 vs q001], proba_pair will be used instead" (they need PDB structures); `mafft_pair` / `muscle_pair`: the program stops with rc 1 on plain sequences.

## Audit round 2 (A2-38 / A2-39 / A2-40 / A2-41 / A2-42, B09 F81-F86): input contract, size guards, parallel exceptions, cancellation
- Sources actually opened: T-Coffee 13.41 `util_dp_suboptimal_nw.c` (Debian `t-coffee` 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, local copy from A1-13b): `forward_proba_pair_wise` / `backward_proba_pair_wise` compute `int l = (seq1Length + 1) * (seq2Length + 1) * NumMatrixTypes` and `vcalloc` / `vrealloc` it - 32-bit arithmetic, no explicit limit; MUSCLE 3.8.1551 (sources.debian.org `muscle3` 3.8.1551-2, `nwsmall.cpp`, `muscle.h`, `globals.cpp`, `glbalign.cpp`): `NWSmall` allocates `char **TB` with one traceback byte per cell and has no size check; `grep -i maxmb` over `muscle.h` and `globals.cpp` finds nothing (the `-maxmb` option of the later 3.8.31 line is accepted and recorded by `MuscleOptions.MaxMb`, no effect). Neither program refuses a long input: both fail by 32-bit wrap-around / out of memory, so the port refuses it with an `ArgumentException` before allocating: T-Coffee `proba_pair` 2^28 floats per matrix (`TCoffeeProbaMaxMatrixFloats`; the program wraps at 2^31), MUSCLE 1e9 cells (`MuscleProfileMatrixMaxCells`; 2 x 30 kb = 9e8 cells = 954 MB measured by the auditor). ClustalW is linear space (Myers-Miller) with O(k^2 L^2) time and no limit of its own: no size guard, the call is cancellable instead.
- Empty members: ClustalW 2.1 cannot read an empty sequence, MUSCLE 3.8 ignores empty FASTA records ("Sequence N is empty"), `t_coffee` stops; the ports keep refusing, with ONE message (`EmptyMemberMessage`) that names the aligners that accept empty members. `ClustalWUserMatrix.Parse` keeps `FormatException` (a parser contract documented in its XML doc; the MCP tool wraps it).
- Parallel pair stage: `Parallel.For` / `ForEach` wrap any body exception in `AggregateException`; the star aligners now run through `SequenceAligner.ParallelFor`, which records the exception of the lowest failing index (the one the sequential path would throw first) and rethrows it with `ExceptionDispatchInfo` (original type and stack).
- Cancellation: `System.Threading.Tasks.ParallelOptions.CancellationToken` is the .NET convention for a cancellable parallel / compute option set; the option records of the aligners follow it (`MsaOptions`, `ClustalWOptions`, `ClustalWProteinOptions`, `MuscleOptions`, `TCoffeeProgramOptions` get `CancellationToken`, ignored by the semantics of the options), `SequenceAligner.WithCancellation` / `BeginCancellationScope` cover every other overload, the dynamic programs poll an `AsyncLocal` token at their row / iteration boundaries (it flows into `Parallel.For` workers).

## Audit A2-14 (B09 F92): ClustalW `-TYPE` guess and `-CLUSTERING=UPGMA`
- **Opened:** Debian `clustalw` 2.1+lgpl-7 sources (sources.debian.org, neutral User-Agent, `src/**.cpp|.h`, 120 files): `fileInput/FileReader.cpp` (`readFile`: the flag is set at `i == 1`), `alignment/Sequence.cpp` (`checkDNAFlag`: A C G T U N >= 85 %), `tree/UPGMA/UPGMAAlgorithm.cpp`, `Node.cpp`, `RootedGuideTree.cpp`, `RootedClusterTree.cpp`, `tree/TreeInterface.cpp` (`getWeightsAndStepsFromDistMatUPGMA`), `general/SymMatrix.h` (`getIndex`, `makeSimilarityMatrix`), `multipleAlign/Iteration.cpp` (head).
- **Numbers:** `-TYPE` omitted: 500 runs of the installed clustalw 2.1 (the "Sequence n: name len bp|aa" unit; 176 DNA, 324 protein incl. 1-sequence sets and sets around the 85 % threshold): 500 / 500 equal. `-CLUSTERING=UPGMA`: the installed binary exits 139 (segmentation fault) on every `-ALIGN` / `-TREE` run (2, 3, 4 sequences); `-ITERATION=TREE | ALIGNMENT -NUMITER=3` runs but leaves a 4-sequence set unchanged. The audit sandbox compiled the unmodified source at `-O0` + AddressSanitizer: the UPGMA alignment ran to the end ("Group 1: Sequences: 2 Score:228", "Group 2: Sequences: 3 Score:199", "Alignment Score 216" on three 13-mers); the next self-built run was refused by the sandbox classifier, so the oracle is the user-run package `handoff/A2-14-clustalw-upgma`.
- **Findings:** the type decision of ClustalW 2.1 looks only at the second sequence; `RootedClusterTree::treeFromDistMatrix` stores no alignment step for two sequences (a defect of the program: the two sequences are not aligned); `Node::merge` overwrites the height of BOTH children with the height of the join (so `doWeightCalc` weights a node by its parent's height); UPGMA's similarities are `100 - 100 d` with no clamp.

## 2026-11-19 audit round 2 (A2-36, B09 F100): MCP tools for the public stage ports

The `multiple_align_*` tools run whole pipelines; the stage ports below had no tool. Each new tool is a thin wrapper that delegates to `SequenceAligner` (no re-implementation) and its tests lock values from the reference-checked library tests: `clustalw_pairwise_align` (ClustalW 2.1 `FullPairwiseAlign`: clustalw prints 40 / 38 / 31 for the globin fragments, `SequenceAligner_ClustalWProtein_Tests`; DNA / protein options `-PWGAPOPEN -PWGAPEXT -TRANSWEIGHT -PWDNAMATRIX -PWMATRIX`), `clustalw_quick_pairwise_align` (`-QUICKTREE`: 92.5926 / 27.7778 hand run), `tcoffee_slow_pair_align`, `tcoffee_lalign_top_alignments` (scores 790 / 500 / 400 / 370), `tcoffee_primary_library`, `tcoffee_guide_tree` (Newick + leaf order), `tcoffee_ktup_similarity` (100 / 70), `sum_of_pairs_score` (`SumOfPairsScore` / `MatrixSumOfPairsScore`, linear / affine, matrix by name). `MuscleGuideTree` has a tool already: `muscle_msa_command` mode `maketree`. `GenomicAnalyzer.FindKnownMotifs` (positions only) is covered by `find_known_motifs` (`matches`).
