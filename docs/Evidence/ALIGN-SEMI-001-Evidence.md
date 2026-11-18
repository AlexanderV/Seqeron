# Evidence: ALIGN-SEMI-001

**Test Unit ID:** ALIGN-SEMI-001
**Algorithm:** Semi-Global Alignment (Fitting / Query-in-Reference)
**Date Collected:** 2026-02-01
**Last Audited:** 2026-03-07

---

## 1. Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Sequence alignment | https://en.wikipedia.org/wiki/Sequence_alignment | 2026-02-01 |
| Wikipedia: Needleman–Wunsch algorithm | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | 2026-02-01 |
| Rosalind: Finding a Motif with Modifications (SIMS) | https://rosalind.info/problems/sims/ | 2026-02-01 |
| Rosalind: Semiglobal Alignment (SMGB) | https://rosalind.info/problems/smgb/ | 2026-02-01 |
| Brudno et al. (2003): Glocal alignment | doi:10.1093/bioinformatics/btg1005 | 2026-02-01 |

---

## 2. Algorithm Definition (from sources)

### 2.1 Definition

Semi-global alignment (also known as "glocal" or "ends-free" alignment) is a **hybrid method** that searches for the best possible partial alignment of two sequences. It combines features of both global and local alignment, allowing alignment to start and/or end at any position in one or both sequences.

**Source:** Wikipedia (Sequence alignment):
> "Hybrid methods, known as semi-global or 'glocal' (short for global-local) methods, search for the best possible partial alignment of the two sequences (in other words, a combination of one or both starts and one or both ends is stated to be aligned)."

The implementation uses the **fitting alignment** variant as defined by Rosalind:

**Source:** Rosalind (SIMS):
> A fitting alignment is "an alignment of a substring of s against all of t."

### 2.2 Primary Use Cases

1. **Overlap alignment**: When the downstream part of one sequence overlaps with the upstream part of another sequence.

2. **Short read alignment (fitting)**: When one sequence is short (e.g., a gene) and the other is long (e.g., a chromosome). The short sequence should be globally aligned, but only a local/partial alignment is desired for the long sequence.

**Source:** Wikipedia (Sequence alignment):
> "Another case where semi-global alignment is useful is when one sequence is short (for example a gene sequence) and the other is very long (for example a chromosome sequence). In that case, the short sequence should be globally (fully) aligned but only a local (partial) alignment is desired for the long sequence."

### 2.3 Algorithm Variants

| Variant | Free End Gaps | Typical Use Case | Rosalind Problem |
|---------|---------------|------------------|------------------|
| Query-in-reference (fitting) | Start and end of reference | Short read mapping, primer alignment | SIMS |
| Overlap | End of seq1, start of seq2 | Sequence assembly | OAP |
| Full semi-global | All four ends free | General substring finding | SMGB |

**Source:** Rosalind (SIMS, SMGB); Wikipedia (Sequence alignment).

**Design Decision:** The implementation uses the query-in-reference (fitting) variant. This is a deliberate design choice selecting one well-defined member of the semi-global family, corresponding to the Rosalind SIMS problem.

**Reference-implementation equivalents (2026-09 review):** Biopython 1.88 `PairwiseAligner(mode="global", end_insertion_score=0)` (target = sequence 1 = query) and parasail 1.3.4 `sg_dx` (s1 = query; "database begin/end gaps free"). Seqeron's `SemiGlobalAlign` (linear) and `SemiGlobalAlignAffine` (Gotoh affine; Biopython open_gap_score = GapOpen + GapExtend, parasail open = −(GapOpen + GapExtend)) matched both on 2428 random/edge cases (scores; every returned alignment is in Biopython's co-optimal set). Only overhanging *reference* residues are free: a query overhang is charged (ACGTAA / CCACGT, 1/−1/−1 → 2; an overlap aligner would give 4).

---

## 3. Algorithm Mechanics (from sources)

### 3.1 Matrix Initialization

| Alignment Type | First Row $F_{0,j}$ | First Column $F_{i,0}$ |
|----------------|---------------------|------------------------|
| Global (NW) | $d \cdot j$ | $d \cdot i$ |
| Local (SW) | $0$ | $0$ |
| Fitting (query-in-ref) | $0$ (free start gaps in reference) | $d \cdot i$ (query fully aligned) |

**Source:** Needleman–Wunsch algorithm (Wikipedia): "First row and first column are subject to gap penalty." For fitting: first row is zero per free end-gap convention (Rosalind SIMS definition).

### 3.2 Recurrence

$$F_{i,j} = \max\left(F_{i-1,j-1} + s(a_i, b_j),\; F_{i-1,j} + d,\; F_{i,j-1} + d\right)$$

No zero floor (unlike Smith–Waterman). Score can be negative. This follows from the Needleman–Wunsch recurrence with linear gap cost.

**Source:** Wikipedia (Needleman–Wunsch algorithm), "Advanced presentation of algorithm" section.

### 3.3 Traceback

For the fitting alignment variant:
- Traceback starts from $\max_j F_{m,j}$ — the **maximum score in the last row** (not the bottom-right cell)
- Traceback proceeds backward to $F_{0,*}$ (top row), ensuring full query coverage
- Remaining reference bases after the traceback endpoint are appended as gaps without penalty
- Ties: the smallest end column j (including j = 0, the query aligned entirely against gaps) wins; traceback ties diagonal > up > left

**Source:** Rosalind (SIMS) — fitting alignment definition implies the optimal score is the maximum over all possible endpoints in the reference. Since the query must be fully consumed (rows 0 to m), the traceback spans the full last row.

### 3.4 Complexity

- **Time**: $O(m \cdot n)$ where $m$ = query length, $n$ = reference length
- **Space**: $O(m \cdot n)$ for the scoring matrix

**Source:** Wikipedia (Needleman–Wunsch algorithm): "The time complexity of the algorithm for two sequences of length n and m is O(mn)."

---

## 4. Documented Corner Cases

| Case | Evidence Source | Expected Behavior |
|------|----------------|-------------------|
| Query embedded in reference | Fitting alignment definition (Rosalind SIMS) | Full query aligned, reference has free leading/trailing gaps |
| Identical sequences | NW global = fitting when lengths equal | Full alignment, score = m × match |
| Query at start/end of reference | Fitting alignment variant | Correct alignment with free trailing/leading reference gaps |
| All mismatches | NW recurrence (no zero floor) | Exact negative score (e.g., m × mismatch) |
| Mixed matches and mismatches | NW recurrence | Score = Σ(match scores) + Σ(mismatch scores) |
| Gap in optimal alignment | NW recurrence (up/left moves) | Score includes gap penalty; gap visible in aligned output |
| Null sequence | .NET API convention | Throws `ArgumentNullException` |

---

## 5. Key Invariants

| ID | Invariant | Evidence |
|----|-----------|----------|
| INV-1 | AlignmentType = SemiGlobal | Implementation contract |
| INV-2 | Aligned sequences have equal length | Alignment theory (all alignment types) |
| INV-3 | Query fully represented | Fitting alignment: `RemoveGaps(aligned1) == query` (Rosalind SIMS) |
| INV-4 | Reference is substring | Fitting alignment: `RemoveGaps(aligned2)` is substring of reference |
| INV-5 | Score = $\max_j F_{m,j}$ | Fitting alignment traceback from max of last row |
| INV-6 | `StartPosition2`/`EndPosition2` = fitted interval of the reference (first/last reference index consumed between the free end gaps; empty interval: End2 = Start2 − 1); Start1/End1 = 0 / m − 1 | Biopython 1.88 `PairwiseAligner(mode='global', end_insertion_score=0)`: leading/trailing gap-in-target columns of each co-optimal alignment (`align` enumeration); parasail 1.3.4 `sg_dx_trace` `end_ref`. Oracle run (2026-11-18, A1-4): 4090 random/edge cases (alphabets A, AC, ACGT; lengths 0..30; 10 scoring schemes incl. GapOpen 0 / GapExtend 0), 1382 semi-global rows: reported interval ∈ set of Biopython co-optimal intervals in 823/823 fully enumerated rows (32 capped at 3000 alignments), parasail score equal in all rows, `Global(query, ref[Start2..End2]).Score` = fitting score in all rows; 0 failures |

---

## 6. References

1. Wikipedia. "Sequence alignment." https://en.wikipedia.org/wiki/Sequence_alignment
2. Wikipedia. "Needleman–Wunsch algorithm." https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
3. Rosalind. "Finding a Motif with Modifications (SIMS)." https://rosalind.info/problems/sims/
4. Rosalind. "Semiglobal Alignment (SMGB)." https://rosalind.info/problems/smgb/
5. Brudno M, Malde S, Poliakov A, Do CB, Couronne O, Dubchak I, Batzoglou S (2003). "Glocal alignment: finding rearrangements during alignment." Bioinformatics 19 Suppl 1:i54-62. doi:10.1093/bioinformatics/btg1005

## Substitution-matrix scoring (audit A1-2, B09 F41)

Matrix-scored fitting alignment is verified against Biopython 1.88 `end_insertion_score = 0` and parasail 1.3.4 `sg_dx` (matrices BLOSUM62/PAM250/NUC.4.4) on 662 semi-global rows (331 linear + 331 affine, integer and fractional gaps): 662/662 equal (parasail on the integer-gap subset). EMBOSS has no one-sided free-end mode. See ALIGN-GLOBAL-001-Evidence.md.

## Ends-free and overlap variants (audit A1-3, B09 F44)

**Sources actually used (2026-11-18):** Biopython 1.88 `Bio.Align.PairwiseAligner` (installed; properties `left/right_insertion_score`, `left/right_deletion_score`, `end_gap_score`; source of `score()` vs `align()` inspected by experiment), parasail 1.3.4 (`nw`, `sg`, `sg_qb`, `sg_qe`, `sg_db`, `sg_de`, `sg_qx`, `sg_dx`, `sg_qb_db`, `sg_qb_de`, `sg_qe_db`, `sg_qe_de`; names mapped to flags by exhaustive search of the flag set that reproduces each variant), EMBOSS needle 6.6.0 default (`-endweight N`), Durbin et al. 1998 §2.3 (overlap/fitting/free end-gap variants, WebSearch snippets), Rosalind SIMS/SMGB (cited since the first review).

**Oracle** (`scripts/alignment/generate_endgaps_oracle.py`, 3 590 rows): Biopython `align().score` = Seqeron on 3 590/3 590 (all 16 policies x linear/affine, scalar and BLOSUM62/EDNAFULL with fractional gaps); co-optimal alignment membership 2 050/2 050; parasail equal on the 6 exact variants (nw, sg_qb, sg_db, sg_qx, sg_dx, sg_qb_db), never above the optimum on the 6 variants with a free end flag (23 lower rows: parasail does not take the border end cells (m,0)/(0,n); an exhaustive reference skipping them reproduces parasail 300/300); EMBOSS needle default 277/287 equal, 10 lower (embaln.c end-gap quirk, F41). Biopython `PairwiseAligner.score()` is not the oracle: it disagrees with its own `align().score` on 25 degenerate rows (one sequence entirely against gaps; `score('AC','T')` = −2 vs the alignment `AC-`/`--T` = 0 with `left_deletion_score = right_insertion_score = 0`), `align().score` equals an independently written exhaustive reference on every row.

**Conventions fixed by the oracle.** Free end gap = costs 0 whatever its length; every other gap is an ordinary gap (affine: GapOpen + k·GapExtend; the adjacent end run of the other sequence is charged); a path that consumes nothing of one sequence (zero-length overlap) is valid; ties: first maximal last-row column, then first maximal last-column row.

## Exact linear-space alignment (audit A1-6, B09 F46)

- **Sources actually used (2026-11-18):** Hirschberg 1975, Myers & Miller 1988 (CABIOS 4:11, PMID 3382986), Huang & Miller 1991 (Adv Appl Math 12:337, DOI 10.1016/0196-8858(91)90017-D), Huang-Hardison-Miller 1990 / SIM (WebSearch records only; papers not opened). Reference tools run locally: EMBOSS 6.6.0 `stretcher` (Myers-Miller linear-space global, end gaps penalised; gap of length k = gapopen + (k-1) gapextend, i.e. `GapOpen = -(G - E)`, `GapExtend = -E`) and `matcher` (SIM local, `-alternatives 1`), Biopython 1.88 `PairwiseAligner`, parasail 1.3.4.
- **Contract verified:** the linear-space result equals the full-matrix aligner row for row (score, rows, coordinates, type, ties): 4 750 random / tie-rich / long-indel / matrix cases with base cells 1-100 over all 16 end-gap policies, local, linear + affine; 2.6-3 kb pairs; the automatic switch on a 5 300 x 5 200 pair for all six aligners.
- **Biopython / parasail** (`scripts/alignment/generate_linear_space_oracle.py`, 1 280 pairs): scores equal on all rows; parasail equal on every row where a variant exists; returned alignments re-score independently and lie in Biopython's co-optimal set whenever <= 40 alignments exist.
- **EMBOSS stretcher / matcher oracle rows (160 + 160):** scores equal on all rows; EMBOSS rows re-score to their score under the Seqeron model; where Biopython finds exactly one optimal alignment (49 stretcher rows, 93 matcher rows) the EMBOSS and Seqeron rows (and matcher's 1-based coordinates) are equal; elsewhere ties are broken differently (55 / 160 and 99 / 160 identical rows).
- **Large EMBOSS cross-check (DNA, 10 % edits + 3 planted indels of 20-300 nt, EDNAFULL, gapopen 16 / gapextend 4):** `stretcher` score = `GlobalAlignAffine` and `matcher` score = `LocalAlignAffine` for 100, 500, 1 000, 2 000, 5 000, 20 000, 50 000 nt (14 / 14): 50 000 x 50 360: 191 903 / 192 525; 20 000 x 19 965: 75 454 / 75 477. Timing: stretcher 1.6 s / 10.8 s for 20 kb / 50 kb, this implementation 15.6 s / 89.5 s.
- **Memory (Release, process peak working set):** 20 000 x 19 994 affine global 111 MB (4.8 GB full matrices), local 118 MB, fitting 113 MB, linear global 57 MB; 100 000 x 100 029 (10^10 cells) affine global 185 MB / 391.7 s, linear global 73 MB / 153.9 s; 10 000 x 10 000 affine: 1 179 MB full matrix vs 103 MB.
- **Limits:** `LinearSpaceCellThreshold` 25 000 000 cells; `GeneralGapMaxCells` 4 000 000; `ProfileMatrixMaxCells` 50 000 000 (progressive / iterative / consistency MSA profile DP); see `docs/algorithms/Alignment/Linear_Space_Alignment.md` §4.
