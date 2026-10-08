# Evidence Artifact: ALIGN-GLOBAL-001

**Test Unit ID:** ALIGN-GLOBAL-001
**Algorithm:** Global Alignment (Needleman–Wunsch)
**Date Collected:** 2026-03-06

---

## Online Sources

### Wikipedia: Needleman–Wunsch Algorithm
**URL:** https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
**Accessed:** 2026-03-06

**Key Extracted Points:**

1. **Purpose:**
   - Algorithm for computing an optimal global alignment between two sequences.
   - One of the first applications of dynamic programming to biological sequences.
   - Aligns the entire sequences end-to-end.

2. **Scoring System (from "Choosing a scoring system"):**
   - Match: +1 (or other positive value depending on substitution matrix)
   - Mismatch: −1 (or negative value from substitution matrix)
   - Indel (gap): −1 (or specified gap penalty)
   - Quote: "The score of the whole alignment candidate is the sum of the scores of all the pairings."

3. **Linear Gap Penalty Model (from "Advanced presentation of algorithm"):**
   - Uses a similarity function $S(a,b)$ and a single linear gap penalty $d$.
   - Basis:
     - $F_{0,j} = d \cdot j$
     - $F_{i,0} = d \cdot i$
   - Recurrence:
     - $F_{i,j} = \max(F_{i-1,j-1} + S(A_i,B_j),\ F_{i,j-1} + d,\ F_{i-1,j} + d)$
   - The entry $F_{n,m}$ gives the maximum score among all possible alignments.

4. **Traceback Rules (from "Tracing arrows back to origin"):**
   - Diagonal path (top-left cell) → match or mismatch.
   - Horizontal path (left cell) → gap in the "side" sequence.
   - Vertical path (top cell) → gap in the "top" sequence.
   - "If there are multiple arrows to choose from, they represent a branching of the alignments."

5. **Complexity:**
   - Time: $O(nm)$ for sequences of length $n$ and $m$.
   - Space: $O(nm)$ for the full matrix.

6. **Example:**
   - Sequences: GCATGCG and GATTACA
   - Scoring: match = +1, mismatch/indel = −1
   - One possible optimal alignment: `GCATG-CG` / `G-ATTACA`
   - Score breakdown: `+−++−−+−` → $4 \times (+1) + 4 \times (-1) = 0$
   - Border initialization shown: first row is 0, −1, −2, −3, −4, −5, −6, −7.

7. **Gap Penalty Extensions (from "Gap penalty"):**
   - Affine gap penalties (gap-open + gap-extend) are mentioned as an *extension* to the basic model.
   - Quote: "The simple and common way to do this is via a large gap-start score for a new indel and a smaller gap-extension score for every letter which extends the indel."
   - The standard NW pseudocode uses **only the linear penalty** $d$.

### Wikipedia: Sequence Alignment (Global vs. Local)
**URL:** https://en.wikipedia.org/wiki/Sequence_alignment
**Accessed:** 2026-02-01

**Key Extracted Points:**

1. **Global Alignment:**
   - Attempts to align the entire length of all sequences (every residue).
   - Most useful when sequences are similar and of roughly equal size.
   - Can start and/or end in gaps.
   - Needleman–Wunsch is the canonical global alignment method.

2. **Multiple Optimal Alignments:**
   - Multiple traceback paths can exist for the same optimal score.
   - Common in sequences with repeated motifs or low complexity regions.

---

## Test Dataset: Wikipedia Example

| Parameter | Value |
|-----------|-------|
| Sequence 1 | GCATGCG |
| Sequence 2 | GATTACA |
| Match Score | +1 |
| Mismatch Score | −1 |
| Gap Penalty (d) | −1 |
| Optimal Score | 0 |
| One Optimal Alignment (seq1) | `GCATG-CG` |
| One Optimal Alignment (seq2) | `G-ATTACA` |
| Score Breakdown | 4 matches (+4) + 2 mismatches (−2) + 2 gaps (−2) = 0 |

### Wikipedia DP Matrix (Border Initialization)

Standard NW border for the example (d = −1):

| | - | G | A | T | T | A | C | A |
|---|---|---|---|---|---|---|---|---|
| **-** | 0 | −1 | −2 | −3 | −4 | −5 | −6 | −7 |
| **G** | −1 | | | | | | | |
| **C** | −2 | | | | | | | |
| **A** | −3 | | | | | | | |
| **T** | −4 | | | | | | | |
| **G** | −5 | | | | | | | |
| **C** | −6 | | | | | | | |
| **G** | −7 | | | | | | | |

Border values: $F(i,0) = -i$, $F(0,j) = -j$ — directly from Wikipedia.

---

## Affine gaps (added 2026-09 review)

- **Gotoh O (1982)** J Mol Biol 162:705-708 — three-matrix DP for w(k) = v + u·k. **Flouri T, Kobert K, Rognes T, Stamatakis A (2015)** "Are all global alignment algorithms and implementations correct?" bioRxiv 10.1101/031500 — Gotoh's P/Q border initialization is wrong (opened via search snippets; bioRxiv PDF blocked). Implemented as `GlobalAlignAffine` with M/X/Y states, X(i,0) = o + i·e, Y(0,j) = o + j·e, all other border cells −∞, state-switching traceback.
- **NCBI BLAST+ User Manual** (NBK279684, blastn options; via search snippet): blastn task default reward 2 / penalty −3 / gapopen 5 / gapextend 2 → `SequenceAligner.BlastDna`.
- Cross-check: Biopython 1.88 `PairwiseAligner` (`open_gap_score = GapOpen + GapExtend`) and parasail 1.3.4 `nw_trace` (`open = -(GapOpen + GapExtend)`) agree with `GlobalAlignAffine` on 2008 cases (random 1–40 nt, 8 scoring schemes, edge cases). A vs CAC (1, −10, −5, −1): −11 (Gotoh's original border: −9).

## Deviations and Assumptions

**None.**

- The implementation follows the standard Needleman–Wunsch linear gap penalty model exactly as given in the Wikipedia pseudocode.
- `ScoringMatrix.GapExtend` acts as the linear gap penalty $d$. `ScoringMatrix.GapOpen` is not used by `GlobalAlign`; it is used by `GlobalAlignAffine` (see above).
- When multiple optimal alignments exist, one is returned deterministically. This is explicitly allowed by the source.
- Empty-input and null-argument handling are API-level contract behaviors, not part of the NW algorithm specification.

---

## References

1. Needleman, Saul B. & Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology. 48(3): 443–53.
2. Wikipedia contributors. "Needleman–Wunsch algorithm". In: Wikipedia, The Free Encyclopedia. Accessed: 2026-03-06.
3. Wikipedia contributors. "Sequence alignment". In: Wikipedia, The Free Encyclopedia. Accessed: 2026-02-01.
