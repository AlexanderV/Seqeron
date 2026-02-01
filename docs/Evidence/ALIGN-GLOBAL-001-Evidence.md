# Evidence Artifact: ALIGN-GLOBAL-001

**Test Unit ID:** ALIGN-GLOBAL-001  
**Algorithm:** Global Alignment (Needleman–Wunsch)  
**Date Collected:** 2026-02-01  

---

## Online Sources

### Wikipedia: Needleman–Wunsch Algorithm
**URL:** https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm

**Key Extracted Points:**

1. **Purpose:**
   - Algorithm for computing an optimal global alignment between two sequences.
   - One of the first applications of dynamic programming to biological sequences.
   - Aligns the entire sequences end-to-end.

2. **Scoring System:**
   - Match: +1 (or other positive value depending on substitution matrix)
   - Mismatch: −1 (or negative value from substitution matrix)
   - Indel (gap): −1 (or specified gap penalty)
   - Recurrence: $F_{i,j} = \max(F_{i-1,j-1} + S(A_i,B_j), F_{i-1,j} + d, F_{i,j-1} + d)$
   - Initialization: $F_{0,j} = d \cdot j$, $F_{i,0} = d \cdot i$ (linear gap cost)

3. **Traceback Rules:**
   - Diagonal path (top-left cell) → match or mismatch.
   - Horizontal path (left cell) → insert gap in first sequence.
   - Vertical path (top cell) → insert gap in second sequence.
   - Multiple traceback paths from the same position → multiple optimal alignments exist.

4. **Complexity:**
   - Time: $O(nm)$ for sequences of length $n$ and $m$.
   - Space: $O(nm)$ for the full matrix.
   - Hirschberg's algorithm variant: $O(nm)$ time, $\Theta(\min(n,m))$ space.

5. **Example:**
   - Sequences: GCATGCG and GATTACA
   - Simple scoring: match = +1, mismatch/indel = −1
   - Optimal score shown in Wikipedia figures.

### Wikipedia: Sequence Alignment (Global vs. Local)
**URL:** https://en.wikipedia.org/wiki/Sequence_alignment

**Key Extracted Points:**

1. **Global Alignment:**
   - Attempts to align the entire length of all sequences (every residue).
   - Most useful when sequences are similar and of roughly equal size.
   - Can start and/or end in gaps.
   - Needleman–Wunsch is the canonical global alignment method.

2. **Local Alignment:**
   - Identifies regions of similarity within long sequences.
   - Smith–Waterman is the canonical local alignment method.
   - More useful for dissimilar sequences.

3. **Semi-Global / Glocal Alignment:**
   - Hybrid method combining one or both starts and one or both ends.
   - Useful when downstream part of one sequence overlaps upstream part of the other.
   - Useful when one sequence is short and the other is very long.

4. **Gap Penalties:**
   - Two-parameter gap cost: gap-open penalty (e.g., −5) and gap-extend penalty (e.g., −1).
   - Biologically motivated: large gaps are more likely as single deletions than multiple small ones.

5. **Multiple Optimal Alignments:**
   - Multiple traceback paths can exist for the same optimal score.
   - Common in sequences with repeated motifs or low complexity regions.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia: Needleman–Wunsch Algorithm

1. **Multiple Optimal Alignments:**
   - When multiple paths from the bottom-right cell to the top-left cell share the same maximum score, each represents an equally valid optimal alignment.
   - The specific alignment returned depends on traceback implementation choices (e.g., priority of diagonal vs. horizontal vs. vertical moves).

2. **Scoring System Dependence:**
   - Alignment quality and the resulting alignment depend critically on the scoring system chosen.
   - Different scoring systems (e.g., BLOSUM, PAM, substitution matrices) yield different alignments for the same sequences.

3. **Gap Penalty Impact:**
   - Large gap penalties favor fewer, longer gaps; small penalties favor many short gaps.
   - The choice of gap-open and gap-extend penalties affects the alignment topology.

### From Wikipedia: Sequence Alignment (General)

1. **Edge Case: Empty Sequences**
   - No explicit behavior specified in the Needleman–Wunsch sources for empty inputs.
   - Typical implementations either throw an error or return an empty/trivial result.

2. **Edge Case: Identical Sequences**
   - Should align perfectly with score = sum of match scores for all positions.

3. **Edge Case: Completely Different Sequences**
   - All mismatches and gaps; score is negative if mismatch/gap penalties are negative.

---

## Test Dataset: Wikipedia Example

| Parameter | Value |
|-----------|-------|
| Sequence 1 | GCATGCG |
| Sequence 2 | GATTACA |
| Match Score | +1 |
| Mismatch Score | −1 |
| Gap Score | −1 |
| Aligned Seq1 (one possible alignment) | GCATGCG or similar (with or without gaps depending on traceback) |
| Aligned Seq2 (one possible alignment) | GATTACA or similar |
| Invariant 1 | Aligned sequences must have equal length |
| Invariant 2 | Removing gaps from aligned sequences yields original sequences |
| Invariant 3 | Alignment score = sum of per-position match/mismatch/gap scores |

---

## Assumptions and Limitations

1. **ASSUMPTION: Empty-Input Behavior**
   - Needleman–Wunsch sources do not specify behavior for empty sequences.
   - Implementation behavior (e.g., return empty result) is assumed acceptable and is not evidence-grounded.

2. **ASSUMPTION: Affine Gap Model**
   - The implementation uses `GapOpen` and `GapExtend` parameters, suggesting a simplified affine gap model.
   - The full Gotoh algorithm (three-matrix affine) is not explicitly referenced in the cited sources; implementation-specific details are not validated against a source.

3. **ASSUMPTION: Single Traceback Path**
   - When multiple optimal alignments exist, the implementation returns a single one.
   - The sources acknowledge multiple optimal alignments exist but do not specify which one should be returned.

---

## Recommendations for Test Coverage

1. **MUST Test:** Reconstruct and verify the Wikipedia GCATGCG/GATTACA example.
2. **MUST Test:** Verify invariants for any valid global alignment (aligned sequence length, gap removal, score recalculation).
3. **MUST Test:** Verify string and DnaSequence overloads produce consistent results.
4. **SHOULD Test:** Verify null arguments throw ArgumentNullException (API contract).
5. **COULD Test:** Verify alternative optimal alignment paths when they exist (requires specific dataset design).

---

## References

1. Needleman, Saul B. & Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology. 48(3): 443–53.
2. Wikipedia contributors. "Needleman–Wunsch algorithm". In: Wikipedia, The Free Encyclopedia. Accessed: 2026-02-01.
3. Wikipedia contributors. "Sequence alignment". In: Wikipedia, The Free Encyclopedia. Accessed: 2026-02-01.
