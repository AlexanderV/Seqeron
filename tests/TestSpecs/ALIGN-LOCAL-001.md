# Test Specification: ALIGN-LOCAL-001

**Test Unit ID:** ALIGN-LOCAL-001
**Area:** Alignment
**Algorithm:** Local Alignment (Smith–Waterman)
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-03-07

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Smith–Waterman algorithm | https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm | 2026-03-06 |

### 1.2 Key Evidence Points

1. Smith–Waterman performs local sequence alignment, finding the highest-scoring pair of subsequences between two sequences. (Smith–Waterman algorithm)
2. The key difference from Needleman–Wunsch is the **zero floor**: negative scores are set to 0, allowing alignment to restart at any position. (Smith–Waterman algorithm, §Comparison)
3. Initialization: first row and first column are set to 0 (no end-gap penalty). (Smith–Waterman algorithm, §Scoring matrix)
4. Traceback: begins at the cell with the **highest score** in the matrix and ends when a cell with score 0 is encountered. (Smith–Waterman algorithm, §Algorithm step 4)
5. Linear gap penalty recurrence: $H_{ij} = \max(0,\; H_{i-1,j-1} + s(a_i,b_j),\; H_{i-1,j} - W_1,\; H_{i,j-1} - W_1)$ (Smith–Waterman algorithm, §Linear gap penalty)
6. Wikipedia example: sequences `TGTTACGG` and `GGTTGACTA` with match +3, mismatch −3, linear gap penalty $W_1 = 2$ yields alignment `GTT-AC` / `GTTGAC` with score 13. (Smith–Waterman algorithm, §Example)

### 1.3 Documented Corner Cases

| Case | Evidence Source | Notes |
|------|-----------------|-------|
| No similarity between sequences | Smith–Waterman zero floor (§Algorithm) | When all comparisons yield negative scores, every $H_{ij} = \max(0, \text{neg}, \text{neg}, \text{neg}) = 0$. Score = 0, alignment is empty. |
| Identical sequences | Derivable from S-W recurrence | For identical sequences of length $n$ with match score $M > 0$ and $M > \|W_1\|$: diagonal dominates at every cell, $H_{i,i} = i \times M$. Maximum $= n \times M$, full alignment, no gaps. |
| Empty sequence input | Smith–Waterman definition (§Scoring matrix) | Matrix dimensions are $1 \times (m+1)$ or $(n+1) \times 1$. Only initialized row/column (all zeros) — no cells to fill. Max score = 0. |

### 1.4 Known Failure Modes / Pitfalls (from sources)

1. Multiple optimal local alignments can exist when several cells share the same maximum score. Implementation returns one. (Smith–Waterman algorithm, §Example)
2. Alignment quality depends on the scoring system and gap penalties. (Smith–Waterman algorithm, §Explanation)
3. Linear gap penalty may not model biological indels accurately; affine gaps (Gotoh) are often preferred. (Smith–Waterman algorithm, §Affine gap penalty)

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `LocalAlign(DnaSequence, DnaSequence, ScoringMatrix?)` | SequenceAligner | **Canonical** | Smith–Waterman local alignment |
| `LocalAlign(string, string, ScoringMatrix?)` | SequenceAligner | Delegate | String wrapper — delegates to same `LocalAlignCore` after `ToUpperInvariant()` |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | Score ≥ 0 (zero floor property) | Yes | Smith–Waterman algorithm: $\max(0, \ldots)$ in recurrence |
| INV-3 | AlignmentType is Local | Yes | Implementation contract |
| INV-4 | Removing gaps from aligned sequences yields substrings of originals | Yes | Smith–Waterman traceback rules |
| INV-5 | StartPosition and EndPosition are within sequence bounds | Yes | Smith–Waterman traceback from matrix positions |

---

## 4. Test Cases

### 4.1 MUST Tests (Required for DoD)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| M1 | Wikipedia example: exact score and alignment (gap in seq1, "left" traceback branch) | seq1=`TGTTACGG`, seq2=`GGTTGACTA`, scoring: match +3, mismatch −3, gap −2 | Score=13, aligned1=`GTT-AC`, aligned2=`GTTGAC`, start1=1, end1=5, start2=1, end2=6, INV-1,3,4,5 hold | Smith–Waterman §Example, hand-verified scoring matrix |
| M2 | Swapped Wikipedia example: gap in seq2, "up" traceback branch | seq1=`GGTTGACTA`, seq2=`TGTTACGG`, same scoring | Score=13, aligned1=`GTTGAC`, aligned2=`GTT-AC`, start1=1, end1=6, start2=1, end2=5, INV-3,4,5 hold | Derived from Wikipedia via S-W transposition: $H^{AB}_{i,j} = H^{BA}_{j,i}$; hand-verified transposed matrix |
| M3 | String overload produces same Wikipedia-expected result | Same Wikipedia sequences as string | Score=13, aligned1=`GTT-AC`, aligned2=`GTTGAC` | Implementation contract — delegates to `LocalAlignCore`; validated against Wikipedia expected values |
| M4 | Empty string input returns `AlignmentResult.Empty` | seq1="", seq2="ACGT" | Empty result | S-W definition: 0-dimensional scoring matrix → maxScore=0 |
| M5 | Null DnaSequence throws `ArgumentNullException` | seq1=null, seq2=valid | ArgumentNullException | .NET convention: `ArgumentNullException.ThrowIfNull` |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Input | Expected | Evidence |
|----|-----------|-------|----------|----------|
| S1 | Identical sequences produce full-length alignment with exact score | seq1=seq2=`ACGTACGT`, scoring: match +3, mismatch −3, gap −2 | Score=24 (8×3), full alignment `ACGTACGT`/`ACGTACGT` | Derivable: identical seqs → diagonal dominates → score = n×match |
| S2 | Completely dissimilar sequences produce score exactly 0 | seq1=`AAAA`, seq2=`TTTT`, scoring: match +3, mismatch −3, gap −2 | Score=0, empty alignment | Derivable: no matches → all $H_{ij} = \max(0, \text{neg}, \text{neg}, \text{neg}) = 0$ |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Input | Expected | Notes |
|----|-----------|-------|----------|-------|
| C1 | Multiple optimal alignment detection | Sequences with tied max scores | Any optimal alignment is acceptable | Evidence: multiple optimal alignments can exist (§Example) |

---

## 5. Coverage Classification

| Area | Status | Tested by | Notes |
|------|--------|-----------|-------|
| Wikipedia example exact values (score, alignment, positions) | ✅ | M1 | Score=13, `GTT-AC`/`GTTGAC`, positions (1,5)/(1,6) |
| Traceback "left" branch (gap in seq1) | ✅ | M1 | `GTT-AC` has gap at position 4 in aligned seq1 |
| Traceback "up" branch (gap in seq2) | ✅ | M2 | `GTT-AC` has gap at position 4 in aligned seq2 |
| Traceback "diagonal" branch (match/mismatch) | ✅ | M1, M2 | Both tests traverse multiple diagonal steps |
| S-W transposition symmetry | ✅ | M1+M2 | Same score=13 for swapped sequences |
| Zero floor (INV-1) | ✅ | M1, S2 | Positive (13) and zero (0) boundary values |
| AlignmentType is Local (INV-3) | ✅ | M1, M2 | Checked in both dataset tests |
| Gaps removal → substring (INV-4) | ✅ | M1, M2 | Exact gap-removed strings + substring containment |
| Position bounds (INV-5) | ✅ | M1, M2 | Exact positions from hand-verified traceback |
| String overload parity | ✅ | M3 | Validated against Wikipedia expected values |
| Empty input handling | ✅ | M4 | Returns `AlignmentResult.Empty` |
| Null argument handling | ✅ | M5 | `ArgumentNullException` |
| Identical sequences (all-match boundary) | ✅ | S1 | Score=24 (8×3), full alignment |
| Dissimilar sequences (all-mismatch boundary) | ✅ | S2 | Score=0, empty alignment |

---

## 6. Open Questions / Decisions

None.

---

## 7. Hand-Verified Scoring Matrices

### 7.1 Wikipedia Example (M1)

Sequences: A = `TGTTACGG`, B = `GGTTGACTA`
Scoring: match = +3, mismatch = −3, linear gap $W_1 = 2$

```
       -   G   G   T   T   G   A   C   T   A
  -  [ 0   0   0   0   0   0   0   0   0   0 ]
  T  [ 0   0   0   3   3   1   0   0   3   1 ]
  G  [ 0   3   3   1   1   6   4   2   1   0 ]
  T  [ 0   1   1   6   4   4   3   1   5   3 ]
  T  [ 0   0   0   4   9   7   5   3   4   2 ]
  A  [ 0   0   0   2   7   6  10   8   6   7 ]
  C  [ 0   0   0   0   5   4   8  13  11   9 ]
  G  [ 0   3   3   1   3   8   6  11  10   8 ]
  G  [ 0   3   6   4   2   6   5   9   8   7 ]
```

**Maximum**: 13 at H[6,7] (row=C, col=C)

**Traceback**: H[6,7]=13 →(diag) H[5,6]=10 →(diag) H[4,5]=7 →(**left**) H[4,4]=9 →(diag) H[3,3]=6 →(diag) H[2,2]=3 →(diag) H[1,1]=0 STOP

**Result**: `GTT-AC` / `GTTGAC`, score = 13

### 7.2 Swapped Wikipedia Example (M2)

Sequences: A = `GGTTGACTA`, B = `TGTTACGG`
Derived from §7.1 via S-W transposition: $H^T[i,j] = H[j,i]$

```
       -   T   G   T   T   A   C   G   G
  -  [ 0   0   0   0   0   0   0   0   0 ]
  G  [ 0   0   3   1   0   0   0   3   3 ]
  G  [ 0   0   3   1   0   0   0   3   6 ]
  T  [ 0   3   1   6   4   2   0   1   4 ]
  T  [ 0   3   1   4   9   7   5   3   2 ]
  G  [ 0   1   6   4   7   6   4   8   6 ]
  A  [ 0   0   4   3   5  10   8   6   5 ]
  C  [ 0   0   2   1   3   8  13  11   9 ]
  T  [ 0   3   1   5   4   6  11  10   8 ]
  A  [ 0   1   0   3   2   7   9   8   7 ]
```

**Maximum**: 13 at H^T[7,6] (row=C, col=C)

**Traceback**: H^T[7,6]=13 →(diag) H^T[6,5]=10 →(diag) H^T[5,4]=7 →(**up**) H^T[4,4]=9 →(diag) H^T[3,3]=6 →(diag) H^T[2,2]=3 →(diag) H^T[1,1]=0 STOP

**Result**: `GTTGAC` / `GTT-AC`, score = 13
