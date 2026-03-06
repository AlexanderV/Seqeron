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

---

## 6. References

1. Wikipedia. "Sequence alignment." https://en.wikipedia.org/wiki/Sequence_alignment
2. Wikipedia. "Needleman–Wunsch algorithm." https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm
3. Rosalind. "Finding a Motif with Modifications (SIMS)." https://rosalind.info/problems/sims/
4. Rosalind. "Semiglobal Alignment (SMGB)." https://rosalind.info/problems/smgb/
5. Brudno M, Malde S, Poliakov A, Do CB, Couronne O, Dubchak I, Batzoglou S (2003). "Glocal alignment: finding rearrangements during alignment." Bioinformatics 19 Suppl 1:i54-62. doi:10.1093/bioinformatics/btg1005
