# Evidence: ALIGN-SEMI-001

**Test Unit ID:** ALIGN-SEMI-001  
**Algorithm:** Semi-Global Alignment  
**Date Collected:** 2026-02-01  

---

## 1. Authoritative Sources

| Source | URL | Accessed |
|--------|-----|----------|
| Wikipedia: Sequence alignment | https://en.wikipedia.org/wiki/Sequence_alignment#Global_and_local_alignments | 2026-02-01 |
| Wikipedia: Needleman–Wunsch algorithm | https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm | 2026-02-01 |
| Brudno et al. (2003): Glocal alignment | doi:10.1093/bioinformatics/btg1005 | 2026-02-01 |

---

## 2. Algorithm Definition (from sources)

### 2.1 Definition

Semi-global alignment (also known as "glocal" or "ends-free" alignment) is a **hybrid method** that searches for the best possible partial alignment of two sequences. It combines features of both global and local alignment, allowing alignment to:
- Start anywhere in one or both sequences (free start gaps)
- End anywhere in one or both sequences (free end gaps)

**Source:** Wikipedia (Sequence alignment, "Global and local alignments" section):
> "Hybrid methods, known as semi-global or 'glocal' (short for global-local) methods, search for the best possible partial alignment of the two sequences (in other words, a combination of one or both starts and one or both ends is stated to be aligned)."

### 2.2 Primary Use Cases

1. **Overlap alignment**: When the downstream part of one sequence overlaps with the upstream part of another sequence. Neither global nor local alignment is entirely appropriate:
   - Global alignment would force extension beyond the region of overlap
   - Local alignment might not fully cover the region of overlap

2. **Short read alignment**: When one sequence is short (e.g., a gene) and the other is long (e.g., a chromosome). The short sequence should be globally aligned, but only a local/partial alignment is desired for the long sequence.

**Source:** Wikipedia (Sequence alignment):
> "This can be especially useful when the downstream part of one sequence overlaps with the upstream part of the other sequence... Another case where semi-global alignment is useful is when one sequence is short (for example a gene sequence) and the other is very long (for example a chromosome sequence)."

### 2.3 Algorithm Variants

Multiple semi-global alignment configurations exist depending on which ends are free:

| Variant | Description | Use Case |
|---------|-------------|----------|
| Query-in-reference | Free gaps at start/end of reference (seq2) | Short read mapping, primer alignment |
| Overlap | Free gaps at end of seq1, start of seq2 | Sequence assembly |
| Full semi-global | Free gaps at all ends | General substring finding |

**Source:** ASSUMPTION based on algorithmic theory. The specific variant depends on which ends are initialized to zero.

---

## 3. Algorithm Mechanics (from sources)

### 3.1 Matrix Initialization

Semi-global alignment modifies the Needleman–Wunsch initialization:

| Alignment Type | First Row | First Column |
|----------------|-----------|--------------|
| Global (NW) | Gap penalties accumulate | Gap penalties accumulate |
| Local (SW) | All zeros | All zeros |
| Semi-global (query-in-ref) | All zeros (free end gaps in ref) | Gap penalties (query must align fully) |

**Source:** Derived from Needleman–Wunsch algorithm description (Wikipedia): "First row and first column are subject to gap penalty."

### 3.2 Traceback

For semi-global alignment (query-in-reference variant):
- Traceback starts from the **maximum score in the last row** (not just the bottom-right cell)
- Traceback proceeds to the top-left cell (ensuring full query coverage)
- Trailing gaps in the reference after alignment end are appended without penalty

**Source:** ASSUMPTION based on implementation analysis and algorithmic theory.

### 3.3 Complexity

- Time: O(n × m) where n and m are sequence lengths
- Space: O(n × m) for the scoring matrix

**Source:** Wikipedia (Needleman–Wunsch algorithm): "The time complexity of the algorithm for two sequences of length n and m is O(mn)."

---

## 4. Documented Corner Cases

| Case | Evidence Source | Expected Behavior |
|------|-----------------|-------------------|
| Query embedded in reference | Sequence alignment theory | Full query aligned, reference has leading/trailing gaps |
| Query overlaps reference end | Sequence alignment (overlap alignment) | Correct overlap detected without penalty for trailing gaps |
| Identical sequences | ASSUMPTION | Full alignment with maximum score |
| Query longer than reference | ASSUMPTION | Query may have internal gaps or poor alignment |
| Empty sequence | ASSUMPTION | Return empty result or throw |
| Null sequence | ASSUMPTION | Throw ArgumentNullException |

---

## 5. Key Invariants

| ID | Invariant | Evidence |
|----|-----------|----------|
| INV-1 | AlignmentType is SemiGlobal | Implementation contract |
| INV-2 | Query sequence (seq1) is fully represented (no clipping at ends) | Semi-global alignment definition |
| INV-3 | Aligned sequences have equal length | Alignment theory |
| INV-4 | Removing gaps from aligned seq1 yields original seq1 | Query must be fully aligned |
| INV-5 | Removing gaps from aligned seq2 yields substring of original seq2 | Reference may have free end gaps |

---

## 6. Fallback Notes

Authoritative online sources for semi-global alignment are limited compared to global and local alignment. The following aspects are marked as ASSUMPTION:
- Specific initialization behavior for the query-in-reference variant
- Traceback from maximum in last row
- Trailing gap handling

These assumptions are validated against the current implementation in Seqeron.Genomics.

---

## 7. References

1. Wikipedia. "Sequence alignment." https://en.wikipedia.org/wiki/Sequence_alignment
2. Wikipedia. "Needleman–Wunsch algorithm." https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm  
3. Brudno M, Malde S, Poliakov A, Do CB, Couronne O, Dubchak I, Batzoglou S (2003). "Glocal alignment: finding rearrangements during alignment." Bioinformatics. 19 Suppl 1:i54-62. doi:10.1093/bioinformatics/btg1005
