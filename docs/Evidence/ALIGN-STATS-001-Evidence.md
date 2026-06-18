# Evidence Artifact: ALIGN-STATS-001

**Test Unit ID:** ALIGN-STATS-001
**Algorithm:** Pairwise Alignment Statistics (Identity / Similarity / Gaps) and Alignment Formatting
**Date Collected:** 2026-06-13

---

## Online Sources

### EMBOSS needle — application documentation (release 6.6)

**URL:** https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation — EMBOSS, the standard global aligner)

**Key Extracted Points:**

1. **Identity definition:** "the percentage of identical matches between the two sequences over the reported aligned region (including any gaps in the length)." Identity counts only exact character matches; the denominator is the full alignment length **including gap columns**.
2. **Similarity definition:** "the percentage of matches between the two sequences over the reported aligned region (including any gaps in the length)." Similarity is broader than Identity: it includes positions that are not identical but score positively under the substitution matrix. In the worked example 65 positions are identical but 90 are "matching" — i.e. 25 non-identical positions score positively and are counted as similar.
3. **Gaps definition:** the count and proportion of gap positions within the aligned sequences; denominator is the alignment length.
4. **Score definition:** "the sum of the matches taken from the scoring matrix, minus penalties arising from opening and extending gaps."
5. **Worked example (HBA_HUMAN vs HBB_HUMAN, EBLOSUM62):** verbatim output block —
   ```
   # Length: 149
   # Identity:      65/149 (43.6%)
   # Similarity:    90/149 (60.4%)
   # Gaps:           9/149 ( 6.0%)
   # Score: 292.5
   ```
   Confirms: percentage = count / Length × 100, denominator = Length (149) which includes the 9 gap columns. 65/149 = 43.6%, 90/149 = 60.4%, 9/149 = 6.0%.

### EMBOSS — Alignment Formats (AlignFormats / srspair markup legend)

**URL:** https://emboss.sourceforge.net/docs/themes/AlignFormats.html
**Accessed:** 2026-06-13
**Authority rank:** 3 (EMBOSS reference documentation)

**Key Extracted Points:**

1. **Standard (pair / srspair) markup line legend:** `|` marks an identical position; `:` marks a similar position (substitution score > 1.0); `.` marks a small positive score; a space marks a mismatch or a gap.
2. **Three-line layout:** the markup/consensus line is placed between the two aligned sequence lines.

### NCBI BLAST — pairwise alignment metric definitions

**URL:** (search) "BLAST pairwise alignment identities positives gaps definition NCBI" — corroborating definitions from NCBI/UC-Berkeley BLAST guide and NCBI BLAST QuickStart (https://www.ncbi.nlm.nih.gov/books/NBK1734/)
**Accessed:** 2026-06-13
**Authority rank:** 3 (NCBI BLAST reference output)

**Key Extracted Points:**

1. **Identities:** percentage of residues with a direct (exact) match in the alignment.
2. **Positives (= Similarity):** number of residues that are either identical or have similar chemical properties — operationally, residues that **score positively in the substitution matrix**. Confirms the "positive substitution score ⇒ similar" rule independently of EMBOSS.
3. **Gaps:** dashes (`-`) indicate gaps — positions where one sequence has no counterpart.

### Percent-identity denominator convention (corroboration)

**URL:** (search) "percent identity alignment number identical positions divided by alignment length including gaps definition"
**Accessed:** 2026-06-13
**Authority rank:** 4 (secondary, corroborates the EMBOSS primary definition)

**Key Extracted Points:**

1. **Denominator includes gaps:** "the percent identity is calculated by counting the number of identical aligning residues, dividing by the total length of the aligned region, including gaps in both sequences, and multiplying by 100." Independently confirms that the alignment length (with gap columns) is the denominator.

### pseqsid — reference implementation of sequence identity/similarity

**URL:** https://github.com/amaurypm/pseqsid
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Identity:** "percentage of identical residues in equivalent positions with respect to the sequence length." With the `alignment` length option the denominator is "total alignment length including gap-only columns."
2. **Similarity:** "percentage of identical or similar residues in equivalent positions." A position counts as similar when the residues score positively under the substitution matrix; gap-only columns are ignored from the similar count.

---

## Documented Corner Cases and Failure Modes

### From EMBOSS needle / AlignFormats

1. **Gap columns counted in the denominator:** the percentage denominator is always the full alignment Length, so identity + (similar-only) + (mismatch-only) + gap fractions partition 100%.
2. **Similarity ≥ Identity always:** the similar set is a superset of the identical set (identical columns always score positively), so Similarity% ≥ Identity% for any scoring model where identical residues score positively.
3. **DNA simple match/mismatch model:** when the substitution score for a non-identical pair is non-positive (Mismatch ≤ 0, as in all DNA models exposed by this class), no non-identical column is similar, so Similarity equals Identity.

### From BLAST

1. **Empty / no-counterpart positions:** dashes represent gaps; gap columns are not identities and not positives.

---

## Test Datasets

### Dataset: EMBOSS needle HBA_HUMAN vs HBB_HUMAN (published worked example)

**Source:** EMBOSS needle documentation, release 6.6 (Rice, Longden & Bleasby 2000).

| Parameter | Value |
|-----------|-------|
| Alignment Length | 149 |
| Identity (count / %) | 65 / 149 = 43.6% |
| Similarity (count / %) | 90 / 149 = 60.4% |
| Gaps (count / %) | 9 / 149 = 6.0% |
| Score | 292.5 |

> This protein example validates the **formula and denominator** (count / Length × 100, Length includes gaps). The repository implementation operates on DNA with a simple match/mismatch model and does not ship the EBLOSUM62 matrix, so the exact 65/90/9 protein counts are reproduced as a formula-level cross-check (see Dataset below), not as a literal API call.

### Dataset: Formula reconstruction from the EMBOSS counts

**Source:** derived directly from the EMBOSS counts above; arithmetic only.

| Quantity | Computation | Expected |
|----------|-------------|----------|
| Identity % | 65 / 149 × 100 | 43.624161… |
| Similarity % | 90 / 149 × 100 | 60.402684… |
| Gap % | 9 / 149 × 100 | 6.040268… |

### Dataset: Hand-constructed DNA alignment (deterministic, exact)

**Source:** constructed for this unit; values derived by direct application of the EMBOSS formulas to a known column layout (no implementation output used).

Aligned columns (SimpleDna model: Match +1, Mismatch −1):

```
seq1: A C G T - A C G T
seq2: A C C T A A C - T
```

| Column | c1 | c2 | Class |
|--------|----|----|-------|
| 1 | A | A | identity |
| 2 | C | C | identity |
| 3 | G | C | mismatch (Mismatch −1 ⇒ not similar) |
| 4 | T | T | identity |
| 5 | - | A | gap |
| 6 | A | A | identity |
| 7 | C | C | identity |
| 8 | G | - | gap |
| 9 | T | T | identity |

Length = 9, Matches = 6, Mismatches = 1, Gaps = 2.
Identity = 6/9 × 100 = 66.6̅%; Similarity = 6/9 × 100 = 66.6̅% (no positive non-identical column under SimpleDna); Gap% = 2/9 × 100 = 22.2̅%.
Markup line (srspair legend, SimpleDna): `||` `space` `|` `space` `||` `space` `|` → `||.|` is **not** used; mismatch and gap are both space: `|| | || |` with column 3 a space and columns 5,8 spaces.

### Dataset: Similarity > Identity case (positive-scoring mismatch)

**Source:** constructed; uses a scoring model with Mismatch = +1 so that every non-identical column scores positively and is counted as similar (direct application of the "positive substitution score ⇒ similar" rule from EMBOSS/BLAST/pseqsid).

```
seq1: A C G T
seq2: A C C T   (column 3 non-identical, Mismatch = +1 ⇒ similar)
```

With `ScoringMatrix(Match: 1, Mismatch: 1, GapOpen: 0, GapExtend: -1)`:
Length = 4, Matches = 3, Mismatches = 1, Gaps = 0.
Identity = 3/4 × 100 = 75%; Similarity = (3+1)/4 × 100 = 100%; Gap% = 0%.
Markup line: `||:|` ( `:` marks the similar column).

---

## Assumptions

1. **ASSUMPTION: srspair `.` (small positive score) collapses to `:` for these models.** EMBOSS distinguishes `:` (score > 1.0) from `.` (small positive score) using the magnitude of the substitution score. The DNA models exposed by this class use a single integer Match/Mismatch scalar with no graded positive scores, so only two outcomes exist: identical (`|`) or non-positive (space); a positive non-identical score, when one is configured, is rendered `:`. The `.` tier is unreachable with the available scoring model and is therefore not emitted. This is a rendering-only choice; it does not affect any counted statistic. Resolved as non-correctness-affecting for the statistics; documented for the formatter.

---

## Recommendations for Test Coverage

1. **MUST Test:** Identity/Similarity/Gaps formula and denominator (count / Length × 100, Length includes gaps) reproduced from the EMBOSS 149-column counts (43.6% / 60.4% / 6.0%). — Evidence: EMBOSS needle worked example.
2. **MUST Test:** DNA SimpleDna alignment ⇒ Similarity equals Identity (no positive non-identical columns). — Evidence: EMBOSS/BLAST "positive substitution score ⇒ similar"; DNA Mismatch < 0.
3. **MUST Test:** Positive-scoring mismatch (Mismatch = +1) ⇒ Similarity > Identity, similar column counted. — Evidence: EMBOSS/BLAST/pseqsid similarity rule.
4. **MUST Test:** Gap counting — column with either side `-` counts as gap, not identity/mismatch. — Evidence: BLAST gap definition; EMBOSS gap count.
5. **MUST Test:** FormatAlignment markup legend `|`/`:`/space per srspair. — Evidence: EMBOSS AlignFormats.
6. **SHOULD Test:** Empty alignment ⇒ AlignmentStatistics.Empty / "" — Rationale: documented empty handling, denominator undefined.
7. **SHOULD Test:** Null alignment ⇒ ArgumentNullException; non-positive lineWidth ⇒ ArgumentOutOfRangeException — Rationale: input validation contract.
8. **COULD Test:** FormatAlignment line wrapping at lineWidth (multi-block output) — Rationale: formatting completeness.
9. **COULD Test:** Invariant Identity ≤ Similarity ≤ 100 and fractions partition to 100% (property test) — Rationale: O(n) but cheap structural invariant.

---

## References

1. Rice P, Longden I, Bleasby A (2000). EMBOSS: The European Molecular Biology Open Software Suite. *Trends in Genetics* 16(6):276–277. https://doi.org/10.1016/S0168-9525(00)02024-2
2. EMBOSS needle application documentation (release 6.6). https://emboss.sourceforge.net/apps/release/6.6/emboss/apps/needle.html (accessed 2026-06-13)
3. EMBOSS Alignment Formats (markup legend). https://emboss.sourceforge.net/docs/themes/AlignFormats.html (accessed 2026-06-13)
4. NCBI. BLAST QuickStart — Comparative Genomics. NCBI Bookshelf NBK1734. https://www.ncbi.nlm.nih.gov/books/NBK1734/ (accessed 2026-06-13)
5. pseqsid — pairwise sequence identity/similarity reference implementation. https://github.com/amaurypm/pseqsid (accessed 2026-06-13)

---

## Change History

- **2026-06-13**: Initial documentation.
