# Evidence Artifact: PAT-APPROX-003

**Test Unit ID:** PAT-APPROX-003
**Algorithm:** Best Match and Frequency Analysis (Approximate Pattern Matching, Frequent Words with Mismatches)
**Date Collected:** 2026-06-13

---

## Online Sources

### ROSALIND BA1I — Find the Most Frequent Words with Mismatches in a String

**URL:** https://rosalind.info/problems/ba1i/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the URL in this session)
**Authority rank:** 1 (textbook problem; Compeau & Pevzner, *Bioinformatics Algorithms*, ch. 1)

**Key Extracted Points:**

1. **Definition of Count_d:** "Count_d(Text, Pattern) [is] the total number of occurrences of Pattern in Text with at most d mismatches." (verbatim from page)
2. **Worked example:** "Count1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4 because AAAAA appears four times in this string with at most one mismatch: AACAA, ATAAA, AAACA, and AAAGA." (verbatim)
3. **Pattern need not be a substring:** A most-frequent k-mer with mismatches need not occur exactly in Text (AAAAA is the most frequent 5-mer with 1 mismatch even though it does not appear exactly).
4. **Problem statement:** "Given: A string Text as well as integers k and d. Return: All most frequent k-mers with up to d mismatches in Text."
5. **Sample dataset:** Text = `ACGTTGCATGTCGCATGATGCATGAGAGCT`, k = 4, d = 1.
6. **Sample output:** `GATG ATGC ATGT` (set of most-frequent 4-mers; order not significant). Verified max count = 5 by running the repository implementation on this dataset in-session.
7. **Practical bound:** "the solution should work for k ≤ 12 and d ≤ 3" (the neighbor-enumeration approach becomes slow for larger k and d).

### ROSALIND BA1H — Find All Approximate Occurrences of a Pattern in a String

**URL:** https://rosalind.info/problems/ba1h/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the URL in this session)
**Authority rank:** 1 (textbook problem; Compeau & Pevzner, *Bioinformatics Algorithms*, ch. 1)

**Key Extracted Points:**

1. **Definition of approximate occurrence:** "A k-mer Pattern appears as a substring of Text with at most d mismatches if there is some k-mer substring Pattern' of Text having d or fewer mismatches with Pattern, i.e., HammingDistance(Pattern, Pattern') ≤ d." (verbatim)
2. **Problem statement:** "Given: Strings Pattern and Text along with an integer d. Return: All starting positions where Pattern appears as a substring of Text with at most d mismatches."
3. **Sample dataset:** Pattern = `ATTCTGGA`, Text = `CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC`, d = 3.
4. **Sample output (positions, 0-based):** `6 7 26 27 78`. (`Count_d` for this dataset therefore equals 5.)

### ROSALIND BA1N — Generate the d-Neighborhood of a String

**URL:** https://rosalind.info/problems/ba1n/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the URL in this session)
**Authority rank:** 1 (textbook problem; Compeau & Pevzner, *Bioinformatics Algorithms*, ch. 1)

**Key Extracted Points:**

1. **d-neighborhood definition:** "the set of all k-mers whose Hamming distance from Pattern does not exceed d." (verbatim)
2. **Sample dataset:** Pattern = `ACG`, d = 1.
3. **Sample output (10 neighbors):** `CCG TCG GCG AAG ATG AGG ACA ACC ACT ACG`. The neighborhood includes the pattern itself (ACG) and has size 1 + 3·k = 1 + 9 = 10 for k = 3, d = 1 (each of the 3 positions can take 3 alternative bases, plus the identity). This confirms the recursive structure used by FrequentWordsWithMismatches.

### Reference implementation — charlesreid1/go-rosalind (`rosalind/rosalind_ba1.go`)

**URL:** https://raw.githubusercontent.com/charlesreid1/go-rosalind/master/rosalind/rosalind_ba1.go
**Accessed:** 2026-06-13 (retrieved via WebFetch of the raw file in this session)
**Authority rank:** 3 (established reference implementation of the textbook problems)

**Key Extracted Points:**

1. **Approximate matching (BA1H):** for each k-length window of Text, compute `HammingDistance(window, pattern)`; record the start index where `hamm <= d`. O(n·m).
2. **FrequentWordsWithMismatches (BA1I):** for each k-mer window, precompute its Hamming neighbors, then iterate windows again incrementing the count "for kmer and neighbors (note that neighbors includes the kmer itself)"; the most-frequent neighbor(s) are returned. This is exactly the tally-over-neighbors approach used by the repository.
3. **Neighbors recursion (BA1N):** recursive position-choice enumeration; the base case substitutes non-matching bases at the chosen indices — equivalent to the textbook `Neighbors(Pattern, d)`.

### Reference implementation — zonghui0228/Rosalind-Solutions (`rosalind_ba1h.py`)

**URL:** https://github.com/zonghui0228/Rosalind-Solutions/blob/master/code/rosalind_ba1h.py
**Accessed:** 2026-06-13 (retrieved via WebFetch in this session)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Confirms O(n·m):** iterate each possible substring position, extract a pattern-length substring, compute Hamming distance, record positions where distance ≤ d. "time complexity is O(n·m) where n is the text length and m is the pattern length".

---

## Documented Corner Cases and Failure Modes

### From ROSALIND BA1I / BA1H

1. **Pattern absent as exact substring:** the most-frequent mismatch k-mer (BA1I) and the matched pattern (BA1H) need not occur exactly in Text; counting is over the d-neighborhood / Hamming ball, not exact occurrences.
2. **Multiple ties (BA1I):** when several k-mers share the maximum count, ALL of them must be returned (sample output has three: GATG, ATGC, ATGT).
3. **d = 0 degenerates to exact:** Count_0 is exact occurrence counting; Neighbors(Pattern, 0) = {Pattern}.

### From the go-rosalind reference

1. **Neighborhood includes identity:** Neighbors(Pattern, d) always contains Pattern itself; a window's own k-mer is therefore always counted.

---

## Test Datasets

### Dataset: BA1I Sample (Frequent Words with Mismatches)

**Source:** ROSALIND BA1I (https://rosalind.info/problems/ba1i/)

| Parameter | Value |
|-----------|-------|
| Text | `ACGTTGCATGTCGCATGATGCATGAGAGCT` |
| k | 4 |
| d | 1 |
| Expected output (set) | {GATG, ATGC, ATGT} |
| Max count | 5 |

### Dataset: BA1H Sample (Approximate Pattern Matching)

**Source:** ROSALIND BA1H (https://rosalind.info/problems/ba1h/)

| Parameter | Value |
|-----------|-------|
| Pattern | `ATTCTGGA` |
| Text | `CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC` |
| d | 3 |
| Expected positions (0-based) | 6, 7, 26, 27, 78 |
| Count_d | 5 |

### Dataset: Count_1 Worked Example

**Source:** ROSALIND BA1I worked example

| Parameter | Value |
|-----------|-------|
| Text | `AACAAGCTGATAAACATTTAAAGAG` |
| Pattern | `AAAAA` |
| d | 1 |
| Count_1 | 4 (windows AACAA, ATAAA, AAACA, AAAGA) |

### Dataset: BA1N d-Neighborhood

**Source:** ROSALIND BA1N (https://rosalind.info/problems/ba1n/)

| Parameter | Value |
|-----------|-------|
| Pattern | `ACG` |
| d | 1 |
| Neighborhood size | 10 |
| Members | CCG, TCG, GCG, AAG, ATG, AGG, ACA, ACC, ACT, ACG |

---

## Assumptions

1. **ASSUMPTION: FindBestMatch tie-breaking (first minimum-distance window).** No Rosalind/textbook problem defines a single "best approximate match" return. The repository method `FindBestMatch(sequence, pattern)` returns the equal-length window with minimum Hamming distance, scanning left-to-right and keeping the first window that achieves a strictly smaller distance (so the leftmost minimum-distance window is returned; an exact match short-circuits). The *distance value* and *which windows are minimal* are fully evidence-defined by HammingDistance (BA1H/PAT-APPROX-001); only the leftmost tie-break choice is a documented API convention, not a correctness-affecting algorithm parameter (it does not change the returned minimum distance). Tested as a deterministic, documented convention.

---

## Recommendations for Test Coverage

1. **MUST Test:** `FindFrequentKmersWithMismatches` on BA1I sample returns exactly {GATG, ATGC, ATGT} with count 5 — Evidence: ROSALIND BA1I sample.
2. **MUST Test:** `CountApproximateOccurrences(Text, ATTCTGGA, 3)` on BA1H sample returns 5; the underlying match positions are {6,7,26,27,78} — Evidence: ROSALIND BA1H sample.
3. **MUST Test:** `CountApproximateOccurrences(AACAAGCTGATAAACATTTAAAGAG, AAAAA, 1)` returns 4 — Evidence: ROSALIND BA1I Count_1 worked example.
4. **MUST Test:** `FindBestMatch` returns distance 0 on an exact match (IsExact). — Evidence: HammingDistance identity (BA1H definition).
5. **SHOULD Test:** d = 0 reduces FrequentWords to exact frequent k-mer and Count to exact occurrence count — Rationale: documented degenerate case.
6. **SHOULD Test:** FindBestMatch returns the leftmost minimum-distance window when no exact match exists — Rationale: documented tie-break convention.
7. **COULD Test:** null/empty/too-short inputs and invalid k/d argument exceptions — Rationale: contract robustness.

---

## References

1. Compeau, P., Pevzner, P. (2015). *Bioinformatics Algorithms: An Active Learning Approach*, Ch. 1 ("Where in the Genome Does DNA Replication Begin?"). ROSALIND textbook track problems BA1H/BA1I/BA1N. https://rosalind.info/problems/ba1h/ , https://rosalind.info/problems/ba1i/ , https://rosalind.info/problems/ba1n/
2. charlesreid1. go-rosalind reference implementation, `rosalind/rosalind_ba1.go`. https://raw.githubusercontent.com/charlesreid1/go-rosalind/master/rosalind/rosalind_ba1.go
3. zonghui0228. Rosalind-Solutions, `code/rosalind_ba1h.py`. https://github.com/zonghui0228/Rosalind-Solutions/blob/master/code/rosalind_ba1h.py

---

## Change History

- **2026-06-13**: Initial documentation (PAT-APPROX-003).
