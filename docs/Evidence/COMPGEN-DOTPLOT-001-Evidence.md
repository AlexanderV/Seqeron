# Evidence Artifact: COMPGEN-DOTPLOT-001

**Test Unit ID:** COMPGEN-DOTPLOT-001
**Algorithm:** Dot Plot Generation (word-match / k-tuple dot matrix)
**Date Collected:** 2026-06-14

---

## Online Sources

### Gibbs & McIntyre (1970) — original dot-matrix method (primary)

**URL:** https://febs.onlinelibrary.wiley.com/doi/10.1111/j.1432-1033.1970.tb01046.x
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary paper)

**Retrieval:** Web search for `Gibbs McIntyre 1970 "The Diagram a Method for Comparing Sequences" European Journal Biochemistry DOI 16 1-11`; the Wiley landing page confirmed the bibliographic record. The full-text PDF is paywalled (HTTP 403 on direct fetch), so the method description below is taken from the secondary sources that cite this primary (Wikipedia, EMBOSS), and only the citation/DOI is attributed to this entry.

**Key Extracted Points:**

1. **Origin:** The dot-matrix ("diagram") method for comparing two sequences was introduced by Gibbs & McIntyre (1970). DOI confirmed from the Wiley record: `10.1111/j.1432-1033.1970.tb01046.x`, *Eur. J. Biochem.* **16**(1):1–11.

### Wikipedia — Dot plot (bioinformatics) (cites Gibbs & McIntyre 1970)

**URL:** https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics)
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing the primary Gibbs & McIntyre 1970)

**Retrieval:** WebFetch of the URL above.

**Key Extracted Points:**

1. **Dot-placement rule:** "When the residues of both sequences match at the same location on the plot, a dot is drawn at the corresponding position." A dot at row i / column j means the residues at positions i and j match.
2. **Diagonals = similarity:** Once plotted, dots combine to form diagonal lines; "Identical proteins will obviously have a diagonal line in the center of the matrix." Indels disrupt the diagonal; repeats add extra diagonals.
3. **Self-comparison main diagonal:** "The main diagonal represents the sequence's alignment with itself."
4. **Noise reduction via tuples:** "One way of reducing this noise is to only shade runs or 'tuples' of residues, e.g. a tuple of 3 corresponds to three residues in a row. This is effective because the probability of matching three residues in a row by chance is much lower than single-residue matches." (This is the word-size mechanism.)

### EMBOSS dottup — reference implementation manual (word-match dot plot)

**URL:** https://www.bioinformatics.nl/cgi-bin/emboss/help/dottup
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation, EMBOSS)

**Retrieval:** WebFetch of the URL above.

**Key Extracted Points:**

1. **Algorithm:** "dottup looks for places where words (tuples) of a specified length have an exact match in both sequences and draws a diagonal line over the position of these words." This is exact word matching (not substitution-matrix scoring — that is the separate `dotmatcher` tool).
2. **Word size trade-off:** "Using a longer word (tuple) size displays less random noise, runs extremely quickly, but is less sensitive. Shorter word sizes are more sensitive to shorter or fragmentary regions of similarity, but also display more random points of similarity (noise) and runs slower."
3. **Characterisation:** "a fast, but not especially sensitive way of creating dotplots" for "displaying regions of substantial similarity between two sequences."

### EMBOSS dottup manpage (default wordsize)

**URL:** https://manpages.ubuntu.com/manpages/xenial/man1/dottup.1e.html
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference implementation man page)

**Retrieval:** WebFetch of the URL above.

**Key Extracted Points:**

1. **Purpose:** "Displays a wordmatch dotplot of two sequences."
2. **Default word size:** wordsize "Default value: 10". This is the source for the implementation's `wordSize` default of 10.

### gavinhuttley — Topics in Bioinformatics, Dotplot (worked example)

**URL:** https://gavinhuttley.github.io/tib/seqcomp/dotplot.html
**Accessed:** 2026-06-14
**Authority rank:** 4 (course material citing the standard method)

**Retrieval:** WebFetch of the URL above.

**Key Extracted Points:**

1. **Match rule (k=1):** "if X[i] == Y[j] then matches[i, j] = 0" — a match (drawn black) is placed at (i, j) exactly when the characters are equal; this is "a k-mer matching algorithm where k=1".
2. **Worked example:** For X = `AGCGT` (rows) and Y = `AT` (columns), matches occur at (i=0, j=0) for the two A's and (i=4, j=1) for the two T's. (Their matrix is indexed matches[first-seq-index, second-seq-index].)
3. **Diagonal interpretation:** "Long stretches of identity form a diagonal."

---

## Documented Corner Cases and Failure Modes

### From EMBOSS dottup / Wikipedia

1. **Word longer than a sequence:** No word of length `wordSize` can be formed if a sequence is shorter than `wordSize`, so no dots are produced.
2. **Random noise with small words:** Short words produce spurious off-diagonal dots (chance matches); this is the documented sensitivity/noise trade-off, not a defect.
3. **Self-comparison:** Comparing a sequence with itself always yields the full main diagonal (every aligned position matches its own word).

---

## Test Datasets

### Dataset: Huttley worked example (k=1 word match)

**Source:** gavinhuttley TIB Dotplot page (course material citing the standard dot-matrix method).

| Parameter | Value |
|-----------|-------|
| sequence1 (x-axis) | `AGCGT` |
| sequence2 (y-axis) | `AT` |
| wordSize | 1 |
| Expected (x, y) dots | (0,0) [A=A], (4,1) [T=T] |

### Dataset: Self-comparison main diagonal

**Source:** Wikipedia — "The main diagonal represents the sequence's alignment with itself."

| Parameter | Value |
|-----------|-------|
| sequence1 = sequence2 | `ACGT` |
| wordSize | 1 |
| Expected | dots include (0,0),(1,1),(2,2),(3,3) (full main diagonal) plus any off-diagonal chance matches |

### Dataset: Exact word match, wordSize > 1

**Source:** EMBOSS dottup — exact word match of specified length; SuffixTree returns all (overlapping) occurrences.

| Parameter | Value |
|-----------|-------|
| sequence1 | `ACGTACGT` |
| sequence2 | `ACGTACGT` |
| wordSize | 4 |
| Expected | every length-4 word start x=0..4: x=0/x=4 "ACGT"→y={0,4}; x=1 "CGTA"→y=1; x=2 "GTAC"→y=2; x=3 "TACG"→y=3 ⇒ (0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4) |

---

## Assumptions

1. **ASSUMPTION: Coordinate orientation (x = sequence1, y = sequence2).** The sources fix the dot at (position-in-seq-A, position-in-seq-B) but the choice of which input is the x-axis is a presentation convention. The implementation maps sequence1→x and sequence2→y; tests assert that orientation. Changing it would transpose the plot but not its information content. Not correctness-affecting for the match set as a relation.
2. **ASSUMPTION: Case-insensitive comparison.** dottup/Gibbs do not mandate case folding; the implementation upper-cases both sequences so `a` matches `A`. Standard for nucleotide/protein dot plots where case is not biologically meaningful. Documented in the algorithm doc.

---

## Recommendations for Test Coverage

1. **MUST Test:** k=1 worked example (`AGCGT` vs `AT`) yields exactly {(0,0),(4,1)}. — Evidence: Huttley worked example.
2. **MUST Test:** exact word match with wordSize=4 (`ACGTACGT` self) yields exactly {(0,0),(0,4),(1,1),(2,2),(3,3),(4,0),(4,4)} over every word start (all overlapping occurrences). — Evidence: EMBOSS dottup exact word match.
3. **MUST Test:** self-comparison main diagonal present (every (i,i) for wordSize=1). — Evidence: Wikipedia main-diagonal statement.
4. **MUST Test:** completely dissimilar sequences (disjoint alphabets) yield no dots. — Evidence: dot-placement rule (dot only on match).
5. **MUST Test:** word longer than a sequence yields no dots; null/empty yields no dots. — Evidence: dottup word formation.
6. **MUST Test:** non-positive wordSize / stepSize throws ArgumentOutOfRangeException. — Evidence: sibling validation convention; a non-positive window is undefined for dottup.
7. **SHOULD Test:** stepSize=2 samples every other x; only even x coordinates appear. — Rationale: documented sampling parameter.
8. **SHOULD Test:** case-insensitivity (`acgt` vs `ACGT`). — Rationale: implementation normalization.
9. **COULD Test:** O(n×m) invariant — number of dots ≤ count of word-start positions × occurrences; property test on shifted-substring sequences. — Rationale: complexity bound.

---

## References

1. Gibbs AJ, McIntyre GA (1970). The Diagram, a Method for Comparing Sequences. Its Use with Amino Acid and Nucleotide Sequences. *Eur. J. Biochem.* 16(1):1–11. https://doi.org/10.1111/j.1432-1033.1970.tb01046.x
2. Rice P, Longden I, Bleasby A (2000). EMBOSS: The European Molecular Biology Open Software Suite. *Trends Genet.* 16(6):276–277 — `dottup` word-match dot plot tool. Manual: https://www.bioinformatics.nl/cgi-bin/emboss/help/dottup ; manpage (default wordsize 10): https://manpages.ubuntu.com/manpages/xenial/man1/dottup.1e.html
3. Wikipedia contributors. Dot plot (bioinformatics). https://en.wikipedia.org/wiki/Dot_plot_(bioinformatics) (accessed 2026-06-14).
4. Huttley G. Topics in Bioinformatics — Dotplot. https://gavinhuttley.github.io/tib/seqcomp/dotplot.html (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
