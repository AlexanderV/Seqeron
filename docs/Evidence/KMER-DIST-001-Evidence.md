# Evidence Artifact: KMER-DIST-001

**Test Unit ID:** KMER-DIST-001
**Algorithm:** K-mer Euclidean Distance (alignment-free word-frequency distance)
**Date Collected:** 2026-06-13

---

## Online Sources

### Zielezinski, Vinga, Almeida & Karlowski (2017) — "Alignment-free sequence comparison: benefits, applications, and tools" (Genome Biology)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full-text HTML; located the worked example in Figure 1)
**Authority rank:** 1 (peer-reviewed review, Genome Biology 18:186, co-authored by Vinga & Almeida who authored the foundational 2003 review)

**Key Extracted Points:**

1. **Word-vector representation:** Each sequence is transformed into a vector by "counting the number of times each particular word (from W₃) appears within the sequences." The example uses sequences x = "ATGTGTG" and y = "CATGTG" with word size k = 3.
2. **Word set (union over both sequences):** W₃ = {ATG, CAT, GTG, TGT}; per-sequence unique words W_X_3 = {ATG, TGT, GTG}, W_Y_3 = {CAT, ATG, TGT, GTG}.
3. **Count vectors (verbatim from Figure 1):** c_X_3 = (1, 0, 2, 2) and c_Y_3 = (1, 1, 1, 1), in the column order (ATG, CAT, GTG, TGT).
4. **Distance rule:** "This difference is very commonly computed by the Euclidean distance, although any metric can be applied."
5. **Zero/identity property:** "identical sequences yield a distance of 0"; higher value ⇒ more distant sequences.

### Lau, Kläne, Leimeister, Morgenstern et al. — "Interpreting alignment-free sequence comparison: what makes a score a good score?" (NAR Genomics and Bioinformatics, 2022)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full-text HTML; section "Overview of alignment-free methods based on k-mer frequencies" and Table 1)
**Authority rank:** 1 (peer-reviewed)

**Key Extracted Points:**

1. **Frequency definition (verbatim):** "The k-mer frequencies are derived from the counts by dividing each k-mer count by the total number of k-mers in the sequence (i.e. the sequence length minus the k-mer length)." For a sequence of length L over a complete alphabet the number of k-mer windows is L − k + 1.
2. **Distance choice (verbatim):** "Once the frequency vectors have been calculated, various distance measures can be used to estimate similarity, which include the Euclidian, Manhattan, Canberra or Chebyshev distances."
3. **Variable length suitability:** Table 1 lists "Euclidian distance (euclid)" as a metric suitable for sequences of variable length.

### Vinga & Almeida (2003) — "Alignment-free sequence comparison—a review" (Bioinformatics)

**URL:** https://academic.oup.com/bioinformatics/article/19/4/513/218529
**Accessed:** 2026-06-13 (abstract/metadata retrieved via WebFetch; full text paywalled — used only for the high-level statement below, not for the numeric formula)
**Authority rank:** 1 (foundational peer-reviewed review)

**Key Extracted Points:**

1. **Word-composition mapping:** word-composition methods map each sequence into a 4^k-dimensional vector according to k-word frequency; the similarity score is obtained by measures including the Euclidean distance, Pearson correlation, Kullback–Leibler discrepancy and cosine distance.

### Boden, Schöneich, Horwege, Lindner, Leimeister & Morgenstern (2014) — "Fast alignment-free sequence comparison using spaced-word frequencies" (Bioinformatics)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC4080745/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full-text HTML; "Benchmark Set-Up" section)
**Authority rank:** 1 (peer-reviewed; reference method)

**Key Extracted Points:**

1. **Euclidean on relative frequencies (verbatim):** "we applied the Euclidean distance to the relative-frequency vectors obtained with our multiple-pattern approach." Confirms that the standard alignment-free Euclidean distance is taken over relative (normalized) word-frequency vectors, not raw counts.

---

## Documented Corner Cases and Failure Modes

### From Zielezinski et al. (2017)

1. **Identical sequences:** distance = 0 (stated explicitly in Figure 1 description).
2. **Union of words:** the comparison vector spans the union of words occurring in either sequence; words absent from a sequence contribute a 0 component (c_X_3 has a 0 for CAT, which occurs only in y).

### From Lau et al. (2022)

1. **Normalization domain:** frequencies are counts divided by (L − k + 1); this requires at least one k-mer window (L ≥ k), otherwise the frequency vector is empty.

---

## Test Datasets

### Dataset: Zielezinski 2017 Figure 1 worked example

**Source:** Zielezinski et al. (2017), Genome Biology 18:186, Figure 1.

| Parameter | Value |
|-----------|-------|
| Sequence x | ATGTGTG |
| Sequence y | CATGTG |
| k | 3 |
| Word order | ATG, CAT, GTG, TGT |
| c_X (counts) | (1, 0, 2, 2) |
| c_Y (counts) | (1, 1, 1, 1) |
| Total windows x | 5 (= 7 − 3 + 1) |
| Total windows y | 4 (= 6 − 3 + 1) |

**Derived count-based Euclidean distance** (verbatim vectors from source): differences (0, −1, 1, 1) ⇒ √(0² + 1² + 1² + 1²) = √3 ≈ 1.7320508075688772.

**Derived frequency-based Euclidean distance** (counts divided by total windows, per Lau et al. 2022 definition; this is what `KmerAnalyzer.KmerDistance` computes):
f_X = (1/5, 0, 2/5, 2/5) = (0.2, 0.0, 0.4, 0.4); f_Y = (1/4, 1/4, 1/4, 1/4) = (0.25, 0.25, 0.25, 0.25).
differences (−0.05, −0.25, 0.15, 0.15) ⇒ squares (0.0025, 0.0625, 0.0225, 0.0225), sum = 0.11 ⇒ √0.11 = 0.33166247903553997.

### Dataset: Single-substitution small case (derivation, k = 1)

**Source:** Direct derivation from Lau et al. (2022) frequency definition; trivially verifiable.

| Parameter | Value |
|-----------|-------|
| Sequence 1 | AAAA |
| Sequence 2 | AAAT |
| k | 1 |
| Words | A, T |
| f_1 | A=4/4=1.0, T=0 |
| f_2 | A=3/4=0.75, T=1/4=0.25 |

Euclidean distance = √((1.0−0.75)² + (0−0.25)²) = √(0.0625 + 0.0625) = √0.125 = 0.3535533905932738.

---

## Assumptions

1. **ASSUMPTION: Case folding / alphabet handling.** Authoritative sources work over a fixed nucleotide alphabet and do not specify case sensitivity. The implementation upper-cases input before counting (via `CountKmers`), so mixed-case inputs are treated as the same k-mer. This is a benign normalization that does not change the source-defined output for canonical (already upper-case) inputs.
2. **ASSUMPTION: Empty / too-short input.** Sources state frequencies require L ≥ k but do not define the distance when a sequence has no k-mer windows. The implementation returns an empty frequency vector for such inputs, so the distance equals the Euclidean norm of the other sequence's frequency vector, and 0 when both are empty. This is the natural extension (a sequence with no words is the zero vector); marked as an assumption because no source defines it explicitly.

---

## Recommendations for Test Coverage

1. **MUST Test:** Zielezinski 2017 Figure 1 example (x="ATGTGTG", y="CATGTG", k=3) ⇒ frequency-based distance √0.11 ≈ 0.3316624790. — Evidence: Zielezinski et al. (2017) Fig. 1 + Lau et al. (2022) frequency definition.
2. **MUST Test:** Identical sequences ⇒ distance exactly 0. — Evidence: Zielezinski et al. (2017) Fig. 1.
3. **MUST Test:** Single-substitution k=1 derivation (AAAA vs AAAT) ⇒ √0.125 ≈ 0.3535533906. — Evidence: Lau et al. (2022) frequency definition.
4. **MUST Test:** Symmetry d(x,y) = d(y,x). — Rationale: Euclidean distance is a metric (symmetric).
5. **SHOULD Test:** Non-overlapping word sets (e.g. all-A vs all-T, k≥2) ⇒ √(1²+1²)=√2 (each sequence is a single distinct k-mer with frequency 1). — Rationale: maximal-disjoint case has a closed form.
6. **SHOULD Test:** k > min(length) for one sequence ⇒ that sequence has an empty vector; distance = norm of the other's frequency vector. — Rationale: documented short-input behavior (ASSUMPTION).
7. **COULD Test:** k ≤ 0 throws ArgumentOutOfRangeException. — Rationale: input validation contract inherited from CountKmers.
8. **COULD Test:** Case-insensitivity (lower vs upper) yields identical distance. — Rationale: documents the normalization assumption.

---

## References

1. Zielezinski A, Vinga S, Almeida J, Karlowski WM. (2017). Alignment-free sequence comparison: benefits, applications, and tools. Genome Biology 18:186. https://pmc.ncbi.nlm.nih.gov/articles/PMC5627421/ (DOI: 10.1186/s13059-017-1319-7)
2. Lau AK, et al. (2022). Interpreting alignment-free sequence comparison: what makes a score a good score? NAR Genomics and Bioinformatics. https://pmc.ncbi.nlm.nih.gov/articles/PMC9442500/
3. Vinga S, Almeida J. (2003). Alignment-free sequence comparison—a review. Bioinformatics 19(4):513–523. https://academic.oup.com/bioinformatics/article/19/4/513/218529 (DOI: 10.1093/bioinformatics/btg005)
4. Boden M, et al. (2014). Fast alignment-free sequence comparison using spaced-word frequencies. Bioinformatics 30(14):1991–1999. https://pmc.ncbi.nlm.nih.gov/articles/PMC4080745/

---

## Change History

- **2026-06-13**: Initial documentation.
