# Evidence Artifact: KMER-STATS-001

**Test Unit ID:** KMER-STATS-001
**Algorithm:** K-mer Statistics (comprehensive k-mer composition statistics)
**Date Collected:** 2026-06-14

---

## Online Sources

### Wikipedia — K-mer

**URL:** https://en.wikipedia.org/wiki/K-mer
**Accessed:** 2026-06-14 (retrieved via WebFetch of the article page)
**Authority rank:** 4 (Wikipedia citing primary sources; its formulas are corroborated below by ranks 1–3)

**Key Extracted Points:**

1. **Total k-mer count:** "a sequence of length L will have L − k + 1 k-mers" (retrieved verbatim from the page). This is the number of overlapping length-k windows.
2. **K-mer universe size:** "there exist n^k total possible k-mers, where n is number of possible monomers (e.g. four in the case of DNA)."
3. **Worked example AGAT:** a sequence AGAT "contains four monomers (A, G, A, and T), three 2-mers (AG, GA, AT), two 3-mers (AGA and GAT), and one 4-mer (AGAT)."
4. **Worked example GTAGAGCTGT (count table):** the article's example table enumerates, for a 10-character sequence:
   - k=1: G,T,A,C — total 10, distinct 4 (G appears 4×, T 3×, A 2×, C 1×).
   - k=2: GT,TA,AG,GA,AG,GC,CT,TG,GT — total 9, distinct 7 (GT and AG each appear twice).
   - k=3: GTA,TAG,AGA,GAG,AGC,GCT,CTG,TGT — total 8, distinct 8 (all distinct).
   These totals match L−k+1 for L=10: 10, 9, 8.

### BioInfoLogics — k-mer counting, part I (Introduction)

**URL:** https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 4 (technical tutorial corroborating the L−k+1 formula and distinct/unique distinction; used only for the worked count table)

**Key Extracted Points:**

1. **Total k-mer count:** "Any sequence of length `L` will contain `L - k + 1` _k-mers_." For ATCGATCAC (L=9), k=3 → 7 total k-mers (retrieved verbatim).
2. **Non-canonical count table for ATCGATCAC, k=3:** ATC=2, TCG=1, CGA=1, GAT=1, TCA=1, CAC=1 → **distinct = 6**, **unique (count==1) = 5**.
3. **Distinct vs unique distinction:** "distinct" counts each different k-mer once; "unique" k-mers are those occurring exactly once. (Used here only for distinct count = 6 in AnalyzeKmers; the unique-count concept is KMER-UNIQUE-001.)

### Spectral concepts in genome informational analysis (Bonnici & Manca / Manca, k-spectrum work)

**URL:** https://arxiv.org/abs/2106.15351 (abstract retrieved via WebFetch); the explicit k-entropy formula was extracted via WebSearch over the same paper's full text on 2026-06-14
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed/preprint on genome k-mer information theory)

**Key Extracted Points:**

1. **k-entropy definition (verbatim, via retrieved text):** "The k-entropy of a genome is defined as the Shannon entropy of the probability distribution assigning to any k-mer its frequency, calculated as: E_k(G) = Σ p(α) log₂ p(α), where p(α) = mult(α) / (|G| − k + 1)." The Shannon entropy is the negation of that sum: E_k(G) = −Σ p(α) log₂ p(α).
2. **Probability of a k-mer α:** p(α) = multiplicity of α divided by the total number of k-mer windows (|G| − k + 1) — i.e. the relative frequency over the L−k+1 overlapping windows.

### arXiv 2511.05300 — Entropy–Rank Ratio (Shannon entropy of a sequence via k-mers)

**URL:** https://arxiv.org/html/2511.05300 (full text retrieved via WebFetch) plus corroborating WebSearch extract on 2026-06-14
**Accessed:** 2026-06-14
**Authority rank:** 1 (preprint, peer-reviewable; used to corroborate the single-sequence k-mer Shannon entropy form)

**Key Extracted Points:**

1. **Single-sequence k-mer Shannon entropy (corroborated form):** H_k(s) = −∑_i p_i log₂(p_i), where p_i is the relative frequency of the i-th distinct k-mer in sequence s (its count divided by the total number of k-mers). The paper uses base-2 logarithm and the convention 0·log(0) = 0.
2. **Base-2 logarithm / bits:** entropy is "a measure of the informational complexity (in bits)"; the log base is 2.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia — K-mer

1. **k > L (window does not fit):** with the L−k+1 formula, when k exceeds the sequence length the count L−k+1 ≤ 0, i.e. there are no k-mers. The statistics over an empty k-mer multiset are therefore empty.
2. **Homopolymer / single distinct k-mer:** when only one distinct k-mer exists (e.g. AAAA, k=2 → AA×3), the frequency distribution has a single component p=1, so Shannon entropy = −1·log₂(1) = 0 (minimum diversity).

### From Spectral concepts paper (k-entropy)

1. **All k-mers distinct (max diversity at fixed total):** when every one of the L−k+1 windows is a distinct k-mer, each p = 1/(L−k+1) and H = log₂(L−k+1) (e.g. GTAGAGCTGT k=3 → 8 distinct → H = log₂8 = 3 bits).

---

## Test Datasets

### Dataset: GTAGAGCTGT (Wikipedia K-mer example)

**Source:** Wikipedia — K-mer, https://en.wikipedia.org/wiki/K-mer (worked example table)

| Parameter | Value |
|-----------|-------|
| Sequence | GTAGAGCTGT (L=10) |
| k=1 → total, distinct, max, min | 10, 4, 4 (G), 1 (C) |
| k=1 average count | 10/4 = 2.5 |
| k=1 Shannon entropy (bits) | −(0.4·log₂0.4 + 0.3·log₂0.3 + 0.2·log₂0.2 + 0.1·log₂0.1) = 1.846439344671… |
| k=2 → total, distinct, max, min | 9, 7, 2 (GT, AG), 1 |
| k=2 average count | 9/7 = 1.285714… (rounded 1.29) |
| k=3 → total, distinct, max, min | 8, 8, 1, 1 |
| k=3 Shannon entropy (bits) | log₂8 = 3.0 (8 equiprobable k-mers) |

### Dataset: ATCGATCAC (BioInfoLogics example)

**Source:** BioInfoLogics — k-mer counting part I, https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/

| Parameter | Value |
|-----------|-------|
| Sequence | ATCGATCAC (L=9) |
| k=3 counts | ATC=2, TCG=1, CGA=1, GAT=1, TCA=1, CAC=1 |
| k=3 → total, distinct, max, min | 7, 6, 2 (ATC), 1 |
| k=3 average count | 7/6 = 1.166… (rounded 1.17) |
| k=3 Shannon entropy (bits) | −(2/7·log₂(2/7) + 5·(1/7·log₂(1/7))) = 2.521640636343… |

### Dataset: AAAA homopolymer (corner case)

**Source:** Derived from the k-entropy definition (single-component distribution; Spectral concepts paper)

| Parameter | Value |
|-----------|-------|
| Sequence | AAAA (L=4), k=2 |
| counts | AA=3 |
| total, distinct, max, min | 3, 1, 3, 3 |
| average count | 3.0 |
| Shannon entropy | −1·log₂1 = 0 |

---

## Assumptions

1. **ASSUMPTION: AverageCount is rounded to 2 decimal places.** The literature defines average k-mer multiplicity as total/distinct = (L−k+1)/distinct, but does not prescribe a display rounding. The repository implementation rounds `AverageCount` to 2 decimals (`Math.Round(averageCount, 2)`). This is a non-correctness-affecting presentation choice (the underlying ratio is exact and verifiable); tests assert the rounded value to remain consistent with the public contract. The mathematically exact ratio is also confirmed in each dataset row.
2. **ASSUMPTION: Entropy is reported unrounded in bits (log base 2).** Sources fix log base 2 (bits); the implementation returns the unrounded double. Tests assert exact values with `.Within(1e-10)`.

---

## Recommendations for Test Coverage

1. **MUST Test:** TotalKmers = L−k+1 for the Wikipedia worked examples (GTAGAGCTGT k=1/2/3) — Evidence: Wikipedia K-mer ("L − k + 1 k-mers"; example table totals 10/9/8).
2. **MUST Test:** UniqueKmers field = distinct k-mer count (GTAGAGCTGT k=2 → 7 distinct; ATCGATCAC k=3 → 6 distinct) — Evidence: Wikipedia example table; BioInfoLogics count table.
3. **MUST Test:** MaxCount / MinCount equal the observed extremes of the multiplicity distribution (GTAGAGCTGT k=1 → max 4, min 1) — Evidence: Wikipedia example table.
4. **MUST Test:** AverageCount = total/distinct (GTAGAGCTGT k=2 → 9/7 ≈ 1.29) — Evidence: derived from totals in the Wikipedia table.
5. **MUST Test:** Entropy = −Σ p log₂ p over k-mer frequencies (GTAGAGCTGT k=1 → 1.84643934…; k=3 → 3.0; ATCGATCAC k=3 → 2.52164064…) — Evidence: k-entropy formula E_k = −Σ p(α) log₂ p(α), p(α)=mult/(L−k+1).
6. **SHOULD Test:** Homopolymer (single distinct k-mer) → entropy 0, max==min==total — Rationale: documented minimum-diversity corner case.
7. **SHOULD Test:** All-distinct sequence → entropy = log₂(distinct); INV max==min==1 — Rationale: documented maximum-diversity corner case.
8. **SHOULD Test:** Empty sequence and k > L → all-zero KmerStatistics — Rationale: L−k+1 ≤ 0 produces no k-mers.
9. **COULD Test:** Case-insensitivity (lower-case input yields identical statistics) — Rationale: implementation upper-cases internally.
10. **COULD Test:** k ≤ 0 throws ArgumentOutOfRangeException — Rationale: k-mer length must be positive.
11. **MUST Test (invariants):** TotalKmers == sum of all k-mer counts; UniqueKmers == distinct count from CountKmers; cross-checked independently — Evidence: definitions of total/distinct.

---

## References

1. Wikipedia contributors. (2026). K-mer. Wikipedia. https://en.wikipedia.org/wiki/K-mer
2. Clavijo, B. (2018). k-mer counting, part I: Introduction. BioInfoLogics. https://bioinfologics.github.io/post/2018/09/17/k-mer-counting-part-i-introduction/
3. Manca, V. et al. (2021). Spectral concepts in genome informational analysis. arXiv:2106.15351. https://arxiv.org/abs/2106.15351
4. Entropy–Rank Ratio: A Novel Entropy–Based Perspective for DNA Complexity and Classification. (2025). arXiv:2511.05300. https://arxiv.org/html/2511.05300

---

## Change History

- **2026-06-14**: Initial documentation.
