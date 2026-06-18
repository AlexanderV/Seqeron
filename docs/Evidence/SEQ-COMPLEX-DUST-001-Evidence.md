# Evidence Artifact: SEQ-COMPLEX-DUST-001

**Test Unit ID:** SEQ-COMPLEX-DUST-001
**Algorithm:** DUST Score (triplet-frequency low-complexity score of Morgulis et al. 2006 SDUST/DUST)
**Date Collected:** 2026-06-14

---

## Online Sources

### Morgulis, Gertz, Schäffer & Agarwala (2006). A fast and symmetric DUST implementation to mask low-complexity DNA sequences. J Comput Biol 13(5):1028–1040.

**URL:** https://pubmed.ncbi.nlm.nih.gov/16796549/
**Accessed:** 2026-06-14 (fetched via WebFetch of the PubMed abstract page)
**Authority rank:** 1 (peer-reviewed primary paper defining the algorithm)

**Key Extracted Points:**

1. **Citation:** WebFetch of the PubMed page returned: Title "A fast and symmetric DUST implementation to mask low-complexity DNA sequences"; Authors Aleksandr Morgulis, E Michael Gertz, Alejandro A Schäffer, Richa Agarwala; Journal of Computational Biology; 13(5):1028-40; 2006; DOI 10.1089/cmb.2006.13.1028.
2. **Method:** WebSearch result block (researchr/PubMed summary) states DUST "is a heuristic algorithm that employs a scoring function based on counting nucleotide triplet frequencies in 64-base windows" and that the new implementation "uses the same function to assign a complexity score to a sequence" while changing only the masking rule (now symmetric and context-insensitive).

### Li, H. (2025). Finding low-complexity DNA sequences with longdust (arXiv:2509.07357).

**URL:** https://arxiv.org/pdf/2509.07357
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 1 (primary preprint by the author of minimap2/samtools; restates the SDUST score verbatim)

**Key Extracted Points:**

1. **Score formula (verbatim from WebFetch):** the SDUST complexity score for a string `x` of length `L` is `∑_t c_x(t)(c_x(t)-1)/2 / (L-2)`, where `c_x(t)` is the count of triplet `t` in `x` and the sum runs over all triplets.
2. **Normalization:** division is by `(L-2)`, i.e. the number of triplets ℓ = L−2 (k = 3 hardcoded).
3. **Thresholding form:** the internal scoring function is written `S_S(c_x) = 1/ℓ(x) ∑_t [c_x(t)(c_x(t)-1)/2] − T` (WebSearch result block), confirming both the `1/ℓ` normalization and the threshold `T` subtracted from the raw score.
4. **Default parameters:** default window size `w = 64`; threshold corresponds to complexity level 20, score 2.0 (WebFetch).
5. **Direction:** repeated triplets raise `∑ c(c-1)/2`, so a HIGH score indicates LOW complexity; low-complexity (high-scoring) regions are the ones masked (WebSearch result block).

### lh3/sdust — reference C implementation (`sdust.c`, master).

**URL:** https://raw.githubusercontent.com/lh3/sdust/master/sdust.c
**Accessed:** 2026-06-14 (fetched via WebFetch of the raw source file)
**Authority rank:** 3 (reference implementation by Heng Li, reimplementing the symmetric DUST algorithm)

**Key Extracted Points:**

1. **Incremental score accumulation (verbatim line):** `++*L, *rw += cw[t]++, *rv += cv[t]++;` — on adding triplet `t`, the *current* count `cw[t]` is added to the running score `rw` before the count is incremented. Summing the pre-increment counts 0+1+2+…+(c−1) over all occurrences of a triplet yields exactly `c(c−1)/2`, so `rw = ∑_t c(c−1)/2` and `L` counts the triplets added.
2. **Threshold comparison (verbatim line):** `if (rw * 10 > L * T)` — equivalent to `rw / L > T / 10`; with default `T = 20` this is `score > 2.0`, confirming the threshold value 2.0.
3. **Default parameters (verbatim):** `int W = 64;` and `int T = 20;` — default window width 64 and threshold level 20 (score 2.0).

---

## Documented Corner Cases and Failure Modes

### From Li (2025) / Morgulis et al. (2006)

1. **Below triplet length (L < 3):** the score `∑ c(c-1)/2 / (L-2)` is undefined when L−2 ≤ 0 (no triplets exist); the operation has no defined complexity for such inputs.
2. **All-distinct triplets:** when every triplet occurs exactly once, every `c(c-1)/2 = 0`, so the score is 0 (maximum complexity).
3. **Maximally repetitive (single triplet):** a homopolymer of length L has one triplet repeated L−2 times, giving the maximal score `(L−2)(L−3)/2 / (L−2) = (L−3)/2`.

---

## Test Datasets

### Dataset: Hand-derived worked examples (k = 3, divisor = number of triplets = L−2)

**Source:** Derived directly from the Li (2025) formula `∑_t c_t(c_t−1)/2 / (L−2)` and the lh3/sdust accumulation.

| Input | L | Triplets (count) | Σ c(c−1)/2 | L−2 | Score |
|-------|---|------------------|------------|-----|-------|
| `ATGC` | 4 | ATG=1, TGC=1 | 0 | 2 | 0.0 |
| `ACGTACGT` | 8 | ACG=2,CGT=2,GTA=1,TAC=1 | 1+1=2 | 6 | 0.3333333333… |
| `AAAAAA` | 6 | AAA=4 | 4·3/2=6 | 4 | 1.5 |
| `ACACACAC` | 8 | ACA=3,CAC=3 | 3+3=6 | 6 | 1.0 |
| `AAAAAAAAAA` | 10 | AAA=8 | 8·7/2=28 | 8 | 3.5 |

### Dataset: Threshold reference

**Source:** lh3/sdust `if (rw * 10 > L * T)` with `T = 20`; Li (2025) score 2.0 / level 20.

| Parameter | Value |
|-----------|-------|
| Default window size | 64 |
| Default threshold (mask if score >) | 2.0 |
| Triplet/word size k | 3 |

---

## Assumptions

1. **ASSUMPTION: General word size `wordSize`** — The paper and reference implementation hardcode k = 3 (triplets). The repository method exposes a `wordSize` parameter; for `wordSize = w` the normalization generalizes to the number of words `L − w + 1` (= `L − 2` when w = 3). This generalization is consistent with the formula but only k = 3 is source-backed; tests assert exact source-derived values only for k = 3.
2. **ASSUMPTION: Input shorter than one word (L < wordSize)** — Neither source defines a score when no word exists. The implementation returns 0 (no repeats ⇒ minimal complexity); this is a defined-output convention, not a source value.

---

## Recommendations for Test Coverage

1. **MUST Test:** Homopolymer `AAAAAA` ⇒ score 1.5 (k=3) — Evidence: Li (2025) formula, hand-derived.
2. **MUST Test:** `ACGTACGT` ⇒ score 0.333… (k=3) — Evidence: Li (2025) formula, hand-derived.
3. **MUST Test:** All-distinct-triplet input `ATGC` ⇒ score 0.0 — Evidence: Li (2025) formula; all c(c−1)/2 = 0.
4. **MUST Test:** Repetitive dinucleotide `ACACACAC` ⇒ score 1.0 — Evidence: Li (2025) formula, hand-derived.
5. **MUST Test:** Longer homopolymer `AAAAAAAAAA` ⇒ score 3.5 — Evidence: Li (2025) formula, hand-derived.
6. **MUST Test:** DnaSequence and string overloads agree (same input ⇒ same score) — Evidence: both wrap one core.
7. **SHOULD Test:** Case-insensitivity (string overload upper-cases) — Rationale: DnaSequence normalizes to upper-case; string overload documents `ToUpperInvariant`.
8. **SHOULD Test:** Null DnaSequence ⇒ ArgumentNullException; null/empty string ⇒ 0 — Rationale: documented validation, matches sibling methods.
9. **COULD Test:** Input shorter than wordSize ⇒ 0 — Rationale: defined-output convention (ASSUMPTION 2).

---

## References

1. Morgulis A, Gertz EM, Schäffer AA, Agarwala R. (2006). A fast and symmetric DUST implementation to mask low-complexity DNA sequences. Journal of Computational Biology 13(5):1028–1040. https://doi.org/10.1089/cmb.2006.13.1028 (abstract retrieved via https://pubmed.ncbi.nlm.nih.gov/16796549/)
2. Li H. (2025). Finding low-complexity DNA sequences with longdust. arXiv:2509.07357. https://arxiv.org/pdf/2509.07357
3. Li H. sdust — Symmetric DUST for finding low-complexity regions in DNA sequences (reference C implementation). https://raw.githubusercontent.com/lh3/sdust/master/sdust.c (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
