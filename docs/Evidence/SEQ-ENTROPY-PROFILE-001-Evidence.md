# Evidence Artifact: SEQ-ENTROPY-PROFILE-001

**Test Unit ID:** SEQ-ENTROPY-PROFILE-001
**Algorithm:** Shannon Entropy Profile (sliding-window per-symbol Shannon entropy)
**Date Collected:** 2026-06-14

---

## Online Sources

### Shannon (1948) — A Mathematical Theory of Communication (primary source)

**URL:** https://en.wikipedia.org/wiki/A_Mathematical_Theory_of_Communication (bibliographic confirmation) ; DOI https://doi.org/10.1002/j.1538-7305.1948.tb01338.x
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary paper)

**Retrieval:** Web search `Shannon 1948 "A Mathematical Theory of Communication" entropy H = - sum p log p Theorem 2 choice uncertainty`, then WebFetch of the Wikipedia article on the paper to confirm full bibliographic details.

**Key Extracted Points:**

1. **Bibliographic details (verbatim from fetched page):** Author Claude E. Shannon; year 1948; title "A Mathematical Theory of Communication"; journal Bell System Technical Journal; Part 1 Vol. 27, Issue 3 (July 1948), pp. 379–423; DOI 10.1002/j.1538-7305.1948.tb01338.x.
2. **Entropy as a measure of choice/uncertainty:** The paper introduces entropy H = −Σ pᵢ log pᵢ as the measure of information, choice and uncertainty of a discrete source (the originating definition cited by all downstream sources below).

### Entropy (information theory) — Wikipedia (citing Shannon 1948)

**URL:** https://en.wikipedia.org/wiki/Entropy_(information_theory)
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing the primary Shannon 1948)

**Retrieval:** WebFetch of the URL with a prompt requesting the formal H definition, log base, units, and the maximum-entropy property.

**Key Extracted Points (closely paraphrased from fetched text):**

1. **Formal definition:** "H(X) = −∑ₓ p(x) log_b p(x)", where p(x) is the probability of each outcome and b is the logarithm base.
2. **Units:** "When b = 2, the entropy is measured in bits (also called shannons)."
3. **Maximum entropy:** "The maximal entropy of an event with n different outcomes is log_b(n): it is attained by the uniform probability distribution." For n equally likely outcomes, entropy equals log₂(n) bits.

### Entropy-Based Biological Sequence Study — IntechOpen (DNA application)

**URL:** https://www.intechopen.com/chapters/75997
**Accessed:** 2026-06-14
**Authority rank:** 3–4 (peer-reviewed open-access book chapter applying Shannon entropy to DNA, citing primaries)

**Retrieval:** Web search `Shannon entropy sequence sliding window DNA bioinformatics local complexity bits formula`, then WebFetch of the chapter URL.

**Key Extracted Points (closely paraphrased from fetched text):**

1. **DNA entropy formula (Eq. 3):** "yᵢ = −Σⱼ pᵢⱼ log pᵢⱼ", where pᵢⱼ is the probability of occurrence of the j-th genetic letter; applied to the four-nucleotide alphabet {A, C, G, T}.
2. **Maximum for DNA:** "For DNA sequences with 4 nucleotides, maximum entropy equals 2 bits (log₂ 4 = 2). This occurs when all nucleotides appear with equal probability."
3. **Sliding window:** The method uses "a sliding 'counter' of width W over [the] DNA sequence", generating per-window symbol frequency counts that become the probability distribution for the entropy calculation.

---

## Documented Corner Cases and Failure Modes

### From Wikipedia (Entropy, information theory) / Shannon (1948)

1. **Zero-probability convention:** Terms with pᵢ = 0 contribute 0 to the sum (0·log 0 ≡ 0 by the limit x log x → 0). Homopolymer (single symbol) windows therefore yield H = 0.
2. **Maximum at uniform distribution:** A window in which every symbol is equally frequent attains H = log₂(k) where k is the number of distinct symbols present (2 bits for the full 4-letter DNA alphabet).

### From IntechOpen chapter

1. **Window narrower than sequence required:** entropy is computed per window of width W slid along the sequence; if W exceeds the sequence length no full window exists.

---

## Test Datasets

### Dataset: Hand-derived per-window entropy values (H = −Σ pᵢ log₂ pᵢ)

**Source:** Derived directly from the Shannon (1948) / Wikipedia definition and verified numerically.

| Window | Symbol counts | H (bits), derivation | H value |
|--------|---------------|----------------------|---------|
| `AAAA` | A=4 | −(1·log₂1) | 0.0 |
| `AATT` | A=2,T=2 | −2·(½·log₂½) | 1.0 |
| `ATGC` | A=T=G=C=1 | −4·(¼·log₂¼) = log₂4 | 2.0 |
| `AAAT` | A=3,T=1 | −(¾·log₂¾ + ¼·log₂¼) | 0.8112781244591328 |
| `AATG` | A=2,T=1,G=1 | −(½·log₂½ + 2·¼·log₂¼) | 1.5 |
| `GCAA` | G=1,C=1,A=2 | −(2·¼·log₂¼ + ½·log₂½) | 1.5 |
| `AAATTC` | A=3,T=2,C=1 | −(½log₂½ + ⅓log₂⅓ + ⅙log₂⅙) | 1.4591479170272448 |

### Dataset: Sliding-window profiles

**Source:** Derived by applying the per-window definition above at each window offset.

| Sequence | windowSize | stepSize | Window offsets | Expected profile (bits) |
|----------|-----------|----------|----------------|--------------------------|
| `AAATGC` | 4 | 1 | 0:`AAAT`, 1:`AATG`, 2:`ATGC` | [0.8112781244591328, 1.5, 2.0] |
| `AAATGCAA` | 4 | 2 | 0:`AAAT`, 2:`ATGC`, 4:`GCAA` | [0.8112781244591328, 2.0, 1.5] |

---

## Assumptions

1. **ASSUMPTION: Per-symbol (k=1) frequencies over the IUPAC letter alphabet** — The implementation computes pᵢ from single-character (mono-nucleotide) frequencies of the letters present in the window (case-folded, non-letters ignored). The cited sources define H over any symbol probability distribution; the choice of the mono-symbol alphabet (rather than k-mer/block entropy) is the implementation's modelling choice and is consistent with the four-letter DNA application (max 2 bits) in the IntechOpen chapter. This does not change the formula, only the alphabet over which pᵢ is taken.

---

## Recommendations for Test Coverage

1. **MUST Test:** Per-window values for uniform (2 bits), two-symbol equal (1 bit), homopolymer (0 bits), and skewed (3:1 → 0.8112781…) windows — Evidence: Wikipedia/Shannon 1948 definition; IntechOpen Eq. 3.
2. **MUST Test:** Full sliding profile with stepSize=1 and stepSize>1, asserting exact per-window values and the number/order of windows — Evidence: IntechOpen sliding-window method.
3. **SHOULD Test:** windowSize > length → empty profile; windowSize == length → single value — Rationale: documented "window narrower than sequence" corner case.
4. **SHOULD Test:** Maximum-entropy invariant (H ≤ log₂ k) and non-negativity (H ≥ 0) — Rationale: Wikipedia maximum-entropy property.
5. **COULD Test:** Case-insensitivity (lowercase input yields same profile) — Rationale: implementation case-folds before counting.

---

## References

1. Shannon, C. E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal, 27(3), 379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x
2. Wikipedia contributors. Entropy (information theory). https://en.wikipedia.org/wiki/Entropy_(information_theory) (accessed 2026-06-14).
3. Wikipedia contributors. A Mathematical Theory of Communication. https://en.wikipedia.org/wiki/A_Mathematical_Theory_of_Communication (accessed 2026-06-14).
4. Entropy-Based Biological Sequence Study. IntechOpen. https://www.intechopen.com/chapters/75997 (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
