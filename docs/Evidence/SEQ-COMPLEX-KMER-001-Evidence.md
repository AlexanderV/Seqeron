# Evidence Artifact: SEQ-COMPLEX-KMER-001

**Test Unit ID:** SEQ-COMPLEX-KMER-001
**Algorithm:** K-mer Entropy (Shannon entropy of the overlapping k-mer frequency distribution)
**Date Collected:** 2026-06-14

---

## Online Sources

### Li, H. (2025). Finding low-complexity DNA sequences with longdust (arXiv:2509.07357)

**URL:** https://arxiv.org/pdf/2509.07357
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 1 (peer-reviewed-style preprint by the author of minimap2/samtools; primary method for k-mer-frequency low-complexity detection)

**Key Extracted Points:**

1. **K-mer extraction:** k-mers are extracted with a sliding window at every position; "For a sequence of length L, there are L − k + 1 overlapping k-mers, where k is the window size." (retrieved via WebFetch query asking whether k-mers are overlapping and how many there are)
2. **Entropy formula:** "Shannon entropy is defined as H = -Σ p_i log₂(p_i), where p_i is the frequency of the i-th k-mer."
3. **Probability estimate:** "If n_i represents the count of k-mer i and N = L − k + 1 is the total number of k-mers, then p_i = n_i/N, and the entropy sums across all observed k-mers."
4. **Complexity interpretation:** low-complexity sequences have skewed k-mer distributions (few k-mers dominate) → low entropy; high-complexity sequences have uniform distributions → high entropy.

### Çakır et al. (2025). Entropy–Rank Ratio: An Entropy-Based Perspective for DNA Complexity (arXiv:2511.05300)

**URL:** https://arxiv.org/html/2511.05300
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 1 (peer-reviewed-style preprint)

**Key Extracted Points:**

1. **Entropy formula (Eq. 17):** S(w) = −Σ_{j=1}^{λ} b_j log(b_j), where b_j is the relative frequency of each element and λ the number of distinct symbols.
2. **Probability (Def. 2.5):** b_j = a_j / M where a_j is the occurrence count and M the total number of tuples, so Σ b_j = 1.
3. **Logarithm base:** "By convention, log denotes the base-2 logarithm," yielding entropy in **bits**; maximum entropy for single nucleotides is log₂(4) = 2 bits.
4. **Saturation:** for very long uniform i.i.d. sequences the entropy converges to log(λ) (the maximum).

### "About Shannon's Entropy" (citing Shannon 1948)

**URL:** https://tcosmo.github.io/2019/04/21/shannon-entropy.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 4 (secondary exposition of the Shannon 1948 primary; used only for the well-established bounds)

**Key Extracted Points:**

1. **Formula:** H(X) = -Σ p_i log(p_i) = Σ p_i log(1/p_i).
2. **Bounds:** 0 ≤ H(X) ≤ log(k); H(X) = 0 when the distribution is deterministic ("one entry with probability one and all the others zero"); H(X) = log(k) when the distribution is uniform over k outcomes.

### Wikipedia, "Entropy (information theory)" (citing Shannon 1948)

**URL:** https://en.wikipedia.org/wiki/Entropy_(information_theory)
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 4 (uses Shannon 1948 primary)

**Key Extracted Points:**

1. **Definition:** H(X) := −Σ_{x∈X} p(x) log p(x), log base depends on application.
2. **Maximum:** H_n(p_1,…,p_n) ≤ H_n(1/n,…,1/n); maximum entropy equals log_b(n) for n equiprobable outcomes.
3. **Minimum:** "The minimum surprise is when p = 0 (impossibility) or p = 1 (certainty) and the entropy is zero bits."
4. **Primary attribution:** "The concept of information entropy was introduced by Claude Shannon in his 1948 paper 'A Mathematical Theory of Communication.'"

---

## Documented Corner Cases and Failure Modes

### From Li (2025), longdust

1. **Single repeated k-mer (homopolymer / tandem repeat):** the distribution collapses to one k-mer with p = 1, giving H = 0 (lowest complexity).
2. **All-distinct k-mers:** each k-mer appears once (p_i = 1/N), giving the maximum H = log₂(N) for that sequence.

### From Shannon 1948 (via Wikipedia / exposition)

1. **Deterministic distribution → H = 0**; uniform distribution → H = log_b(n) (n distinct symbols).

### Spec-undefined / contract-resolved (not in sources; see Assumptions)

1. **L < k (sequence shorter than k):** no k-mers exist; resolved to return 0 (consistent with siblings in `SequenceComplexity`).
2. **k < 1:** invalid window; resolved to throw `ArgumentOutOfRangeException` (matches sibling guards).
3. **null sequence:** `DnaSequence` overload throws `ArgumentNullException`; `string` overload returns 0 for null/empty (matches sibling string overloads).

---

## Test Datasets

### Dataset: Hand-derived worked examples from the formula H = −Σ (n_i/N) log₂(n_i/N), N = L−k+1

**Source:** derivation from Li (2025) formula H = -Σ p_i log₂ p_i with p_i = n_i/(L−k+1).

| Input | k | k-mers (overlapping) | Counts | N=L−k+1 | H (bits) |
|-------|---|----------------------|--------|---------|----------|
| `ACGT` | 1 | A,C,G,T | each 1 | 4 | log₂(4) = 2.0 |
| `ACGT` | 2 | AC,CG,GT | each 1 | 3 | log₂(3) = 1.5849625007211562 |
| `ATATAT` | 2 | AT,TA,AT,TA,AT | AT=3,TA=2 | 5 | −(0.6·log₂0.6 + 0.4·log₂0.4) = 0.9709505944546686 |
| `AAAA` | 2 | AA,AA,AA | AA=3 | 3 | 0.0 |
| `AAACGT` | 2 | AA,AA,AC,CG,GT | AA=2,AC=1,CG=1,GT=1 | 5 | −(0.4·log₂0.4 + 3·0.2·log₂0.2) = 1.9219280948873623 |
| `AC` | 5 | (none, L<k) | — | — | 0.0 |

Derivation of `ATATAT`,k=2: this is the binary entropy of p=0.6: H = −0.6·log₂0.6 − 0.4·log₂0.4 = 0.6·0.7369655942 + 0.4·1.3219280949 = 0.4421793565 + 0.5287712380 = 0.9709505945.
Derivation of `AAACGT`,k=2 (N=5; p=2/5,1/5,1/5,1/5): H = −[0.4·log₂0.4 + 3·(0.2·log₂0.2)] = 0.4·1.3219280949 + 0.6·2.3219280949 = 0.5287712380 + 1.3931568569 = 1.9219280949 (exact 1.9219280948873623, = log₂5 − 0.4).

---

## Assumptions

1. **ASSUMPTION: L < k returns 0** — No source numerically specifies the L < k case (no k-mers exist). Resolved by the contract used across `SequenceComplexity` siblings (`CalculateLinguisticComplexity`, `EstimateCompressionRatio`) which return 0 for empty/too-short input. Non-correctness-affecting beyond this boundary: the entropy of an empty multiset is conventionally 0.
2. **ASSUMPTION: invalid k (< 1) throws `ArgumentOutOfRangeException`; null `DnaSequence` throws `ArgumentNullException`; null/empty string returns 0** — Failure modes are not specified by the entropy literature (they are library-API contract). Resolved to match sibling method guards in the same class. API-shape only; does not change entropy values for valid inputs.

---

## Recommendations for Test Coverage

1. **MUST Test:** `ACGT`,k=1 → 2.0 (uniform, H = log₂4). — Evidence: Çakır 2025 (max entropy = log₂4); Shannon uniform bound.
2. **MUST Test:** `ACGT`,k=2 → log₂3 ≈ 1.5849625 (all-distinct k-mers, H = log₂N). — Evidence: Li 2025 all-distinct case; Shannon uniform bound.
3. **MUST Test:** `ATATAT`,k=2 → 0.9709505945 (non-uniform; binary entropy of 0.6). — Evidence: Li 2025 formula H = −Σ p_i log₂ p_i, p_i = n_i/(L−k+1).
4. **MUST Test:** `AAAA`,k=2 → 0.0 (deterministic distribution). — Evidence: Shannon H=0 for certainty; Li 2025 skewed-distribution → low entropy.
5. **MUST Test:** `AAACGT`,k=2 → 1.9219280949 (mixed counts; = log₂5 − 0.4). — Evidence: Li 2025 formula.
6. **SHOULD Test:** L < k returns 0. — Rationale: documented boundary (no k-mers).
7. **SHOULD Test:** invalid k / null throw; null/empty string → 0. — Rationale: documented failure modes (contract).
8. **SHOULD Test:** string and DnaSequence overloads agree (case-insensitive). — Rationale: API consistency; DnaSequence upper-cases input.
9. **COULD Test (invariant):** 0 ≤ H ≤ log₂(L−k+1) for any valid input. — Rationale: Shannon bounds.

---

## References

1. Li, H. (2025). Finding low-complexity DNA sequences with longdust. arXiv:2509.07357. https://arxiv.org/pdf/2509.07357
2. Çakır, et al. (2025). Entropy–Rank Ratio: A Novel Entropy-Based Perspective for DNA Complexity and Classification. arXiv:2511.05300. https://arxiv.org/html/2511.05300
3. Shannon, C. E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal 27:379–423, 623–656 — as exposited at https://tcosmo.github.io/2019/04/21/shannon-entropy.html and https://en.wikipedia.org/wiki/Entropy_(information_theory) (Shannon 1948 primary not directly machine-readable; bounds taken from these citing secondaries).

---

## Change History

- **2026-06-14**: Initial documentation.
