# Evidence Artifact: SEQ-COMPLEX-KMER-001

**Test Unit ID:** SEQ-COMPLEX-KMER-001
**Algorithm:** K-mer Entropy (Shannon entropy of the overlapping k-mer frequency distribution)
**Date Collected:** 2026-06-14

---

## Review 2026-09 correction (B04)

The original collection attributed the k-mer Shannon-entropy formula to Li (2025, longdust). Re-opening
the longdust primary material shows this is **wrong**: the longdust README ("The longdust algorithm",
raw.githubusercontent.com/lh3/longdust/master/README.md) and its math notes
(raw.githubusercontent.com/lh3/longdust/master/tex/notes.tex) define the score
S_L(x) = Σ_{t∈κ(x)} log c_x(t)! − f(ℓ(x)/4^k), a Poisson composite-likelihood score, **not**
H = −Σ p_i log₂ p_i. Li (2025) is kept only for the notation ℓ(x) = |x| − k + 1 (number of overlapping
k-mers), which the README does state. The "quoted" Li 2025 points 2–4 below were WebFetch summariser
output, not text of the paper, and must not be cited. Replacement sources actually opened / used:

| Source | What was opened | What it confirms |
|---|---|---|
| Shannon (1948) | via textbook formula (primary PDF not reachable) | H = −Σ p log p; 0 ≤ H ≤ log n; H = 0 iff deterministic |
| Herzel, Ebeling & Schmitt (1994) Phys. Rev. E 50:5061; Schmitt & Herzel (1997) J. Theor. Biol. 188:369 | WebSearch snippets only (publisher/arXiv blocked) | block entropy H_n = −Σ p⁽ⁿ⁾(A₁…A_n) log p⁽ⁿ⁾(A₁…A_n) of n-mers of DNA; finite-sample underestimation when N is small vs 4ⁿ |
| BBMap/BBDuk `EntropyTracker.java` (raw.githubusercontent.com/BioInfoTools/BBMap/master/current/structures/EntropyTracker.java) | source code opened | reference implementation: windowKmers = windowBases − k + 1 (overlapping), pk = count/windowKmers, eSum = Σ −pk·log(pk); BBDuk then multiplies by 1/ln(windowKmers) to get a 0–1 score |
| scipy.stats.entropy(counts, base=2) | executed (python3) | numerical cross-check of every dataset below plus a seeded random 200-mer, k ∈ {1,2,3,5,8} (see tests R1) |

Pastore et al. (2025) could not be re-opened (arXiv blocked); its log-base-2 convention is standard and
not load-bearing.

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

### Pastore et al. (2025). Entropy–Rank Ratio: An Entropy-Based Perspective for DNA Complexity (arXiv:2511.05300)

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

1. **MUST Test:** `ACGT`,k=1 → 2.0 (uniform, H = log₂4). — Evidence: Pastore et al. 2025 (max entropy = log₂4); Shannon uniform bound.
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
2. Pastore, E. P., Passarino, G., Sapia, P., De Rango, F. (2025). Entropy–Rank Ratio: A Novel Entropy-Based Perspective for DNA Complexity and Classification. arXiv:2511.05300. https://arxiv.org/html/2511.05300
3. Shannon, C. E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal 27:379–423, 623–656 — as exposited at https://tcosmo.github.io/2019/04/21/shannon-entropy.html and https://en.wikipedia.org/wiki/Entropy_(information_theory) (Shannon 1948 primary not directly machine-readable; bounds taken from these citing secondaries).

---

## Change History

- **2026-06-14**: Initial documentation.
- **2026-09-28**: Review 2026-09 (B04): corrected the Li 2025 misattribution (see top); added BBDuk
  EntropyTracker reference implementation and scipy cross-check values (ATGCATGCAT k=2 → 1.974937501201927;
  ATGCGATCGATCG k=2 → 2.4591479170272446, k=3 → 2.7321588913645702; seeded random 200-mer
  k=1/2/3/5/8 → 1.9964735194730474 / 3.937571048725419 / 5.6801547658649625 / 7.382517114501296 /
  7.582094342967564).
- **2026-10-01**: completeness audit WP13 (B04 F54): bias corrections and normalisation — see the section below.

## Revision 2026-10-01 — bias-corrected and normalised k-mer entropy (B04 F54)

**Sources opened:** R package `entropy` 1.3.2 (`raw.githubusercontent.com/cran/entropy/master/R/entropy.MillerMadow.R`,
`entropy.empirical.R`, `entropy.plugin.R`: `H = entropy.empirical(y, "log") + (m − 1)/(2n)`, `m = sum(y > 0)`, then
÷ log 2); `ndd` 1.10.6 sdist (`ndd/estimators.py`: `MillerMadow` = plug-in + 0.5(k − 1)/n; `Grassberger` docstring
"equation 35 in arXiv physics/0307138", G series G(1) = −γ − ln 2, G(2) = 2 + G(1), G(2m) = ψ(m + ½) + ln 2, odd = previous
even; its `fit` returns `log(n) − estimate/n` with `estimate = −Σ x G(x)`, i.e. ln N **+** Σ n G/N — a sign error that its own
`check.py` locks as 6.221 for counts whose plug-in entropy is 2.635; eq. 35 gives 2.734); infomeasure 0.6.3
(`estimators/entropy/grassberger.py` = Grassberger **1988**: ln N − ψ(n) − (−1)ⁿ/(n + 1); `miller_madow.py`); entropart
(`R/Shannon.R`: `Grassberger` = 1988 form, `Grassberger2003` = ψ(N) − Σ n/N (ψ(n) + (−1)ⁿ ∫₀¹ t^{n−1}/(1+t) dt)); WebSearch
snippets of arXiv:2310.07547 / Entropy 24:680 quoting Ĥ^G = ln N − (1/N) Σ n_i G_{n_i} with G₁ = −γ − ln 2, G₂ = 2 − γ − ln 2,
G_{2n+1} = G_{2n}, G_{2n+2} = G_{2n} + 2/(2n+1). arXiv / MDPI themselves are blocked (proxy 403 / EGRESS_BLOCKED). BBMap 40.02
`tracker/EntropyTracker.java` static `calcEntropy(bytes, counts, k)` / `calcEntropyFromCounts` (multiplier 1/ln(windowKmers)).

**Cross-check (harness `scratchpad/wp13/x_ent.py`, C# harness `xc ent`):** 3 000 strings (seed 1313; 60 % random ACGT,
25 % low-complexity, 15 % `ACGTNacgtRY`; L 1–3 000; k 1–10), 2 832 with N ≥ 1, 168 with L < k (all outputs 0):
R `entropy.empirical` / `entropy.MillerMadow` (log2) → 0 mismatches (max |Δ| 9.1e-13); Grassberger vs mpmath 40-digit eq. 35,
vs the recurrence (mpmath) and vs ndd's G series with eq. 35's sign → 0 (max 6.1e-14 / 6.1e-14 / 5.3e-15); normalised
(÷ log₂ N) plug-in / MM / Grassberger vs the oracles → 0 (2 801 cases with N > 1; 31 cases N = 1 → 0); normalised plug-in vs
BBTools `EntropyTracker.calcEntropy` (Java harness on `bbtools.jar`, 2 412 ACGT cases with N ≥ 2) → 0 float differences
(max |Δ| 5.8e-8 = float rounding). Digamma: `StatisticsHelper.Digamma` = mpmath/scipy on 9 points (≤ 2e-15 rel).

## WP15 — entropy-rank ratio R (B04 completeness audit, 2026-10-01): BLOCKED

- **Attribution corrected.** The arXiv listing for 2511.05300 (shown in a WebSearch result) names the authors as
  E. P. Pastore, G. Passarino, P. Sapia and F. De Rango (University of Calabria). The earlier "Çakır et al." was wrong.
- **What is known**, from search snippets of the abstract and §Methods: R is the share of all length-T blocks
  (under non-overlapping n-tuples) whose Shannon entropy is ≤ the target's. It is computed by turning frequency
  vectors into integer partitions and weighting them with multinomial coefficients. R lies in [0, 1] for fixed (T, n).
- **What is missing:**
  - how a remainder T mod n is handled;
  - the tie tolerance;
  - how non-ACGT symbols are treated;
  - any published R value (table or worked example);
  - the authors' code, which no search found.
- **URLs tried:**

  | Tool | Target | Result |
  |---|---|---|
  | `curl` | arxiv.org/abs/2511.05300, arxiv.org/html/2511.05300 and its v1, arxiv.org/pdf/2511.05300, export.arxiv.org/abs/2511.05300, api.semanticscholar.org/graph/v1/paper/arXiv:2511.05300 | each `curl: (56) CONNECT tunnel failed, response 403` |
  | `curl` | export.arxiv.org/api/query?id_list=2511.05300 | HTTP 403 `Host not in allowlist: export.arxiv.org` |
  | WebFetch | arxiv.org/html/2511.05300, www.alphaxiv.org/abs/2511.05300, api.semanticscholar.org | `EGRESS_BLOCKED` |
  | `gh search repos` / `gh search code` | "entropy rank ratio", "ratio-guided cropping" | HTTP 403 (session bound to its repositories) |
  | WebSearch | 4 queries | abstract and definition snippets only |

- **Decision:** not implemented. The definition's edge conventions cannot be confirmed, and there is no reference value
  to cross-check against (campaign rule 3). The spec line in K-mer_Entropy.md now records this.
