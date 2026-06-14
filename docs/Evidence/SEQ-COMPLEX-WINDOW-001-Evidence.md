# Evidence Artifact: SEQ-COMPLEX-WINDOW-001

**Test Unit ID:** SEQ-COMPLEX-WINDOW-001
**Algorithm:** Windowed Sequence Complexity (sliding-window complexity profile)
**Date Collected:** 2026-06-14

---

## Online Sources

### Linguistic sequence complexity (Wikipedia, citing Trifonov 1990; Gabrielian & Bolshoy 1999; Troyanskaya et al. 2002)

**URL:** https://en.wikipedia.org/wiki/Linguistic_sequence_complexity
**Accessed:** 2026-06-14
**Retrieved by:** WebSearch query "linguistic complexity sequence Gabrielian Bolshoy sliding window vocabulary usage"; then WebFetch of the article URL above.
**Authority rank:** 4 (Wikipedia citing primaries; the primaries Trifonov 1990 / Troyanskaya 2002 are ranked 1).

**Key Extracted Points:**

1. **Vocabulary usage Uᵢ:** "the ratio of the actual vocabulary size of a given sequence to the maximal possible vocabulary size for a sequence of that length." (extracted from fetched article)
2. **Maximum possible vocabulary at level i:** `min(4^i, N - i + 1)` where `4^i` is the DNA-alphabet theoretical maximum and `N - i + 1` is the positional maximum for a sequence of length N. (extracted from fetched article)
3. **Range:** "The value of C provides a measure of sequence complexity in the range 0<C<1 for various DNA sequence fragments." (extracted from fetched article)
4. **Sliding-window use:** repo's linguistic-complexity unit (SEQ-COMPLEX-001) uses the *summation* variant `LC = (Σ Vᵢ) / (Σ Vmax,i)`; this windowed unit applies that exact per-window LC plus Shannon entropy across a sliding window to produce a complexity profile. (repo doc `docs/algorithms/Sequence_Composition/Linguistic_Complexity.md`, §2.2)

### Linguistic complexity profiles — fast algorithm (Troyanskaya et al. 2002)

**URL:** https://pubmed.ncbi.nlm.nih.gov/12050064/
**Accessed:** 2026-06-14
**Retrieved by:** WebSearch (same query as above) surfaced the PubMed entry; then WebFetch of the URL.
**Authority rank:** 1 (peer-reviewed, Bioinformatics 18(5):679–688).

**Key Extracted Points:**

1. **Subword-count definition:** "Linguistic complexity corresponds to repetitiveness of a genomic text"; the program "utilizes suffix trees to compute the number of subwords present in genomic sequences, thereby allowing calculation of linguistic complexity in time linear in genome size." (extracted from fetched abstract) — i.e. complexity is driven by counts of distinct subwords, the summation form used by the repo's LC method.
2. **Profile concept:** the paper computes *complexity profiles* (complexity as a function of position along a genome), which is precisely the sliding-window output of this unit. (title: "Sequence complexity profiles of prokaryotic genomic sequences")

### Entropy (information theory) — Shannon entropy (Wikipedia, citing Shannon 1948)

**URL:** https://en.wikipedia.org/wiki/Entropy_(information_theory)
**Accessed:** 2026-06-14
**Retrieved by:** WebSearch query "Shannon entropy formula H = -sum p log2 p information content bits"; then WebFetch of the article URL.
**Authority rank:** 4 (Wikipedia citing the primary Shannon 1948, ranked 1).

**Key Extracted Points:**

1. **Definition:** `H(X) = -Σ p(x) log_b p(x)`; base 2 yields units of bits. (extracted from fetched article)
2. **Maximum:** "Entropy reaches its maximum value of log_b(n) when the distribution is uniform across n equally likely outcomes." For DNA (n=4 bases) max = log₂(4) = 2 bits. (extracted from fetched article)
3. **Minimum:** "For a deterministic distribution (certainty), entropy equals zero." (extracted from fetched article)
4. **Zero convention:** "the value of the corresponding summand 0 log_b(0) is taken to be 0." (extracted from fetched article)
5. **Primary:** Claude Shannon, "A Mathematical Theory of Communication," Bell System Technical Journal, 1948. (extracted from fetched article)

---

## Documented Corner Cases and Failure Modes

### From Wikipedia (Entropy) / Shannon 1948

1. **Uniform window:** every base equally frequent ⇒ Shannon entropy = log₂(4) = 2.0 bits (maximum for DNA).
2. **Homopolymer window:** a single base ⇒ deterministic distribution ⇒ Shannon entropy = 0.

### From Wikipedia (Linguistic sequence complexity) / Troyanskaya 2002

3. **Highly repetitive window:** low linguistic complexity because few distinct subwords are observed relative to the maximum possible.
4. **Window shorter than full sequence near the 3′ end:** the sliding driver only emits windows fully contained in the sequence (`i + w ≤ L`); a trailing partial fragment is not emitted.

---

## Test Datasets

### Dataset: Uniform repeat window (Shannon max + LC)

**Source:** Shannon 1948 (entropy max for uniform distribution); Linguistic sequence complexity (Wikipedia), Uᵢ = observed/Vmax with Vmax,i = min(4^i, N-i+1).

Window = `ACGTACGT` (length 8, equal base counts A=C=G=T=2), LC computed with maxWordLength = min(6, windowSize) = 6:

| Length i | distinct observed Vᵢ | Vmax,i = min(4^i, 8-i+1) |
|----------|----------------------|---------------------------|
| 1 | 4 (A,C,G,T) | 4 |
| 2 | 4 (AC,CG,GT,TA) | 7 |
| 3 | 4 (ACG,CGT,GTA,TAC) | 6 |
| 4 | 4 (ACGT,CGTA,GTAC,TACG) | 5 |
| 5 | 4 (ACGTA,CGTAC,GTACG,TACGT) | 4 |
| 6 | 3 (ACGTAC,CGTACG,GTACGT) | 3 |

| Parameter | Value |
|-----------|-------|
| ShannonEntropy | 2.0 (uniform, log₂4) |
| Σ Vᵢ (observed) | 23 |
| Σ Vmax,i (possible) | 29 |
| LinguisticComplexity | 23/29 = 0.7931034482758621 |

### Dataset: Homopolymer window (poly-A)

**Source:** Shannon 1948 (deterministic ⇒ 0); Linguistic sequence complexity (Wikipedia).

Window = `AAAAAAAA` (length 8), LC with maxWordLength = 6:

| Length i | distinct observed Vᵢ | Vmax,i = min(4^i, 8-i+1) |
|----------|----------------------|---------------------------|
| 1 | 1 | 4 |
| 2 | 1 | 7 |
| 3 | 1 | 6 |
| 4 | 1 | 5 |
| 5 | 1 | 4 |
| 6 | 1 | 3 |

| Parameter | Value |
|-----------|-------|
| ShannonEntropy | 0.0 (single base, deterministic) |
| Σ Vᵢ (observed) | 6 |
| Σ Vmax,i (possible) | 29 |
| LinguisticComplexity | 6/29 = 0.20689655172413793 |

### Dataset: Window enumeration geometry

**Source:** Sliding-window definition (Troyanskaya 2002 complexity profile) + repository contract.

| Parameter | Value |
|-----------|-------|
| Sequence length L | 24 |
| windowSize w | 8 |
| stepSize s | 8 |
| Number of windows | floor((24-8)/8)+1 = 3 |
| Window 0 | WindowStart=0, WindowEnd=7, Position=0+8/2=4 |
| Window 1 | WindowStart=8, WindowEnd=15, Position=8+4=12 |
| Window 2 | WindowStart=16, WindowEnd=23, Position=16+4=20 |

---

## Assumptions

1. **ASSUMPTION: Center-position convention** — The reported `Position` is the window center computed as `WindowStart + windowSize/2` using integer division. No cited source mandates the center label (sources specify the *profile values*, not the position label); the center convention is the repository's choice and is non-correctness-affecting for the complexity values themselves. It is the documented contract of `ComplexityPoint`.
2. **ASSUMPTION: Default windowSize=64, stepSize=10** — The repository defaults follow common low-complexity windowing practice but are not fixed by a single cited source; they are caller-overridable parameters and do not affect the value computed for an explicitly-specified window.
3. **ASSUMPTION: per-window LC uses maxWordLength = min(6, windowSize)** — The cap of 6 mirrors Gabrielian & Bolshoy's "limit vocabulary assessment to W rather than all N-1 values" efficiency choice; the specific cap (6) is a repository parameter consistent with the existing LC unit.

---

## Recommendations for Test Coverage

1. **MUST Test:** Window count = floor((L-w)/s)+1 for L≥w, and 0 when L<w. — Evidence: sliding-window geometry (Troyanskaya 2002 profile; repository contract).
2. **MUST Test:** ShannonEntropy of a uniform window = 2.0 and of a homopolymer window = 0.0. — Evidence: Shannon 1948 (uniform max = log₂4; deterministic = 0).
3. **MUST Test:** LinguisticComplexity of `ACGTACGT` window = 23/29 and of poly-A window = 6/29. — Evidence: Wikipedia LC (Uᵢ = observed/min(4^i,N-i+1)) hand-derivation above.
4. **MUST Test:** WindowStart/WindowEnd/Position coordinates per window (0-based, end inclusive, center = start+w/2). — Evidence: repository ComplexityPoint contract.
5. **SHOULD Test:** Null DnaSequence ⇒ ArgumentNullException; windowSize<1 and stepSize<1 ⇒ ArgumentOutOfRangeException. — Rationale: documented failure modes.
6. **SHOULD Test:** L < w ⇒ empty profile (no partial trailing window). — Rationale: corner case 4.
7. **COULD Test:** Invariant 0 ≤ ShannonEntropy ≤ 2 and 0 ≤ LinguisticComplexity ≤ 1 across windows. — Rationale: Shannon/LC value bounds.

---

## References

1. Shannon, C.E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal, 27(3):379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x (formula/properties accessed via https://en.wikipedia.org/wiki/Entropy_(information_theory), 2026-06-14)
2. Troyanskaya, O.G., Arbell, O., Koren, Y., Landau, G.M., Bolshoy, A. (2002). Sequence complexity profiles of prokaryotic genomic sequences: a fast algorithm for calculating linguistic complexity. Bioinformatics, 18(5):679–688. https://doi.org/10.1093/bioinformatics/18.5.679 (https://pubmed.ncbi.nlm.nih.gov/12050064/, accessed 2026-06-14)
3. Gabrielian, A., Bolshoy, A. (1999). Sequence complexity and DNA curvature. Computers & Chemistry, 23(3–4):263–274. https://doi.org/10.1016/S0097-8485(99)00007-8 (cited via https://en.wikipedia.org/wiki/Linguistic_sequence_complexity, accessed 2026-06-14)
4. Trifonov, E.N. (1990). Making sense of the human genome. In: Structure & Methods, Vol 1, Adenine Press. (cited via https://en.wikipedia.org/wiki/Linguistic_sequence_complexity, accessed 2026-06-14)
5. Wikipedia. Linguistic sequence complexity. https://en.wikipedia.org/wiki/Linguistic_sequence_complexity (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
