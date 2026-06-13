# Evidence Artifact: PROTMOTIF-LC-001

**Test Unit ID:** PROTMOTIF-LC-001
**Algorithm:** Low-Complexity Region Detection (SEG, Wootton & Federhen 1993)
**Date Collected:** 2026-06-14

---

## Online Sources

### NCBI `ncbi-seg` man page (official program specification)

**URL:** https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 2 (official specification / standard distribution of the SEG program)

**Key Extracted Points:**

1. **Trigger window length [W]:** "An integer greater than zero [Default 12]." — the sliding window used to evaluate local complexity.
2. **Trigger complexity [K1]:** "The maximum complexity of a trigger window in units of bits. K1 must be equal to or greater than zero. The maximum value is 4.322 (log[base 2]20) for amino acid sequences [Default 2.2]." — windows with complexity ≤ K1 trigger a low-complexity segment.
3. **Extension complexity [K2]:** "The maximum complexity of an extension window in units of bits. Only values greater than K1 are effective in extending triggered windows. Range of possible values is as for K1 [Default 2.5]."
4. **Complexity definition:** "Complexity" is "defined by equation (3) of Wootton & Federhen (1993)"; the permissible range is 0 to 4.322 bits per residue, i.e. 0 to log₂(20) (maximum entropy for 20 amino acids).

### NCBI BLAST `blast_seg.c` (reference implementation)

**URL:** https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html and …/blast__seg_8c_source.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation in an established library, NCBI C++ Toolkit)

**Key Extracted Points:**

1. **Defaults (verbatim constants):** `const int kSegWindow = 12; const double kSegLocut = 2.2; const double kSegHicut = 2.5;` — confirms W=12, K1(locut)=2.2, K2(hicut)=2.5.
2. **Entropy function:** `double s_Entropy(Int4 *sv)` "Calculates entropy of an integer array" by computing Shannon entropy across residue frequency counts, normalized by the logarithm of alphabet size. The `s_LnPerm()` helper "calculates 'K2' entropy per Wootton and Federhen's equation 3."
3. **Log-factorial table:** the source contains a precomputed `lnfact[]` array (natural logarithms of factorials), used by the exact compositional-complexity / permutation form of equation (3).
4. **Two-stage algorithm:** "first, identification of approximate raw segments of low-complexity" (all overlapping trigger windows of length W with complexity ≤ K1), "second local optimization" extending while complexity ≤ K2.

### SeqComplex reference implementation (caballero, `SeqComplex.pm`)

**URL:** https://raw.githubusercontent.com/caballero/SeqComplex/master/SeqComplex.pm (retrieved with `curl`)
**Accessed:** 2026-06-14
**Authority rank:** 3 (open reference implementation)

**Key Extracted Points:**

1. **Shannon-entropy complexity (`sub ce`), verbatim formula in the doc comment:** `H(X) = SUM( P(Xi) * Log2( P(Xi) ) )` over symbols i; "The unit of entropy for log2 is shannon or more commonly as 'bits'." Code: for each symbol with count > 0, `r = count/tot; ce -= r * log_k(2, r)` → K = −Σ pᵢ·log₂(pᵢ).
2. **Wootton-Federhen complexity (`sub cwf`):** `up = Σ_{1..W} log_k(N, W)`; `dw = Σ_{distinct b} log_k(N, n_b)`; `complexity = (up - dw)/tot`. This is the Stirling-style per-residue form of K = (1/L)·log_N(L^L / Πnᵢ) using alphabet base N.
3. **`log_k(base, num)`:** `res = log(num)/log(base)` — change-of-base logarithm.

### universalmotif `sequence_complexity` documentation (Bioconductor R package)

**URL:** https://rdrr.io/github/bjmt/universalmotif/man/sequence_complexity.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 3 (reference implementation documentation), citing Wootton & Federhen 1993 and Orlov & Potapov 2004 as primaries.

**Key Extracted Points:**

1. **WF interpretation:** "the Wootton-Federhen complexity score is a reflection of the numbers of each unique letter found in the window"; less complex sequences score closer to 0, more complex closer to 1 (normalized variant).
2. **Provenance:** "The Wootton-Federhen (Wootton and Federhen, 1993) and Trifonov (Trifonov, 1990) algorithms … are well described within Orlov and Potapov (2004)."

### Pei & Grishin, *Bioinformatics* 21(2):160 — "A new algorithm for detecting low-complexity regions"

**URL:** https://academic.oup.com/bioinformatics/article/21/2/160/187330
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **SEG mechanism:** "SEG detects LCRs based on an information measure of the complexity state vector, which reflects residue composition appearing on a sliding window, with no regard of the patterns or periodicity of sequence repetitiveness."
2. **Two-pass + defaults:** "a two-pass algorithm" using "a sliding window" with "default parameters W = 12, K2(1) = 2.2 bit and K2(2) = 2.5 bits."

### Shashidhara/Mier et al., *Briefings in Bioinformatics* 21(2):458 — Shannon entropy form

**URL:** https://academic.oup.com/bioinformatics/article/22/24/2980/208627
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Shannon entropy complexity:** "One of the most popularly used complexity measures for sequences, Shannon entropy (Shannon, 1951), is defined as −∑ᵢ₌₁²⁰ pᵢ log pᵢ" where pᵢ is the fractional composition of amino acid type i and the sum is over the 20 amino acids.

---

## Documented Corner Cases and Failure Modes

### From NCBI `ncbi-seg` man page

1. **Window longer than sequence:** SEG operates on windows of length W; sequences shorter than W contain no complete trigger window, so no low-complexity segment can be triggered.
2. **K1 ≥ K2 effectiveness:** "Only values greater than K1 are effective in extending triggered windows" — extension complexity K2 must exceed K1 to have an effect.
3. **Parameter bounds:** complexity cutoffs must be in [0, log₂(N)]; for 20 amino acids the maximum meaningful cutoff is 4.322.

### From SeqComplex (`cwf`/`ce`)

1. **Single-symbol denominator guard:** entropy/complexity is defined only when total count `tot > 1`; for a degenerate window the contribution is 0 (a homopolymer has complexity 0).

---

## Test Datasets

### Dataset: SEG default parameters (official spec)

**Source:** NCBI `ncbi-seg` man page; NCBI `blast_seg.c` (`kSegWindow`, `kSegLocut`, `kSegHicut`).

| Parameter | Value |
|-----------|-------|
| Window length W | 12 |
| Trigger complexity K1 (locut) | 2.2 bits/residue |
| Extension complexity K2 (hicut) | 2.5 bits/residue |
| Max complexity (20 aa) | log₂(20) = 4.321928 |

### Dataset: Worked complexity values (derived from Shannon-entropy formula K = −Σ pᵢ·log₂ pᵢ)

**Source:** Formula from Mier et al. / SeqComplex `ce`; values computed independently (Python `math.log2`).

| Window composition (L=12) | K (bits/residue) |
|---------------------------|------------------|
| 12 × A (homopolymer) | 0.000000 |
| 11 × A, 1 × B | 0.413817 |
| 10 × A, 2 × B | 0.650022 |
| 8 × A, 4 × B | 0.918296 |
| 6 × A, 6 × B | 1.000000 |
| 12 distinct residues (1 each) | log₂(12) = 3.584963 |

---

## Assumptions

1. **ASSUMPTION: Complexity measure choice (Shannon entropy in bits/residue).** Equation (3) of Wootton & Federhen has two interconvertible operational forms in the retrieved reference implementations: the exact compositional complexity K = (1/L)·log_N(L!/Πnᵢ!) (NCBI `lnfact[]` permutation form) and the Shannon entropy −Σ pᵢ·log₂ pᵢ (NCBI `s_Entropy`, SeqComplex `ce`, Mier et al.). The official man page specifies complexity "in units of bits … maximum 4.322 = log₂20," which is the Shannon-entropy form in bits per residue. This implementation uses the Shannon-entropy bits/residue form, matching the man-page units and the `s_Entropy`/`ce` reference code. Justification: it is the form whose range (0 … log₂20) and units ("bits") are stated verbatim in the official specification.
2. **ASSUMPTION: Empty / short-sequence behavior.** The spec defines complexity only for full windows of length W; the original sources do not prescribe a return value for an input shorter than W or empty. This implementation returns no regions (empty result) for sequences shorter than the window, consistent with "no complete trigger window exists."

---

## Recommendations for Test Coverage

1. **MUST Test:** Homopolymer window has complexity 0 and is reported as a low-complexity region. — Evidence: K = −Σ pᵢ log₂ pᵢ = 0 for a single symbol (SeqComplex `ce`, Mier et al.).
2. **MUST Test:** A maximally diverse window (12 distinct residues, L=12) has complexity log₂(12) = 3.584963 > K2 = 2.5 and is NOT a low-complexity region. — Evidence: formula + SEG defaults.
3. **MUST Test:** Per-window complexity equals the exact Shannon entropy for biased windows (11A/1B → 0.413817; 6A6B → 1.0; 10A/2B → 0.650022). — Evidence: derived values table.
4. **MUST Test:** SEG defaults are W=12, K1=2.2, K2=2.5. — Evidence: man page + `blast_seg.c` constants.
5. **MUST Test:** A poly-Q tract embedded in a diverse protein is detected as a single low-complexity region with correct inclusive boundaries. — Evidence: SEG detects compositionally biased segments; boundary = window span.
6. **SHOULD Test:** Sequence shorter than the window yields no regions. — Rationale: no complete trigger window (man page).
7. **SHOULD Test:** Two separated biased tracts yield two distinct regions. — Rationale: independent triggered segments.
8. **COULD Test:** Region complexity is in [0, log₂20]. — Rationale: value-bound invariant.

---

## References

1. Wootton JC, Federhen S. 1993. Statistics of local complexity in amino acid sequences and sequence databases. Computers & Chemistry 17(2):149–163. https://doi.org/10.1016/0097-8485(93)85006-X
2. NCBI. SEG program man page (`ncbi-seg`). https://manpages.ubuntu.com/manpages/focal/man1/ncbi-seg.1.html
3. NCBI C++ Toolkit. `blast_seg.c` (s_Entropy, s_LnPerm, kSegWindow/kSegLocut/kSegHicut). https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/doxyhtml/blast__seg_8c.html
4. Caballero J. SeqComplex (`SeqComplex.pm`, subs `cwf`, `ce`, `log_k`). https://github.com/caballero/SeqComplex
5. Tremblay BJM. universalmotif: `sequence_complexity`. Bioconductor. https://rdrr.io/github/bjmt/universalmotif/man/sequence_complexity.html
6. Pei J, Grishin NV. 2005. A new algorithm for detecting low-complexity regions in protein sequences (CARD/SEG comparison). Bioinformatics 21(2):160–166. https://academic.oup.com/bioinformatics/article/21/2/160/187330
7. Shannon CE. 1951. Prediction and entropy of printed English. Bell System Technical Journal. (Shannon entropy −Σ pᵢ log pᵢ as cited by Mier et al., Bioinformatics 22(24):2980.) https://academic.oup.com/bioinformatics/article/22/24/2980/208627

---

## Change History

- **2026-06-14**: Initial documentation.
