# Evidence Artifact: SEQ-ATSKEW-001

**Test Unit ID:** SEQ-ATSKEW-001
**Algorithm:** AT Skew — (A − T) / (A + T)
**Date Collected:** 2026-06-14

---

## Online Sources

### Lobry, J. R. (1996) — Asymmetric substitution patterns in the two DNA strands of bacteria (primary source)

**URL:** https://pubmed.ncbi.nlm.nih.gov/8676740/ (DOI: 10.1093/oxfordjournals.molbev.a025626)
**Accessed:** 2026-06-14 (retrieved via WebFetch of the PubMed record)
**Authority rank:** 1 (peer-reviewed primary paper)

**Key Extracted Points:**

1. **Origin of the skew concept:** PubMed abstract retrieved verbatim — "There was a departure from intrastrand equifrequency between A and T or between C and G, showing that the substitution patterns of the two strands of DNA were asymmetric." This is the founding observation that AT skew quantifies (deviation of A from T within one strand).
2. **Citation:** J R Lobry. Mol Biol Evol. 1996 May;13(5):660–5. PMID 8676740, DOI 10.1093/oxfordjournals.molbev.a025626. Analyzed E. coli, B. subtilis, H. influenzae.

### Charneski et al. (2011) — Atypical AT Skew in Firmicute Genomes (PLOS Genetics)

**URL:** https://journals.plos.org/plosgenetics/article?id=10.1371%2Fjournal.pgen.1002283
**Accessed:** 2026-06-14 (retrieved via WebFetch of the article HTML)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Formula (verbatim):** AT skew is defined as **(A − T) / (A + T)**. Stated in the abstract and reused throughout when discussing leading-strand nucleotide asymmetry.
2. **Citation:** Charneski CA, Honti F, Bryant JM, Hurst LD, Feil EJ (2011). PLoS Genet 7(9): e1002283.

### Wikipedia "GC skew" (cited primary: Lobry 1996)

**URL:** https://en.wikipedia.org/wiki/GC_skew
**Accessed:** 2026-06-14 (retrieved via WebFetch)
**Authority rank:** 4 (Wikipedia citing the Lobry 1996 primary, used here only to corroborate formula + range)

**Key Extracted Points:**

1. **Formulas (verbatim):** `GC skew = (G − C)/(G + C)` and `AT skew = (A − T)/(A + T)`.
2. **Value range (verbatim):** "The nucleotide composition skew spectra ranges from −1 ... to +1". For AT skew, −1 corresponds to A = 0 and +1 corresponds to T = 0.
3. **Primary attribution:** the definitions trace to Lobry, J. R. (1996), Mol Biol Evol 13(5):660–665.

### Biopython `Bio.SeqUtils.GC_skew` source (reference implementation)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-14 (retrieved via WebFetch of the raw source file)
**Authority rank:** 3 (reference implementation in an established library)

**Key Extracted Points:**

1. **Counting (verbatim from source):** `g = s.count("G") + s.count("g")` and `c = s.count("C") + s.count("c")` — base counting is **case-insensitive**. By analogy the AT-skew of the library convention counts A/a and T/t.
2. **Zero-denominator handling (verbatim):** `try: skew = (g - c) / (g + c) except ZeroDivisionError: skew = 0.0` — when the relevant base pair is absent (denominator 0) the skew is **0.0**, not NaN/exception.
3. **Ambiguous bases (verbatim docstring):** "Does NOT look at any ambiguous nucleotides." Only the canonical bases are counted; all other symbols are ignored in the numerator and denominator.

---

## Documented Corner Cases and Failure Modes

### From Biopython `GC_skew`

1. **Empty window / no A or T:** denominator A + T = 0 ⇒ skew = 0.0 (ZeroDivisionError caught and replaced by 0.0). The AT-skew analog returns 0 when the sequence contains no A and no T.
2. **Non-ACGT symbols ignored:** characters that are not the counted bases (gaps, N, IUPAC ambiguity) do not contribute to numerator or denominator.

### From Wikipedia / Lobry (1996)

1. **Bounds:** skew ∈ [−1, +1]. AT skew = +1 ⇔ T = 0 (all A among A/T); AT skew = −1 ⇔ A = 0 (all T among A/T).

---

## Test Datasets

### Dataset: Hand-derived worked values from the formula (A − T)/(A + T)

**Source:** Formula from Charneski et al. (2011) and Wikipedia/Lobry (1996); values computed by direct substitution (no library run needed — they are arithmetic consequences of the definition).

| Input sequence | A | T | (A − T)/(A + T) | Expected AT skew |
|----------------|---|---|------------------|------------------|
| `AAAA` | 4 | 0 | 4/4 | 1.0 |
| `TTTT` | 0 | 4 | −4/4 | −1.0 |
| `ATAT` | 2 | 2 | 0/4 | 0.0 |
| `AAAT` | 3 | 1 | 2/4 | 0.5 |
| `ATTT` | 1 | 3 | −2/4 | −0.5 |
| `GGCC` | 0 | 0 | 0/0 → 0 | 0.0 (no A/T) |
| `AAATGGGCCC` (A=3, T=1, G/C ignored) | 3 | 1 | 2/4 | 0.5 |
| `aaat` (lowercase) | 3 | 1 | 2/4 | 0.5 (case-insensitive) |

---

## Assumptions

1. **ASSUMPTION: Lowercase + non-ACGT handling for the AT-skew analog.** AT skew is not given its own source-code line in Biopython (only GC_skew is shipped), so the case-insensitive counting and "ignore everything that is not A/T" behavior are taken from the directly analogous `GC_skew` reference implementation rather than from an AT-skew-specific source. The formula itself is fully sourced (Charneski 2011, Lobry 1996); only the symbol-handling convention is inferred by analogy. This matches the repository implementation, which uppercases input and counts only 'A'/'T'.

---

## Recommendations for Test Coverage

1. **MUST Test:** AT skew of a pure-A sequence = +1.0 and pure-T = −1.0 (the bounds). — Evidence: Wikipedia/Lobry range [−1,+1].
2. **MUST Test:** balanced A = T ⇒ 0.0, and asymmetric A:T ⇒ exact fraction (e.g. AAAT ⇒ 0.5). — Evidence: formula (A−T)/(A+T), Charneski 2011.
3. **MUST Test:** sequence with no A and no T ⇒ 0.0 (zero-denominator). — Evidence: Biopython ZeroDivisionError → 0.0.
4. **MUST Test:** G/C and other non-A/T symbols are ignored (denominator counts only A+T). — Evidence: Biopython "does NOT look at ambiguous nucleotides" + A/T counting.
5. **SHOULD Test:** lowercase input gives the same value as uppercase. — Rationale: case-insensitive counting (Biopython); repo normalizes via ToUpperInvariant.
6. **SHOULD Test:** null string ⇒ 0; empty string ⇒ 0; null DnaSequence ⇒ ArgumentNullException. — Rationale: documented input validation / failure modes.
7. **COULD Test:** DnaSequence overload agrees with the string overload on the same sequence. — Rationale: delegate equivalence.

---

## References

1. Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. Molecular Biology and Evolution 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626 (PMID 8676740)
2. Charneski, C. A., Honti, F., Bryant, J. M., Hurst, L. D., Feil, E. J. (2011). Atypical AT Skew in Firmicute Genomes Results from Selection and Not from Mutation. PLoS Genetics 7(9): e1002283. https://doi.org/10.1371/journal.pgen.1002283
3. Wikipedia contributors. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-14).
4. Biopython project. Bio.SeqUtils (GC_skew). https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/__init__.py (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
