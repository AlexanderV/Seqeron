# Evidence Artifact: SEQ-COMPOSITION-001

**Test Unit ID:** SEQ-COMPOSITION-001
**Algorithm:** Sequence Composition (nucleotide composition + amino-acid composition)
**Date Collected:** 2026-06-14

---

## Online Sources

### Biopython — `Bio.SeqUtils` source (`gc_fraction`, `GC_skew`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-14 (fetched the raw source file in this session)
**Authority rank:** 3 (reference implementation in an established bioinformatics library)

**Key Extracted Points:**

1. **GC content (`gc_fraction`):** Core computation `gc = sum(seq.count(x) for x in "CGScgs")`; returns `gc / length`, a float in [0, 1]. Default `ambiguous="remove"` "will only count GCS and will only include ACTGSWU when calculating the sequence length."
2. **Empty sequence:** the function returns zero for an empty sequence (length-zero denominator handled).
3. **GC skew (`GC_skew`):** `skew = (g - c) / (g + c)` where `g`/`c` are the per-window G/C counts; returns a list of floats over sequential windows.
4. **Zero-denominator handling (`GC_skew`):** returns `0.0` for a window with no G or C (catches the zero-division case).
5. **Case-insensitivity:** GC counting explicitly includes lowercase letters (`"CGScgs"`), so composition is case-insensitive.

### Wikipedia — "GC skew" (using the cited primary, Lobry 1996)

**URL:** https://en.wikipedia.org/wiki/GC_skew
**Accessed:** 2026-06-14 (fetched in this session)
**Authority rank:** 4 (Wikipedia; the formula traces to the Lobry 1996 primary below)

**Key Extracted Points:**

1. **GC skew formula:** `GC skew = (G − C)/(G + C)`.
2. **AT skew formula:** `AT skew = (A − T)/(A + T)`.
3. **Primary citation:** Lobry, J. R. (1996), *Molecular Biology and Evolution* 13:660–665.
4. **Cumulative skew / replication origin:** in cumulative GC-skew plots the peaks correspond to the switch points (terminus / origin of replication).

### IUPAC codes — nucleotide and amino-acid single-letter symbols

**URL:** https://www.bioinformatics.org/sms/iupac.html
**Accessed:** 2026-06-14 (fetched in this session)
**Authority rank:** 2 (IUPAC nomenclature standard, reproduced)

**Key Extracted Points:**

1. **Canonical nucleotides:** A, C, G, T, U. Degenerate codes include R, Y, S, W, K, M, B, D, H, V and **N = any base**.
2. **20 standard amino-acid single-letter codes:** A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y.

---

## Documented Corner Cases and Failure Modes

### From Biopython `Bio.SeqUtils`

1. **Empty sequence:** `gc_fraction` returns 0 for an empty sequence.
2. **No G or C present (zero denominator):** `GC_skew` returns 0.0 when `g + c == 0`. By the same reasoning AT skew returns 0 when `a + t == 0`.
3. **Mixed case:** counting includes lowercase, so the result is case-insensitive.

### From IUPAC

1. **Non-canonical letters:** symbols outside {A,C,G,T,U} (e.g. `N`, degenerate codes, `X`) are not standard nucleotides and are tracked separately from the four/five canonical bases.

---

## Test Datasets

### Dataset: Hand-derived worked examples (from the formulas above)

**Source:** Biopython `Bio.SeqUtils` and Wikipedia "GC skew" (cited above).

| Input | A | T | G | C | U | GC content (G+C)/total | GC skew (G−C)/(G+C) | AT skew (A−T)/(A+T) |
|-------|---|---|---|---|---|------------------------|---------------------|---------------------|
| `ATGC` | 1 | 1 | 1 | 1 | 0 | 2/4 = 0.5 | 0/2 = 0 | 0/2 = 0 |
| `GGGC` | 0 | 0 | 3 | 1 | 0 | 4/4 = 1.0 | 2/4 = 0.5 | 0 (a+t=0) |
| `AAUUGGCC` | 2 | 0 | 2 | 2 | 2 | 4/8 = 0.5 | 0/4 = 0 | 2/2 = 1.0 |

### Dataset: Amino-acid composition worked example

**Source:** IUPAC single-letter amino-acid codes (cited above).

| Input | Residues | Length | Counts |
|-------|----------|--------|--------|
| `MKVLWA` | M,K,V,L,W,A | 6 | each = 1 |

---

## Assumptions

1. **ASSUMPTION: Degenerate IUPAC codes (S, W, R, Y, …) are not counted toward composition totals.** Biopython's `gc_fraction` counts `S` toward GC and `W` toward the denominator. The repository implementation counts only A/T/G/C/U toward GC/AT totals and routes other letters to `CountN`/`CountOther`. For sequences over the standard {A,T,G,C,U} alphabet (this unit's scope) the two agree exactly; the difference manifests only on degenerate symbols. Documented as an intentional simplification, not an invented constant.

---

## Recommendations for Test Coverage

1. **MUST Test:** Exact A/T/G/C/U/N/Other counts and Length partition — Evidence: IUPAC alphabet + definition of nucleotide composition.
2. **MUST Test:** GC content = (G+C)/(A+T+G+C+U) — Evidence: Biopython `gc_fraction`.
3. **MUST Test:** GC skew = (G−C)/(G+C) incl. a negative case — Evidence: Wikipedia "GC skew" / Biopython `GC_skew`.
4. **MUST Test:** AT skew = (A−T)/(A+T) — Evidence: Wikipedia "GC skew".
5. **MUST Test:** Empty/null → all-zero composition — Evidence: Biopython `gc_fraction` empty handling.
6. **MUST Test:** Amino-acid composition exact residue counts and length — Evidence: IUPAC amino-acid codes.
7. **SHOULD Test:** Case-insensitivity — Rationale: Biopython counts lowercase.
8. **SHOULD Test:** Zero-denominator skews → 0 — Rationale: Biopython zero-division handling.

---

## References

1. Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution* 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
2. Cock, P. J. A. et al. Biopython, `Bio.SeqUtils` (`gc_fraction`, `GC_skew`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py (accessed 2026-06-14)
3. Wikipedia contributors. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-14)
4. IUPAC nucleotide and amino-acid single-letter codes. https://www.bioinformatics.org/sms/iupac.html (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation. Records that SEQ-COMPOSITION-001 is a duplicate Registry entry for the two composition methods already delivered under SEQ-STATS-001; consolidated rather than re-implemented (see TestSpec §7).
