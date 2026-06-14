# Evidence Artifact: SEQ-GC-PROFILE-001

**Test Unit ID:** SEQ-GC-PROFILE-001
**Algorithm:** GC Content Profile (sliding-window GC content)
**Date Collected:** 2026-06-14

---

## Online Sources

### Wikipedia — GC-content (citing primary literature)

**URL:** https://en.wikipedia.org/wiki/GC-content
**Accessed:** 2026-06-14 (fetched in-session via WebFetch of the article URL)
**Authority rank:** 4 (Wikipedia citing primaries; used for the closed-form definition)

**Key Extracted Points:**

1. **GC-content percentage formula (verbatim):** the article gives GC-content as
   `(G + C) / (A + T + G + C) × 100%` — numerator is the count of guanine and
   cytosine bases; denominator is the count of the four standard bases A, T, G, C;
   the result is multiplied by 100 to express it as a percentage.
2. **Denominator definition:** the denominator is the total of the four standard
   bases (A+T+G+C), i.e. ambiguous/other symbols (e.g. N) are not part of the
   standard-base denominator in this closed form.
3. **AT/GC ratio (verbatim):** `(A + T) / (G + C)` — an alternative compositional ratio
   (not used by this unit, recorded for completeness).

### Biopython — Bio.SeqUtils.gc_fraction (reference implementation)

**URL:** https://biopython.org/docs/latest/api/Bio.SeqUtils.html ;
source file https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-14 (both the API doc page and the raw source file fetched in-session via WebFetch)
**Authority rank:** 3 (reference implementation; Biopython, `Bio.SeqUtils`, function `gc_fraction`)

**Key Extracted Points:**

1. **gc_fraction definition:** "Calculate G+C percentage in seq (float between 0 and 1)."
   `gc_fraction(seq, ambiguous="remove")`. Returns a fraction in [0, 1]; multiplying by
   100 yields the Wikipedia percentage.
2. **Numerator:** G and C are always counted as full GC contributors; S (ambiguity code
   for G or C) is also counted as 1.0.
3. **Denominator — `ambiguous="remove"` (default):** only ACTGSWU are included when
   computing the length; other symbols (e.g. N) are removed from the denominator. This
   matches the Wikipedia A+T+G+C denominator for standard bases.
4. **Denominator — `ambiguous="ignore"`:** all symbols are included in the length;
   only G/C/S count toward GC.
5. **Worked doctest values (verbatim, copied from source):**
   - `gc_fraction("ACTGN", "ignore")` → `0.40` (GC=2, length=5).
   - `gc_fraction("ACTGN", "remove")` → `0.50` (GC=2, length=4 after removing N).
   - `gc_fraction("ACTGN", "weighted")` → `0.50`.
   - `gc_fraction("ACTG")` → `0.50`; RNA `"GGAUCUUCGGAUCU"` → `0.50`.
6. **GC123:** "Calculate G+C content: total, for first, second and third positions",
   returning "a tuple of four floats (percentages between 0 and 100)". Example
   `GC123("ACTGTN")` → `(40.0, 50.0, 50.0, 0.0)`. Confirms the percentage (×100) form is
   the standard external presentation of GC content.

---

## Documented Corner Cases and Failure Modes

### From Biopython (gc_fraction)

1. **Ambiguous bases (N):** under the default `remove` mode, N is excluded from the
   denominator (`"ACTGN" remove → 0.50`); under `ignore` it inflates the denominator
   (`→ 0.40`). The choice changes the value, so it is correctness-affecting.
2. **RNA (U):** U is treated as a non-GC base equivalent to T (`"GGAUCUUCGGAUCU" → 0.50`).

### From Wikipedia (GC-content)

1. **Empty / no standard bases:** the formula is undefined when A+T+G+C = 0 (division by
   zero); a windowed implementation must define a convention (this unit returns 0 for a
   window with no A/T/U/G/C base — see Assumptions).

---

## Test Datasets

### Dataset: Biopython gc_fraction doctests (scaled to percentage ×100)

**Source:** Biopython `Bio.SeqUtils.gc_fraction` doctests (retrieved source file, 2026-06-14)

| Input | Mode | gc_fraction (0–1) | GC% (×100) |
|-------|------|-------------------|-----------|
| `ACTG` | remove | 0.50 | 50.0 |
| `ACTGN` | remove (default) | 0.50 | 50.0 |
| `ACTGN` | ignore | 0.40 | 40.0 |
| `GGAUCUUCGGAUCU` (RNA) | remove | 0.50 | 50.0 |

### Dataset: Hand-derived window values from the Wikipedia formula (×100)

**Source:** derived from `(G+C)/(A+T+G+C)×100` (Wikipedia GC-content, 2026-06-14)

| Window | G+C | A+T+G+C | GC% |
|--------|-----|---------|-----|
| `GGGG` | 4 | 4 | 100.0 |
| `AAAA` | 0 | 4 | 0.0 |
| `ATGC` | 2 | 4 | 50.0 |
| `GGGA` | 3 | 4 | 75.0 |
| `GCAT` | 2 | 4 | 50.0 |
| `GGAN` (N excluded, remove) | 2 | 3 | 66.66666666666666 |

---

## Assumptions

1. **ASSUMPTION: Empty-window convention (window with no standard base)** — Wikipedia
   and Biopython leave GC content undefined when A+T+G+C = 0 (division by zero). The
   repository returns `0` for such a window. This mirrors the sibling unit
   SEQ-GC-ANALYSIS-001 (`GcSkewCalculator`), which returns `GcContent = 0` for "no G/C
   bases" / zero-division, so the convention is consistent within the repository, but it
   is not dictated by the external sources. Used only for the degenerate all-N window case.

2. **ASSUMPTION: Denominator excludes non-standard symbols (N etc.)** — resolved by
   evidence: Biopython default `ambiguous="remove"` and the Wikipedia A+T+G+C denominator
   both exclude N from the denominator. The implementation counts only A/T/U/G/C in the
   denominator, so this matches the `remove` convention; it is therefore source-backed,
   not an open assumption. (Recorded here for traceability of the design choice.)

---

## Recommendations for Test Coverage

1. **MUST Test:** each window's GC% equals `(G+C)/(A+T+G+C)×100` with exact values
   (50, 75, 100, 0) — Evidence: Wikipedia GC-content formula; Biopython gc_fraction ×100.
2. **MUST Test:** N is excluded from the denominator (window `GGAN` → 200/3 = 66.6…%) —
   Evidence: Biopython `ACTGN` remove → 0.50 (N removed from length).
3. **MUST Test:** window count = ⌊(n − w)/step⌋ + 1 and window offsets 0, step, 2·step, …
   — Rationale: sliding-window definition (matches sibling SEQ-ENTROPY-PROFILE-001 INV-05).
4. **MUST Test:** RNA U counts as non-GC (window with U gives same value as T) —
   Evidence: Biopython RNA `GGAUCUUCGGAUCU` → 0.50.
5. **SHOULD Test:** window with no standard base (all-N) → 0 — Rationale: zero-division
   convention (Assumption 1).
6. **SHOULD Test:** windowSize > length, null, empty → empty profile — Rationale:
   no full window exists; guarded input.
7. **COULD Test:** case-insensitivity (lowercase input equals uppercase) — Rationale:
   implementation case-folds before counting.
8. **COULD Test:** every window value is bounded in [0, 100] — Rationale: GC% bound
   from the formula (numerator ≤ denominator).

---

## References

1. Wikipedia. 2026. GC-content. https://en.wikipedia.org/wiki/GC-content (accessed 2026-06-14).
2. Cock P.J.A. et al. 2009. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 ; reference implementation `Bio.SeqUtils.gc_fraction` source: https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
