# Validation Report: SEQ-GC-PROFILE-001 — GC Content Profile (sliding-window GC content)

- **Validated:** 2026-06-16   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateGcContentProfile(string sequence, int windowSize = 100, int stepSize = 1)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES

## Stage A — Description

### Sources opened (retrieved this session, 2026-06-16)

1. **Wikipedia — GC-content** (https://en.wikipedia.org/wiki/GC-content), fetched via WebFetch.
   - Verbatim formula: numerator `G + C`, denominator `A + T + G + C`, `× 100 %`.
   - Confirms positional/local variation: "These variations in GC-ratio within the genomes
     of more complex organisms result in a mosaic-like formation with islet regions called
     isochores." and "genes are often characterised by having a higher GC-content in contrast
     to the background GC-content for the entire genome." This justifies evaluating GC on a
     sliding window to produce a positional profile.
2. **Biopython `Bio.SeqUtils.gc_fraction`** source
   (https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py),
   fetched via WebFetch.
   - Docstring: "Calculate G+C percentage in seq (float between 0 and 1)."
   - Doctests (verbatim values this session): `gc_fraction("ACTG") → 0.50`;
     `gc_fraction("GGAUCUUCGGAUCU") → 0.50`; `gc_fraction("ACTGN","remove") → 0.50`;
     `gc_fraction("ACTGN","ignore") → 0.40`; `gc_fraction("ACTGN","weighted") → 0.50`.
   - Numerator always counts **G, C, S**. Denominator by mode:
     - `remove` (default): G+C+S+A+T+W+U — only unambiguous nucleotides; **N excluded**.
     - `ignore`: full sequence length (N inflates denominator → 0.40).
     - `weighted`: full length, ambiguous codes contribute fractional GC.
3. **Sliding-window window count** (WebSearch, general DS/algorithms references): the count for
   window `w`, step `s` over length `n` is the standard `⌊(n − w)/s⌋ + 1` for `w ≤ n`.

### Formula check
The unit's formula `GC% = (G + C) / (A + T + G + C) × 100` matches the Wikipedia formula
exactly (numerator, denominator, ×100). Biopython returns the same quantity as a fraction in
[0,1]; ×100 yields the percentage form, consistent with Biopython `GC123` ("percentages
between 0 and 100"). The percentage (0–100) presentation is a documented, accepted deviation
from `gc_fraction`'s [0,1] and is consistent with sibling SEQ-GC-ANALYSIS-001.

### Edge-case semantics
- **N excluded from denominator** — source-backed by Biopython default `remove`
  (`ACTGN → 0.50`, N removed from length) and the Wikipedia A+T+G+C standard-base denominator.
- **U is non-GC, equivalent to T** — Biopython RNA `GGAUCUUCGGAUCU → 0.50`.
- **Window with no standard base (all-N) → 0** — NOT dictated by the external sources (the
  formula is undefined, division by zero). Documented as repository Assumption A1, matching
  SEQ-GC-ANALYSIS-001 (`GcContent = 0` for no-G/C / zero-division). Acceptable convention.
- **windowSize > length, null, empty → empty profile** — standard guarded input; consistent
  with the sliding-window definition (no full window exists).

### Independent cross-check (numbers, hand-derived from the Wikipedia formula + Biopython)
| Window | G+C | denom (A+T+U+G+C) | GC% | Source |
|--------|-----|-------------------|-----|--------|
| GGGG | 4 | 4 | 100.0 | Wikipedia |
| ATGC | 2 | 4 | 50.0 | Wikipedia / Biopython ACTG→0.50 |
| GGGA | 3 | 4 | 75.0 | Wikipedia |
| AAAT | 0 | 4 | 0.0 | Wikipedia |
| TGCC | 3 | 4 | 75.0 | Wikipedia |
| GGAN | 2 | 3 (N excluded) | 66.66666666666666 | Biopython remove ACTGN→0.50 |
| GGAU | 2 | 4 (U non-GC) | 50.0 | Biopython GGAUCUUCGGAUCU→0.50 |
Window counts (n=10, w=4): step 1 → 7, step 2 → ⌊6/2⌋+1 = 4, step 3 → ⌊6/3⌋+1 = 3.

### Findings / divergences
- **Scope note (not a defect against the test set):** Biopython's `remove` denominator also
  includes the ambiguity codes **S** (G or C) and **W** (A or T), and counts **S** toward the
  GC numerator. The repository implementation counts only A/T/U/G/C — S and W are treated as
  non-standard and excluded from the denominator. The unit's documentation explicitly scopes
  itself to "N-style exclusion only; other IUPAC degenerate codes are excluded rather than
  fractionally weighted" (GC_Content_Profile.md §5.3, §6.2). No test case in the spec uses S
  or W, so no expected value is contradicted; the divergence is documented and within scope.
  This is the basis for the Stage-A PASS (with the divergence noted) rather than a FAIL.

**Stage A verdict: PASS** — formula and edge conventions are correct and externally sourced;
the S/W scope limitation is documented and does not affect any sourced expected value.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs:905-936`.

```
guard: null/empty OR windowSize > length → yield break
upper = ToUpperInvariant()
for i = 0; i <= length - windowSize; i += stepSize:
    count gc (G/C) and total (A/T/U/G/C) over [i, i+windowSize)
    yield total > 0 ? gc/total*100 : 0
```

### Formula realised correctly?
Yes. `gc` increments only on G/C; `total` increments on G/C/A/T/U (N and any other symbol
excluded → matches Biopython `remove` for N). `gc/total*100` is the exact Wikipedia
percentage. `total == 0 → 0` realises Assumption A1. Case-folded via `ToUpperInvariant`.
The loop `i += stepSize` over `i <= length - windowSize` yields exactly `⌊(n−w)/step⌋ + 1`
windows at offsets `0, step, 2·step, …`.

### Cross-verification table recomputed vs code (full suite run, Failed: 0)
Every value in the Stage-A table is reproduced by the code via the passing tests:
M1 GGGGGGGGGG→[100.0]; M2 ATGC×3 (w4 s4)→[50,50,50]; M3 GGGAAATGCC (w4 s3)→[75,0,75];
M4 GGAN→[200/3]; M5 GGAU→[50]; M6 counts 7/4/3; S1 empty; S2 GGCC→[100]; S3 empty; S4 NNNN→[0];
C1 lowercase==uppercase; C2 bounded [0,100]. All match the external references.

### Variant/delegate consistency
Single public method, single overload. The only other caller is the MCP passthrough
`AnalysisTools.cs:219`. No `*Fast`/instance variant to reconcile.

### Test quality audit (against external sources, not the code)
- **Sourced, not code-echoes:** expected values (100.0, 50.0, 75.0, 0.0, 200/3, 50.0) are
  derived by hand from the Wikipedia formula / Biopython doctests, computed independently of
  the implementation. A wrong implementation returning fractions (1.0/0.75) or counting N
  (50.0 for GGAN) or counting U as GC (75.0 for GGAU) would FAIL these tests — they
  discriminate against the documented pitfalls.
- **No green-washing:** exact values pinned with `.Within(1e-10)`; no Greater/AtLeast/range
  where an exact value is known (M2/C2 use `All(...)` but each is paired with an exact length
  or exact bound check). No skipped/ignored tests. M4 uses exact `200.0/3.0`.
- **Coverage:** the single public method's happy path, all Stage-A branches (GC count, N
  exclusion, U-as-non-GC, zero-division→0), and edge/error cases (null, empty, window>length,
  window==length, all-N, case-folding, bound) are all exercised. Window-count invariant
  covered for steps 1/2/3; offsets implicitly verified by M3's ordered [75,0,75] under step 3.
- **Honest green:** FULL unfiltered suite `dotnet test` = **Failed: 0, Passed: 6617**;
  `dotnet build` = 0 errors. (4 pre-existing NUnit2007 warnings are in the unrelated file
  ApproximateMatcher_EditDistance_Tests.cs — not touched by this unit.)

**Test-quality gate: PASS.**

### Findings / defects
No defect. One documented scope limitation (S/W ambiguity codes not handled like Biopython
`remove`) — already disclosed in the algorithm doc; out of scope for this unit's spec and not
contradicting any sourced value. This is the reason Stage B is PASS-WITH-NOTES rather than PASS.

## Verdict & follow-ups
- **Stage A:** PASS
- **Stage B:** PASS-WITH-NOTES (documented S/W scope limitation; no incorrect values)
- **End-state:** CLEAN — no defect required fixing; build + full suite green (Failed: 0).
- **Follow-up (optional, non-blocking):** if full Biopython `remove` parity is later desired,
  count S toward GC and S/W toward the denominator, and add S/W test cases sourced from
  `gc_fraction("ACTGSSSS","remove") → 0.75`.
