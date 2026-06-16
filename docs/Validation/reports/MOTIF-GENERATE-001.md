# Validation Report: MOTIF-GENERATE-001 — IUPAC-Degenerate Consensus Generation

- **Validated:** 2026-06-16   **Area:** Matching
- **Canonical method(s):** `MotifFinder.GenerateConsensus(IEnumerable<string>)` (+ private `GetIupacCode`)
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session (retrieved, not just cited)

| Source | URL | What it confirms |
|--------|-----|------------------|
| UCSC Genome Browser — IUPAC ambiguity codes | https://genome.ucsc.edu/goldenPath/help/iupac.html | Full single-letter table (verbatim, fetched 2026-06-16) |
| Wikipedia — Nucleic acid notation (Table 1, cites NC-IUB 1984) | https://en.wikipedia.org/wiki/Nucleic_acid_notation | Same table; primary ref = NC-IUB 1984 (Cornish-Bowden, NAR) |
| Bioconductor DECIPHER `ConsensusSequence` | https://rdrr.io/bioc/DECIPHER/man/ConsensusSequence.html | Threshold-consensus mechanism; equal-abundance → degeneracy code; default threshold 0.05 |

### Formula / mapping check (set → IUPAC symbol)

Retrieved verbatim from UCSC and Wikipedia this session; the two agree exactly and Wikipedia
attributes the table to the NC-IUB 1984 primary standard:

| Base set | Symbol | UCSC | Wikipedia |
|----------|:------:|:----:|:---------:|
| {A,G} | R | ✓ | ✓ |
| {C,T} | Y | ✓ | ✓ |
| {C,G} | S | ✓ | ✓ |
| {A,T} | W | ✓ | ✓ |
| {G,T} | K | ✓ | ✓ |
| {A,C} | M | ✓ | ✓ |
| {C,G,T} | B (not-A) | ✓ | ✓ |
| {A,G,T} | D (not-C) | ✓ | ✓ |
| {A,C,T} | H (not-G) | ✓ | ✓ |
| {A,C,G} | V (not-T) | ✓ | ✓ |
| {A,C,G,T} | N (any) | ✓ | ✓ |

Singletons {A},{C},{G},{T} → standard base. The `GetIupacCode` switch matches this table
character-for-character. The core biology (set→symbol) is fully source-backed.

### Edge-case semantics

- Empty collection → `""` (guard contract).
- Null collection → `ArgumentNullException` (guard).
- Case-insensitive: input upper-cased before counting; non-ACGT chars ignored.
- Output length = first sequence length (per-column construction).

### Independent cross-check (hand computation, threshold = n×0.25, include iff count > threshold)

| Column (n) | thr | passing bases | symbol | source for symbol |
|------------|----:|---------------|:------:|-------------------|
| A,G (2) | 0.5 | {A,G} | R | NC-IUB/UCSC |
| C,T (2) | 0.5 | {C,T} | Y | NC-IUB/UCSC |
| C,G,T (3) | 0.75 | {C,G,T} | B | NC-IUB/UCSC |
| A,A,G,G,C (5) | 1.25 | {A(2),G(2)}; C(1) dropped | R | NC-IUB + threshold filter |
| A,C,G,T (4) | 1.0 | none → fallback most-frequent (tie→A) | A | implementation contract |

All trace to the retrieved table / documented threshold.

### Findings / divergences (PASS-WITH-NOTES)

1. **Threshold parameterisation differs from DECIPHER.** DECIPHER *removes* characters that
   represent *less than* `threshold` fraction (default 0.05). This implementation *includes* a
   base iff `count > total × 0.25` (strict `>`, fixed 25%). This is a documented, named design
   constant (`IupacInclusionThreshold = 0.25`), explicitly disclosed in code XML docs, Evidence
   §Assumption-1 and TestSpec §6 — the threshold-consensus *family* is authoritative; the exact
   cut is implementation-specific. Correctness-affecting but disclosed, not invented-untraceable.
2. **N (four-base) is unreachable** under strict `>25%`: four bases cannot each exceed 25%, so
   four-equal columns fall through to the most-frequent fallback rather than emitting N. This is a
   direct consequence of the threshold semantics and is documented (TestSpec §1.4, INV-3 note).
3. **Fallback when no base passes** = single most-frequent base, alphabetical tie-break (`MaxBy`
   over A<C<G<T insertion order). No authoritative spec defines this corner; it is an
   implementation contract, honestly labelled as such.

None of these contradict the authoritative biology; the set→symbol output is dictated by the
NC-IUB table for every case where ≥1 base passes. → **PASS-WITH-NOTES**.

## Stage B — Implementation

### Code path reviewed

`src/.../Seqeron.Genomics.Analysis/MotifFinder.cs:339–486`
(`GenerateConsensus` + `GetIupacCode`, threshold const `IupacInclusionThreshold = 0.25`).

### Formula realised correctly?

Yes. Per-column counts over {A,C,G,T}; `threshold = total*0.25`; bases with `count > threshold`
are ordered alphabetically, joined, and mapped through a switch identical to the NC-IUB/UCSC table
I retrieved this session. Zero-passing → `counts.MaxBy(Value).Key` (first max = alphabetically
earliest given A,C,G,T insertion order).

### Cross-verification table recomputed vs code (test run)

| Test | Input | Expected (sourced) | Code | Match |
|------|-------|--------------------|------|:-----:|
| M1–M6 | 2-base columns | R,Y,S,W,K,M | same | ✓ |
| M7–M10 | 3-base columns | B,D,H,V | same | ✓ |
| M11 | ATGC×3 | ATGC | same | ✓ |
| M12 | ATGC,GTGC | RTGC | same | ✓ |
| M13 | AAAA,AAGT,AACT,AATT | A,A,A,T | same | ✓ |
| M14 | A,A,G,G,C | R | same | ✓ |
| M15 | A,C,G,T | A (fallback) | same | ✓ |
| S1 | atgc,gtgc | RTGC | same | ✓ |
| S2 | [] | "" | same | ✓ |
| S3 | ACGTACG×2 | len 7 | same | ✓ |
| C1 | null | ArgumentNullException | same | ✓ |

### Variant / delegate consistency

`CreateConsensusFromAlignment` is a *separate* most-common-base (Rosalind CONS) method, correctly
distinguished in XML docs; not part of this unit's degenerate-consensus contract. No delegate drift.

### Test quality audit (HARD gate)

- **Sourced expectations:** M1–M12 lock exact NC-IUB symbols retrieved this session; a wrong
  mapping (e.g. AG→Y) fails M1. No code-echo assertions.
- **No green-washing:** all assertions are exact `Is.EqualTo`; no Greater/AtLeast/Contains/ranges;
  no widened tolerance; no skipped/ignored tests. Old weak tests (`.Or.EqualTo`, `Does.Match`,
  `BeOneOf`) in `MotifFinderTests.cs` / `MutationKillerTests.cs` were already removed and replaced
  by exact-value cases (verified: only NOTE comments remain).
- **Coverage:** every 2-base and 3-base symbol branch, singletons, multi-column, strict-25%
  boundary (M13), minority-drop (M14), fallback (M15), case (S1), empty (S2), length (S3), null
  (C1). The `_ => 'N'` branch is unreachable under the threshold (documented), so the fallback is
  exercised instead — no live path left untested.
- **Honest green:** FULL unfiltered suite = **Failed: 0, Passed: 6606**. Build: 0 errors
  (4 pre-existing warnings in unrelated files).

→ **Test-quality gate: PASS.**

### Findings / defects

None. No code or test changes required.

## Verdict & follow-ups

- **Stage A:** PASS-WITH-NOTES (threshold is an implementation-specific design constant; N
  unreachable; fallback is an implementation contract — all documented, biology is source-exact).
- **Stage B:** PASS. Code faithfully realises the validated mapping; tests are exact and sourced.
- **End-state:** ✅ CLEAN — no defect found; algorithm fully functional.
- **Test-quality gate:** PASS.
- No findings logged in FINDINGS_REGISTER (no defect).
