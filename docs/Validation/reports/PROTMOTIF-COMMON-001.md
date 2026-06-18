# Validation Report: PROTMOTIF-COMMON-001 — Common Motif Finding

- **Validated:** 2026-06-16   **Area:** ProteinMotif
- **Canonical method(s):** `ProteinMotifFinder.FindCommonMotifs(string)` (delegates per-pattern to `FindMotifByPattern`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Scope note (prompt framing vs. actual algorithm)

The session prompt framed this unit as "finding common/shared protein motifs **across multiple
sequences**" and pointed at the LCSM / longest-common-substring family (Rosalind LCSM, Biopython).
That is **not** what `FindCommonMotifs` computes. The repo's `FindCommonMotifs` takes **one** protein
sequence and scans it against a curated library of well-known PROSITE-style functional motifs
(N-glycosylation, PKC/CK2 phosphorylation, P-loop, RGD, …), returning every occurrence of every
library pattern. "Common" here means "common/well-known functional motifs", not "common substring
shared by a set of sequences". I validated the algorithm the code actually implements (and that the
TestSpec/Evidence/algorithm doc consistently describe); the LCSM interpretation does not apply.

## Stage A — Description

### Sources opened & what they confirm (retrieved this session, 2026-06-16)

| Source | URL | Extracted (verbatim) |
|--------|-----|----------------------|
| PROSITE PS00001 | https://prosite.expasy.org/PS00001 | name `ASN_GLYCOSYLATION`, pattern `N-{P}-[ST]-{P}` |
| PROSITE PS00005 | https://prosite.expasy.org/PS00005 | name `PKC_PHOSPHO_SITE`, pattern `[ST]-x-[RK]` |
| PROSITE PS00006 | https://prosite.expasy.org/PS00006 | name `CK2_PHOSPHO_SITE`, pattern `[ST]-x(2)-[DE]` |
| PROSITE PS00016 | https://prosite.expasy.org/PS00016 | name `RGD`, pattern `R-G-D.` |
| PROSITE PS00017 | https://prosite.expasy.org/PS00017 | name `ATP_GTP_A`, pattern `[AG]-x(4)-G-K-[ST]` |
| ScanProsite doc | https://prosite.expasy.org/scanprosite/scanprosite_doc.html | syntax `x`, `[..]`, `{..}`, `x(n)`, `x(n,m)`, `-`; default matching "greedy, allows overlaps but not included matches — two overlapping matches are rejected if one is entirely contained within the other" |

All five bundled PROSITE patterns in `CommonMotifs` match their official ExPASy entries **exactly**
(accession, entry name, and PA-line). The PROSITE→regex translations in the code are faithful:
`{P}`→`[^P]`, `[ST]`→`[ST]`, `x`→`.`, `x(2)`→`.{2}`, `x(4)`→`.{4}`.

### Formula / semantics check

- Element semantics (`[..]` allowed set, `{..}` excluded set, `x` wildcard, `x(n)`/`x(n,m)`
  repetition, `-` separator) match the ScanProsite documentation verbatim.
- Overlap reporting ("overlaps, no includes") matches the ScanProsite default; the lookahead
  `(?=(...))` engine reproduces it (a fully-contained match cannot arise from a fixed/min-width
  pattern advanced one position at a time, so "no includes" is automatically satisfied for these
  bundled patterns).
- Coordinate convention: PROSITE reports **1-based inclusive**; the repo `MotifMatch.Start/End` are
  **0-based inclusive**. This is a documented API-shape assumption, not a correctness defect — matched
  substring content and relative positions are identical (shift of −1 vs a ScanProsite report). The
  ScanProsite page I fetched did not print the 1-based statement explicitly, but 1-based inclusive is
  the universal PROSITE/UniProt convention; the divergence is correctly disclosed in spec & doc.

### Edge-case semantics

Null/empty → empty (guard); Pro at a `{P}` position rejects the N-glyco window (PS00001); lowercase
input upper-cased before matching; no-motif sequence → empty. All defined and sourced.

### Independent cross-check (hand-computed this session, not from the code)

Using a standalone Python `re` implementation of the same PROSITE→regex + `(?=(...))` overlap
translation, I recomputed every dataset:

| Sequence | Pattern | Independent result | TestSpec expected | Match |
|----------|---------|--------------------|-------------------|-------|
| `AAAANFTAAAA` | PS00001 | 4..7 `NFTA` | 4..7 `NFTA` | ✓ |
| `AAAANPSAAAAANPTAAA` | PS00001 | (none) | (none) | ✓ |
| `AAAAASARKAAA` | PS00005 | 5..7 `SAR` | 5..7 `SAR` | ✓ |
| `AAAASAAEASDEDAAA` | PS00006 | 4..7 `SAAE`; 9..12 `SDED` | same | ✓ |
| `AAAAAGXXXXGKSAAAA` | PS00017 | 5..12 `GXXXXGKS` | same | ✓ |
| `AARGDKK` | PS00016 | 2..4 `RGD` | 2..4 `RGD` | ✓ |
| `RGDRGD` | PS00016 | 0..2; 3..5 `RGD` | same | ✓ |
| `RGDNFTA` (whole library) | all | exactly RGD 0..2 + ASN_GLYCOSYLATION 3..6 | (was Contains) | ✓ |
| `STRK` (overlap) | PS00005 | 0..2 `STR`; 1..3 `TRK` | (new) | ✓ |

Findings: none. Description is biologically and mathematically correct.

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/ProteinMotifFinder.cs`
- `CommonMotifs` library — lines 79–141 (patterns + regex)
- `FindCommonMotifs` — lines 150–164 (null/empty guard, upper-case, iterate library, delegate)
- `FindMotifByPattern` — lines 171–208 (lookahead overlap engine, 0-based Start/End)

### Formula realised correctly?

Yes. Each `CommonMotifs` entry stores both the PROSITE PA-line and its equivalent .NET regex; the five
audited regexes are exact translations of the official patterns. `FindMotifByPattern` wraps the regex
in `(?=(...))` with `RegexOptions.IgnoreCase`, yielding 0-based inclusive `Start`/`End` and the
captured substring — reproducing PROSITE "overlaps, no includes" reporting.

### Cross-verification vs code

Ran the canonical test file: every M/S/C dataset above produces exactly the independently hand-computed
values. The strengthened M7 (exact 2-hit assertion) and new M8b (overlapping windows) both pass against
the real `FindCommonMotifs`, confirming the C# `Regex` engine agrees residue-for-residue with the
independent Python computation.

### Variant/delegate consistency

`FindCommonMotifs` is the sole aggregator; per-pattern matching delegates to `FindMotifByPattern`
(canonical under PROTMOTIF-PATTERN-001). Identity propagation verified: each hit carries the library
entry's `Name`→`MotifName` and `Accession`→`Pattern` (S2, and now M7 also asserts both fields).

### Numerical robustness

No precision/overflow concern: matching is regex over a bounded-width pattern set; coordinates are ints
within sequence length. Null/empty guarded.

### Test quality audit (HARD gate)

| Check | Result |
|-------|--------|
| Sourced expectations, not code echoes | PASS — every expected value traces to a PROSITE PA-line I fetched this session and/or an independent Python recomputation, not to the implementation's output. |
| No green-washing (no weak asserts where exact value known) | **One weakness found & fixed:** M7 used `Does.Contain` although `RGDNFTA` yields *exactly* two library hits — strengthened to an exact ordered 2-hit assertion (count + name + accession + start/end + sequence). |
| Cover all logic / Stage-A branches | **Gap found & filled:** INV-03 "overlapping occurrences reported" was only exercised by *adjacent* (non-overlapping) RGDs (M8). Added M8b: genuinely overlapping `[ST]-x-[RK]` windows `STR`(0..2) and `TRK`(1..3) over `STRK`, both reported. All other branches (allowed/excluded set, wildcard, fixed/variable gap, multi-pattern, multi-occurrence, identity, null/empty, case, negative control, determinism, substring invariant) already covered. |
| Honest green (full unfiltered suite, Failed: 0) | PASS — `dotnet build` 0 errors; full suite **6584 passed, 0 failed**, 1 pre-existing skipped benchmark unrelated to this unit. |

Defects in production code: **none**. Test defects: two test-quality weaknesses (weak M7 assertion;
missing genuinely-overlapping case), both fixed in-session.

## Verdict & follow-ups

- **Stage A:** PASS — all five bundled PROSITE patterns and the ScanProsite syntax/overlap semantics
  independently confirmed; every dataset hand-computed and matches.
- **Stage B:** PASS — implementation faithfully realises the validated description; cross-check exact.
- **Test-quality gate:** PASS after fixes (strengthened M7 to exact; added overlapping-occurrence M8b).
- **End-state:** CLEAN — no production-code defect; the two test weaknesses were completely fixed and
  the full unfiltered suite is green.
- **Note:** the prompt's "common substring across multiple sequences / LCSM" framing does not match
  this unit; `FindCommonMotifs` is a one-sequence PROSITE-library scan, validated as such.
