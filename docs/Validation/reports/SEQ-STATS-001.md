# Validation Report: SEQ-STATS-001 — Sequence Composition Statistics (+ opt-in Biopython conventions)

- **Validated:** 2026-06-24   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateNucleotideComposition(string)`; delegates `SummarizeNucleotideSequence(string?)`, `CalculateAminoAcidComposition(string)`. Opt-in convention surfaces (commit `6e900e92`): `SequenceStatistics.CalculateGcContentProfile(..., bool fraction=false)`, `SequenceExtensions.CalculateGcFraction(this string/ReadOnlySpan<char>, GcAmbiguityMode)`.
- **Stage A verdict:** PASS-WITH-NOTES
- **Stage B verdict:** PASS

## Scope of this revalidation

SEQ-STATS-001 was reset by the limitations campaign (commit `6e900e92`, "opt-in
Biopython/VCF compatibility modes"). The core composition statistics (length, base counts,
GC content, GC/AT skew) were already validated CLEAN on 2026-06-15; this session re-confirms
them and validates the **new opt-in convention behaviour** added by that commit:

1. `CalculateGcContentProfile(..., fraction:true)` — emit GC as a [0,1] fraction instead of
   the default percentage [0,100].
2. `CalculateGcFraction(GcAmbiguityMode {Remove,Ignore,Weighted})` — Biopython `gc_fraction`
   IUPAC-ambiguity parity overload.

Defaults are unchanged: GC profile defaults to percentage, the parameterless
`CalculateGcFraction()` keeps the existing A/T/G/C/U-only behaviour.

## Stage A — Description

### Sources opened this session (with extracted numbers)

1. **Biopython `Bio.SeqUtils.gc_fraction` + `_gc_values`** — raw master source fetched
   2026-06-24 (https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py):
   - Numerator (all modes start): `gc = sum(seq.count(x) for x in "CGScgs")`.
   - **remove** (default): `length = gc + sum(seq.count(x) for x in "ATWUatwu")`.
   - **ignore**: `length = len(seq)` (no extra numerator).
   - **weighted**: `gc += Σ count(x)·_gc_values[x]` for `x in "BDHKMNRVXY"`; `length = len(seq)`.
   - `_gc_values = {G:1, C:1, A:0, T:0, U:0, S:1, W:0, M:0.5, R:0.5, Y:0.5, K:0.5, V:0.667,
     B:0.667, H:0.333, D:0.333, X:0.5, N:0.5}`.
   - "Note that this will return zero for an empty sequence." Returns a float in [0,1].
2. **Wikipedia "GC content"** — GC content = (G+C)/(A+T+G+C); standard fraction in [0,1].
3. (Carried, unchanged from 2026-06-15) Wikipedia "GC skew" + Lobry 1996 for skew formulas.

### Formula check
- GC content (default) = (G+C)/(A+T+G+C+U), [0,1] — matches Wikipedia and Biopython remove
  mode over the standard alphabet. ✓
- Opt-in fraction profile: identical windows to the percentage profile, divided by 100 (scale
  1.0 vs 100). A correct, source-cited representation toggle. ✓
- `CalculateGcFraction(mode)` reproduces Biopython `gc_fraction` exactly per the three modes
  and `_gc_values` (see cross-check below). ✓
- GC skew = (G−C)/(G+C), AT skew = (A−T)/(A+T), zero-denominator → 0 — unchanged, sourced. ✓

### Edge-case semantics
- Empty/null → 0 in every GC-fraction mode (sourced to Biopython "return zero for an empty
  sequence"). ✓
- Zero-length denominator (e.g. `VH` in Remove: no A/T/G/C/S/W/U) → 0 (guarded). ✓
- Default (percentage / A,T,G,C,U-only) is explicitly preserved; opt-in is by-design. ✓

### Independent cross-check (hand computation vs Biopython gc_fraction)

| Input | default frac | Remove | Ignore | Weighted | Biopython reasoning |
|-------|-------------|--------|--------|----------|---------------------|
| `GCAT`   | 0.5 | 0.5    | 0.5    | 0.5    | gc=2, no ambiguity, len=4 |
| `GCATN`  | 0.5 | **0.5** (2/4) | **0.4** (2/5) | **0.5** (2.5/5) | N: dropped / in-len / +0.5 |
| `GCVBHD` | 1.0 | **1.0** (2/2) | **0.333** (2/6) | **0.667** (4/6) | V,B=2/3; H,D=1/3 |
| `GCSW`   | 1.0 | 0.75 (3/4) | 0.75 | 0.75 | S→num+len, W→len only |
| `VH`     | 0.0 | 0.0 | 0.0 | **0.5** (1/2) | weighted V=2/3,H=1/3 |
| `""`     | 0.0 | 0.0 | 0.0 | 0.0 | empty → 0 |

Default (percentage) for the profile: windows `GC,AT,GC` → `100,0,100`; with `fraction:true`
→ `1.0,0.0,1.0`. `GCATGCAT` overall → 50% (default) vs 0.5 (fraction). All confirmed.

### Findings / divergences (Stage A → PASS-WITH-NOTES)
- **N1 (carried, documented):** the *default* composition counts degenerate IUPAC codes as
  Other and excludes them from GC/AT totals (differs from Biopython on degenerate symbols).
  This is now explicitly addressed by the opt-in `CalculateGcFraction(GcAmbiguityMode)`
  overload, which reproduces Biopython's degenerate handling exactly. Documented intentional
  convention, not a defect.
- **N2 (carried, documented):** `AtContent` uses (A+T+U)/total for RNA; AtSkew uses the
  DNA-specific (A−T)/(A+T). Defensible documented convention. No correctness issue.

## Stage B — Implementation

### Code path reviewed
- `src/Seqeron/Algorithms/Seqeron.Genomics.Core/SequenceExtensions.cs`:
  - `CalculateGcFraction(ReadOnlySpan<char>)` lines 81–111 (default, A/T/G/C/U only).
  - `enum GcAmbiguityMode` lines 136–155; `_gc_values` constants 160–164.
  - `CalculateGcFraction(ReadOnlySpan<char>, GcAmbiguityMode)` lines 181–214; string overload
    220–223; `WeightedGcValue` 225–233.
- `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs`:
  `CalculateGcContentProfile(..., bool fraction=false)` lines 905–944 — `scale = fraction ? 1.0
  : PercentScale`, applied at the single `yield return` (line 940).

### Formula realised correctly? (evidence)
- Remove: numerator counts G/C/S (case-insensitive via `char.ToUpperInvariant`); A/T/U/W add to
  `strongWeakCount` (the Remove denominator); everything else excluded from both. `length =
  strongWeakCount`. Matches Biopython remove. ✓
- Ignore/Weighted: `length = totalLength` (full sequence). Weighted adds `WeightedGcValue` for
  the default (non-G/C/S/A/T/U/W) branch only. Matches Biopython ignore/weighted. ✓
- `WeightedGcValue`: V/B=2/3, H/D=1/3, M/R/Y/K/X/N=0.5, else 0 — verbatim `_gc_values`. Note S
  is handled in the strong branch (=1) and W in the weak branch (=0), matching Biopython's
  split (S in "CGS", W in "ATWU", and W also has _gc_values[W]=0 so the omission is exact). ✓
- Empty/null → 0 (both overloads guard); zero denominator → 0. ✓
- Profile fraction toggle: pure scale factor; defaults `false` → unchanged percentage. ✓

### Cross-verification table recomputed vs code
Executed the actual compiled `Seqeron.Genomics.Core` against the table above (small driver,
since deleted). Output matched every hand value exactly: `GCATN` Remove/Ignore/Weighted =
0.5 / 0.4 / 0.5; `GCVBHD` = 1.0 / 0.333 / 0.667; `GCSW` = 0.75; `VH` Weighted = 0.5; empty = 0.
Default parameterless `GCAT` = 0.5 (unchanged). ✓

### Variant/delegate consistency
- Default `CalculateGcFraction()` and `CalculateGcFractionFast(string)` agree (the latter
  forwards to the span overload). The mode overload is a separate opt-in surface and does not
  alter the default path. ✓
- `CalculateGcContentProfile(fraction:true)` yields the same windows as `fraction:false` scaled
  by 1/100 (test asserts `[100,0,100]` vs `[1.0,0.0,1.0]`). ✓

### Test quality audit (HARD gate)
- `ConventionCompatibility_OptIn_Tests.cs` (the relevant SEQ-STATS-001 cases): Remove
  unambiguous 0.5, Remove S/W split 0.75, Remove N-excluded 0.5, Ignore 0.4, Weighted N=0.5,
  Weighted V/H 0.5, empty/null 0, default-overload-unchanged 0.5, profile fraction `[1,0,1]`
  vs `[100,0,100]`. All `Is.EqualTo(...).Within(1e-10)` against externally-sourced Biopython
  values — would fail a wrong implementation. ✓
- Canonical `SequenceStatistics_CalculateNucleotideComposition_Tests.cs` (core stats) retains
  exact-value assertions (counts partition, GC content 0.5/1.0, GC skew ±0.5, AT skew 0.5/1.0,
  empty/null 0). ✓
- Filtered run: **Failed: 0, Passed: 31** (canonical composition tests + convention tests).

### Findings / defects
None. The opt-in surfaces reproduce Biopython `gc_fraction` exactly; defaults are unchanged.

## Verdict & follow-ups
- **Stage A:** PASS-WITH-NOTES (carried documented conventions N1/N2; the new opt-in API
  directly addresses N1 with exact Biopython parity).
- **Stage B:** PASS — implementation faithfully realises the validated formulas and the
  Biopython ambiguity modes; tests assert exact externally-sourced values and cover all modes
  + the default-unchanged guards.
- **End-state:** CLEAN — fully functional; no defect found; no code changed this session.
- Build: 0 errors / 0 warnings. Filtered SEQ-STATS-001 + ConventionCompatibility tests: 31/0.
