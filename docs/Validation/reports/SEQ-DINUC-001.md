# Validation Report: SEQ-DINUC-001 — Dinucleotide Analysis

- **Validated:** 2026-06-15   **Area:** Statistics
- **Canonical method(s):** `SequenceStatistics.CalculateDinucleotideFrequencies(string)`,
  `SequenceStatistics.CalculateDinucleotideRatios(string)`
  (codon-frequency `CalculateCodonFrequencies` migrated to SEQ-CODON-FREQ-001 on 2026-06-14 — out of scope here)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

## Stage A — Description

### Sources opened this session
- **Karlin S., "Pervasive properties of the genomic signature", PMC126251** —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC126251/ (WebFetch).
  Verbatim: odds ratio `r_xy = f_xy / (f_x f_y)`, where `f_x` is the *(normalized)* frequency of base `x`
  and `f_xy` the frequency of dinucleotide `xy` "in the leading strand". `r=1` ⇒ frequency matches the
  value expected under statistical independence (no bias); `>1` over-, `<1` under-representation.
  Symmetrized `r*_xy` is "computed from frequencies of the sequence concatenated with its inverted
  complement."
- **Karlin & Burge (1995) thresholds via MBE 19(6):964** — https://academic.oup.com/mbe/article/19/6/964/1095097
  (WebSearch). Verbatim restatement: "the XY dinucleotide was considered to be underrepresented if
  ρ*XY ≤ 0.78 and overrepresented if ρ*XY ≥ 1.23", attributed to Karlin & Burge (1995).
- **Gardiner-Garden & Frommer (1987) CpG O/E** — same odds-ratio shape with the dinucleotide count
  normalized by N rather than N−1 (the N/(N−1) difference). Repository follows the Karlin N−1 convention.

### Formula check
`ρ_XY = f_XY / (f_X · f_Y)` matches Karlin PMC126251 exactly. Single-base frequency `f_X` is normalized
over the count of valid bases (N); dinucleotide frequency `f_XY` is normalized over the N−1 dinucleotide
positions. This mixed normalization is the standard Karlin odds-ratio convention. `ρ = 1` no-bias
baseline confirmed. The 0.78/1.23 classification thresholds are documentation-only (not applied in code),
which is the correct choice — the methods return the raw ρ.

### Edge-case semantics
- null / empty / length < 2 → empty dictionary (input guard) — defined.
- A constituent base absent ⇒ expected = 0. **Mathematical note:** any dinucleotide that *appears* in the
  output necessarily has both constituent bases present, so `f_X·f_Y > 0` always holds for output keys. The
  `expected > 0 ? … : 0` branch is therefore unreachable defensive code, not an observable behaviour. The
  description's "ρ = 0 when expected = 0" is harmless but never triggers; recorded as such (no defect).
- Non-alphabet units ({A,T,G,C,U} for dinucleotides) excluded from counts — defined and standard.
- RNA U treated as a fifth base — not contradicted by any source; RNA signatures are defined identically.

### Independent cross-check (exact rationals, recomputed in Python `fractions`, not from the code)
Sequence `ATGCGCGT` (A=1,T=2,G=3,C=2; N=8; 7 dinucleotide positions):
- ρ_GC = ρ_CG = (2/7)/((3/8)(2/8)) = **64/21 ≈ 3.047619047619**
- ρ_AT = (1/7)/((1/8)(2/8)) = **32/7 ≈ 4.571428571429**
- ρ_TG = ρ_GT = (1/7)/((2/8)(3/8)) = **32/21 ≈ 1.523809523810**
- f_GC = f_CG = **2/7**, f_AT = f_TG = f_GT = **1/7**, Σf = **1.0**
- Homopolymer `AAAA`: f_AA=1, f_A=1 ⇒ ρ_AA = **1.0** (no-bias baseline)
- `ATAT`: ρ_AT = **8/3**, ρ_TA = **4/3** (G,C absent → those keys never appear; all finite)
- `AUGCGC` (RNA): f_AU = **1/5**, f_GC = **2/5**
- `ATGNCG` (ambiguous N): valid dinucleotides AT, TG, CG each f = **1/3**; ρ_AT = **25/3**, ρ_TG = ρ_CG = **25/6**

Every value the tests assert traces to this independent recomputation, not to the implementation output.

### Findings / divergences
None material. N−1 vs N normalization is a documented, authoritative modeling choice (Karlin). The
single-strand (non-symmetrized) ρ is correctly documented as an intentional simplification.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs`:
- `CalculateDinucleotideFrequencies` (lines 602–631): null/empty/<2 guard; scans adjacent pairs; counts
  only `ATGCU` pairs; `freq = count/total` where total = number of valid dinucleotide positions. Realises
  `f_XY = count/(N−1)` exactly when all positions are valid.
- `CalculateDinucleotideRatios` (lines 641–673): null/empty/<2 guard; single-base freqs from
  `CalculateNucleotideComposition` over A,T,G,C,U with total = N; dinucleotide freqs from the method above;
  `ratios[d] = expected > 0 ? observed/expected : 0`. Realises `ρ_XY = f_XY/(f_X f_Y)` exactly.

### Formula realised correctly?
Yes. Traced and confirmed against the independent recomputation above — all values match.

### Cross-verification table recomputed vs code
Ran the fixture: M1/M2/M3/M7/S4/S5/C2/C1 all green; values equal the externally recomputed rationals
(within 1e-10 tolerance for exact double divisions).

### Variant/delegate consistency
`CalculateDinucleotideRatios` consumes `CalculateDinucleotideFrequencies` and
`CalculateNucleotideComposition`; the composition's A,T,G,C,U counts give N (not N−1), which is correct
for single-base frequencies. No `*Fast`/instance variants for these methods.

### Test quality audit (HARD gate)
- **Sourced, not code-echoed:** all expected values are exact rationals derived from Karlin's definition
  and verified by an independent Python computation this session; a wrong implementation would fail them.
- **No green-washing:** no weakened assertions, no widened tolerances, no skipped tests.
- **Coverage fixes applied this session (Stage B defects in the tests, not the code):**
  1. **S4 was mislabeled** — its comment claimed to exercise the "expected = 0 ⇒ ρ = 0" guard, but that
     branch is unreachable (a present dinucleotide always has both bases present). Rewrote the comment to
     state what it actually verifies, and strengthened it: added exact `ρ_TA = 4/3` and an
     `Is.EquivalentTo` key-set assertion (absent bases G, C produce no keys).
  2. **Missing ambiguous-base exclusion coverage** — neither method had a test for non-ACGT exclusion.
     Added `CalculateDinucleotideRatios_AmbiguousBase_ExcludedFromCounts` (`ATGNCG`: keys {AT,TG,CG};
     ρ_AT=25/3, ρ_TG=ρ_CG=25/6) and `CalculateDinucleotideFrequencies_AmbiguousBase_ExcludedFromCounts`
     (`ATGNCG`: each f=1/3, Σ=1, exact key set), values externally recomputed.
- **All branches covered:** guards (S1/S2), valid-pair counting + frequency normalization (M2/M3/C2),
  ambiguous exclusion (new S5/S5f), ratio formula (M1/M7/S4), case-insensitivity (C1), RNA U (C2).
  The only uncovered line is the unreachable `: 0` defensive branch (documented above).
- **Honest green:** FULL unfiltered suite `Failed: 0, Passed: 6523` (1 benchmark skipped by design);
  `dotnet build` 0 errors; changed test file produces no warnings.

### Findings / defects
- **DEF (test-only, FIXED):** S4 comment over-claimed guard coverage; two methods lacked ambiguous-base
  exclusion tests. Fixed in this session (see audit). No implementation defect found.

## Verdict & follow-ups
- **Stage A:** PASS — formula, baseline, thresholds, and all numeric examples confirmed against Karlin
  PMC126251 + Karlin & Burge (MBE restatement) and an independent recomputation.
- **Stage B:** PASS — code faithfully realises the validated formula; tests strengthened to lock sourced
  values and cover the ambiguous-base branch.
- **Test-quality gate:** PASS.
- **End-state:** CLEAN — no implementation defect; the test gaps found were completely fixed and the full
  suite is green.
- **Note:** the `expected = 0 ⇒ ρ = 0` guard is unreachable defensive code (kept; harmless); documented so
  no future session mistakes it for testable behaviour.
