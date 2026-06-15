# Validation Report: ASSEMBLY-STATS-001 — Assembly Statistics (N50/L50/Nx/Lx/auN, gaps)

- **Validated:** 2026-06-15   **Area:** Assembly
- **Canonical method(s):** `GenomeAssemblyAnalyzer.CalculateStatistics`, `CalculateNx` (3-arg core + 2-arg delegate), `CalculateN50`, `CalculateAuN`, `FindGaps` (plus the trivial `CalculateNxCurve` wrapper).
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN
- **Test-quality gate:** PASS

## Stage A — Description

### Sources opened this session (live)
- **QUAST `quast_libs/N50.py`** (https://raw.githubusercontent.com/ablab/quast/master/quast_libs/N50.py) — fetched the verbatim `NG50_and_LG50` loop:
  ```python
  s = reference_length
  limit = reference_length * (100.0 - percentage) / 100.0
  lg50 = 0
  for l in numlist:
      s -= l
      lg50 += 1
      if s <= limit:
          ng50 = l
          return ng50, lg50
  ```
  `au_metric` = sum of squared contig lengths / reference length. Empty list → `N50_and_L50` returns `(None, None)`.
- **Wikipedia "N50, L50, and related statistics"** (https://en.wikipedia.org/wiki/N50,_L50,_and_related_statistics) — N50 = "sequence length of the shortest contig at 50% of the total assembly length"; L50 = "count of smallest number of contigs whose length sum makes up half of genome size"; N90 = length for which contigs ≥ that length contain at least 90% of total. Worked examples: **Assembly A {80,70,50,40,30,20}, total 290 → N50=70, L50=2**; **Assembly B {…,10,5}, total 305 → N50=50, L50=3**.
- **Heng Li (2020) auN blog** (https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity) — Nx = "Contigs no shorter than Nx cover x% of the assembly"; **auN = ∑ᵢ Lᵢ² / ∑ⱼ Lⱼ** (no sorting needed).
- Miller, Koren & Sutton (2010) §1.2 (cited in Evidence) — "the smallest contig … whose combined length represents at least 50% of the assembly" (inclusive ≥).

### Formula check
- **Nx/Lx inclusive boundary.** QUAST stops at `s <= limit` with `s = total − cumulative`, `limit = total·(100−x)/100`. Algebra: `total − cumulative ≤ total·(100−x)/100` ⇔ `100·cumulative ≥ total·x`. This is exactly the integer test in the code (`cumulative*100 >= totalLength*threshold`). Confirmed inclusive "at least x%", largest-first.
- **auN = Σl²/Σl** confirmed by both lh3 and QUAST `au_metric`.

### Edge-case semantics
- Empty input: authoritative tools return `None`; the repo deliberately returns zeros (Nx=Lx=0, auN=0, all-zero `AssemblyStatistics`). Documented API choice, non-correctness-affecting. Accepted.
- Inclusive boundary (cumulative exactly = x%) selects that contig — sourced.
- Monotonicity N90 ≤ N50, L90 ≥ L50 — sourced (Wikipedia).
- FindGaps coordinate contract (0-based inclusive [Start,End], Length=End−Start+1, maximal N/n runs) is an implementation contract, not externally standardized — verified by hand.

### Independent cross-check (hand-computed this session)
| Input | Metric | Hand value | Source |
|---|---|---|---|
| A {80,70,50,40,30,20}, tot 290 | N50 / L50 | 80+70=150 ≥ 145 → **70 / 2** | Wikipedia (live) |
| A, x=90 (261) | N90 / L90 | 240<261, 270≥261 → **30 / 5** | Wikipedia def |
| A | auN | 16700/290 ≈ 57.586 | lh3 / QUAST |
| B {…,10,5}, tot 305 | N50 / L50 | 200 ≥ 152.5 → **50 / 3** | Wikipedia (live) |
| {100,80,60,40,20}, tot 300 | auN | 22000/300 = 73.333… | lh3 / QUAST |
| {50,50}, tot 100 | N50 / L50 | cum 50 = 50% (inclusive) → **50 / 1** | Miller / QUAST |
| "ACGTNNACGTNNNNNNACGT", minGap=5 | gaps | [4,5] len2 filtered; [10,15] len6 kept | hand |
| "ACGTNNNACGTNNNNNNACGT" | gaps | [4,6] len3; [11,16] len6 | hand |

All match the spec/Evidence and the live external sources.

### Findings / divergences
None. Description is mathematically and biologically correct, every formula traces to a source retrieved this session.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs`
- `CalculateNx` core (L237-264): integer-exact inclusive test `cumulative*100 >= totalLength*threshold` — provably equal to QUAST's `s <= limit`. `long` accumulation avoids overflow.
- `CalculateNx(IEnumerable,int)` (L276-281) and `CalculateN50` (L294-297) delegate correctly (sort descending + total → core; N50 = `.Nx`).
- `CalculateAuN` (L321-332): `Σ((double)l*l)/total` = Σl²/Σl; empty/zero-total → 0.
- `CalculateStatistics` (L121-165): aggregates N50/L50 (x=50), N90/L90 (x=90), largest=`lengths[0]`, smallest=`lengths[^1]`, totals, GC over non-N bases, gap stats; field order matches the record (L29-44) and the tests.
- `FindGaps` (L341-369): sentinel-terminated scan over maximal N/n runs, yields 0-based inclusive [Start,End], `Length = i − gapStart`, with `minGapLength` filter. Matches the contract.

### Formula realised correctly?
Yes. The inclusive boundary, the auN formula, and the largest-first cumulative selection all match the validated description exactly. Post-loop fallback in `CalculateNx` is dead code for thresholds ≤ 100 (the full sum always satisfies the test) and harmless.

### Cross-verification vs code
The full unfiltered suite (which runs all 23 canonical tests) recomputes every row above against the actual code — all green (see below).

### Variant/delegate consistency
2-arg `CalculateNx`, `CalculateN50`, and `CalculateStatistics` all produce the same N50/L50 as the core overload on shuffled and contig-realised inputs (M7, M8, M9). Confirmed.

### Test quality audit
- File `tests/Seqeron/Seqeron.Genomics.Tests/GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs`, 23 `[Test]` methods.
- Sourced expectations, not code echoes: exact `Is.EqualTo` for N50/L50/N90/L90/coords; `Within(1e-10)` (tight) for auN/GC/gap%. Values trace to Wikipedia/lh3/QUAST, not to current output.
- No green-washing: no `GreaterThan`/`AtLeast`/ranges where an exact value is known; M4 monotonicity is a genuine mathematical property (`<=`/`>=`) with the exact bounds named in the message — appropriate for a property test.
- Coverage: all 6 canonical methods + all Stage-A branches (inclusive boundary M5, single contig S1, empty C1/C2/C3, all-N C4, leading/trailing/interior/multiple gaps, minGap filter). Both worked examples (A and B) present.
- Honest green: full unfiltered suite **Failed: 0, Passed: 6497** (1 benchmark skipped via attribute, not an assertion skip); `dotnet build` 0 errors (4 pre-existing warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`, untouched).

### Minor note (not a defect)
`CalculateNxCurve` (public wrapper that maps `CalculateNx` over a threshold list) and `AssemblyStatistics.MedianLength` (upper-median, auxiliary, ASSUMPTION 2) are outside the unit's canonical N50/L50/Nx/auN contract and not asserted as canonical values. The curve wrapper is covered transitively by the `CalculateNx` tests; no source-defined behaviour is left unverified.

### Findings / defects
None. No code or test change required this session.

## Verdict & follow-ups
- **Stage A: PASS · Stage B: PASS · State: ✅ CLEAN · Test-quality gate: PASS.**
- No defects logged. Implementation faithfully realises the independently-sourced N50/L50/Nx/Lx/auN/gap definitions; tests lock exact sourced values; full suite green.
