---
type: source
title: "Validation report: ASSEMBLY-STATS-001 (assembly statistics — N50/L50/Nx/Lx/auN, gaps)"
tags: [validation, assembly, governance]
doc_path: docs/Validation/reports/ASSEMBLY-STATS-001.md
sources:
  - docs/Validation/reports/ASSEMBLY-STATS-001.md
source_commit: d584af4da843a888434b5c54e7277e8f3085b085
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ASSEMBLY-STATS-001

The two-stage **validation write-up** for test unit **ASSEMBLY-STATS-001** (assembly statistics — the
contiguity/QC summary computed over a set of contig lengths: N50, L50, Nx, Lx, N90, L90, auN, plus
totals, largest/smallest, GC, and an N-run gap summary), validated 2026-06-15 in a fresh context. This
is the *report* artifact that feeds one row of the [[validation-ledger]]; it records the validator's
independent **verdict** on both the algorithm description (Stage A) and the shipped code (Stage B). The
metric definitions, the inclusive-≥ threshold, the auN formula and the published oracles are summarized
in the concept [[assembly-statistics]] (anchor of the assembly STATS family), and the wider campaign is
[[validation-and-testing]]. Distinct from [[assembly-stats-001-evidence]] (the pre-implementation
evidence artifact, sourced from `docs/Evidence/`) — this is the independent re-validation verdict.

Canonical methods under test (in
`src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/GenomeAssemblyAnalyzer.cs`):
`CalculateStatistics`, `CalculateNx` (3-arg core + 2-arg delegate), `CalculateN50`, `CalculateAuN`,
`FindGaps`, plus the trivial `CalculateNxCurve` wrapper.

## Verdict

**Stage A: ✅ PASS · Stage B: ✅ PASS · State: ✅ CLEAN · Test-quality gate: PASS.** No defect found;
no production code or test changed this session. The implementation faithfully realises the
independently-sourced N50/L50/Nx/Lx/auN/gap definitions and the tests lock exact sourced values. Full
unfiltered suite **Passed: 6497, Failed: 0** (1 benchmark skipped via attribute, not an assertion
skip); `dotnet build` 0 errors (4 pre-existing warnings in an unrelated `ApproximateMatcher` test file,
untouched).

## Stage A — description (algorithm faithfulness)

Theory re-checked against primary sources opened live this session, independent of the repo artifacts:

- **QUAST `quast_libs/N50.py`** — the verbatim `NG50_and_LG50` cumulative loop (`limit =
  reference_length·(100−percentage)/100`, stop when `s <= limit` with `s = total − cumulative`); the
  `au_metric` = sum of squared contig lengths / reference length; empty list → `N50_and_L50` returns
  `(None, None)`.
- **Wikipedia "N50, L50, and related statistics"** — N50 = "sequence length of the shortest contig at
  50% of the total assembly length"; L50 = the smallest **count** of contigs whose lengths sum to half
  the total; N90 defined at 90% (N90 ≤ N50). Worked examples: **Assembly A {80,70,50,40,30,20}, total
  290 → N50=70, L50=2**; **Assembly B {…,10,5}, total 305 → N50=50, L50=3**.
- **Heng Li (2020) auN blog** — Nx: "Contigs no shorter than Nx cover x% of the assembly";
  **auN = ∑ᵢ Lᵢ² / ∑ⱼ Lⱼ** (no sorting needed).
- **Miller, Koren & Sutton (2010) §1.2** — "the smallest contig … whose combined length represents at
  least 50% of the assembly" (the inclusive **≥** boundary).

**Formula check.** The Nx/Lx inclusive boundary was proven equivalent by algebra: QUAST's `s <= limit`
with `s = total − cumulative`, `limit = total·(100−x)/100` ⇔ `100·cumulative ≥ total·x`, which is
exactly the integer test the repo uses (`cumulative*100 >= totalLength*threshold`). Confirmed inclusive
"at least x%", largest-first. `auN = Σl²/Σl` confirmed by both lh3 and QUAST `au_metric`.

**Edge-case semantics (sourced).** Empty input: authoritative tools return `None`; the repo
deliberately returns all-zero `AssemblyStatistics` (Nx=Lx=0, auN=0) — a documented API-shape choice,
non-correctness-affecting. Inclusive boundary (cumulative exactly = x% selects that contig) — sourced.
Monotonicity **N90 ≤ N50, L90 ≥ L50** — sourced (Wikipedia). `FindGaps` coordinate contract (0-based
inclusive `[Start,End]`, `Length = End−Start+1`, maximal `N`/`n` runs) is an implementation contract,
not externally standardized — verified by hand.

**Independent cross-check (hand-computed this session, vs sources not code).**

| Input | Metric | Hand value | Source |
|---|---|---|---|
| A {80,70,50,40,30,20}, tot 290 | N50 / L50 | 80+70=150 ≥ 145 → **70 / 2** | Wikipedia |
| A, x=90 (261) | N90 / L90 | 240<261, 270≥261 → **30 / 5** | Wikipedia |
| A | auN | 16700/290 ≈ **57.586** | lh3 / QUAST |
| B {…,10,5}, tot 305 | N50 / L50 | 200 ≥ 152.5 → **50 / 3** | Wikipedia |
| {100,80,60,40,20}, tot 300 | auN | 22000/300 = **73.333…** | lh3 / QUAST |
| {50,50}, tot 100 | N50 / L50 | cum 50 = 50% (inclusive) → **50 / 1** | Miller / QUAST |
| "ACGTNNACGTNNNNNNACGT", minGap=5 | gaps | [4,5] len2 filtered; [10,15] len6 kept | hand |

All match the spec/Evidence and the live external sources. **Stage A = PASS** — no divergence.

## Stage B — implementation (code review)

Code path `GenomeAssemblyAnalyzer.cs`:

- `CalculateNx` core (L237-264): integer-exact inclusive test `cumulative*100 >= totalLength*threshold`
  — provably equal to QUAST's `s <= limit`; `long` accumulation avoids overflow.
- `CalculateNx(IEnumerable,int)` (L276-281) and `CalculateN50` (L294-297) delegate correctly (sort
  descending + total → core; N50 = `.Nx`).
- `CalculateAuN` (L321-332): `Σ((double)l*l)/total` = Σl²/Σl; empty/zero-total → 0.
- `CalculateStatistics` (L121-165): aggregates N50/L50 (x=50), N90/L90 (x=90), largest=`lengths[0]`,
  smallest=`lengths[^1]`, totals, GC over non-N bases, gap stats; field order matches the record.
- `FindGaps` (L341-369): sentinel-terminated scan over maximal `N`/`n` runs, yields 0-based inclusive
  `[Start,End]`, `Length = i − gapStart`, with the `minGapLength` filter.

**Formula realised correctly.** The inclusive boundary, the auN formula and the largest-first
cumulative selection all match the validated description exactly. The post-loop fallback in
`CalculateNx` is dead code for thresholds ≤ 100 (the full sum always satisfies the test) and harmless.
Variant/delegate consistency: 2-arg `CalculateNx`, `CalculateN50` and `CalculateStatistics` all produce
the same N50/L50 as the core overload on shuffled and contig-realised inputs (M7–M9).

**Test-quality audit (HARD gate: PASS).** File
`tests/Seqeron/Seqeron.Genomics.Tests/GenomeAssemblyAnalyzer_AssemblyStatistics_Tests.cs`, 23 `[Test]`
methods. Sourced expectations, not code echoes: exact `Is.EqualTo` for N50/L50/N90/L90/coords,
`Within(1e-10)` for auN/GC/gap%; values trace to Wikipedia/lh3/QUAST, not to current output. No
green-washing (no `GreaterThan`/ranges where an exact value is known; the M4 monotonicity `<=`/`>=` is
a genuine property test with exact bounds). Coverage spans all 6 canonical methods plus every Stage-A
branch (inclusive boundary, single contig, empty, all-N, leading/trailing/interior/multiple gaps,
minGap filter); both worked examples A and B present.

## Findings

- **No algorithm defect. No production code or test changed. End-state ✅ CLEAN.** Description, formulas
  and every edge case are independently confirmed against QUAST `N50.py`, Wikipedia, Heng Li (2020) and
  Miller 2010; tests assert exact source-traced values.
- **No contradictions** among the sources — all give the identical largest-first, inclusive-≥
  cumulative-threshold definitions, and QUAST's `au_metric` matches Heng Li's ΣL²/ΣL exactly.
- **Non-defect notes (out of the canonical N50/L50/Nx/auN contract, not asserted as canonical):**
  `CalculateNxCurve` (a wrapper that maps `CalculateNx` over a threshold list, covered transitively) and
  `AssemblyStatistics.MedianLength` (upper median, auxiliary, Assumption 2). The empty-input all-zero
  return (vs QUAST `None`) is Assumption 1 — API-shape only.
- **Follow-ups:** none.
