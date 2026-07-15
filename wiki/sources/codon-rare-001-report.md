---
type: source
title: "Validation report: CODON-RARE-001 (rare codon detection — FindRareCodons + %MinMax + Sherlocc clusters, CodonOptimizer)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/CODON-RARE-001.md
sources:
  - docs/Validation/reports/CODON-RARE-001.md
source_commit: 8ce0af79a29f9fbddc217026508abc6f2c572e61
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CODON-RARE-001

The two-stage **validation write-up** for test unit **CODON-RARE-001** (rare codon analysis —
detecting low-frequency codons in a CDS and the rare-codon clusters/runs that matter biologically),
validated **2026-06-25** in a fresh re-validation that expanded the unit from the single per-codon
method to **three** public methods. This is the *report* artifact that feeds one row of the
[[validation-and-testing|validation ledger]]; it records the validator's independent **verdict** on
both the algorithm description (Stage A) and the shipped code (Stage B). [[test-unit-registry]] defines
the unit; the algorithm, its threshold rule, its cluster methods, invariants and worked oracles are
synthesized in the concept [[rare-codon-analysis]]. Distinct from [[codon-rare-001-evidence]] — the
pre-implementation evidence artifact sourced from `docs/Evidence/` — this page is the independent
two-stage re-validation verdict sourced from `docs/Validation/reports/CODON-RARE-001.md`.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: ✅ CLEAN.** No code defect, no code change. The
two fixtures `CodonOptimizer_FindRareCodons_Tests` (21) + `CodonOptimizer_RareCodonClusters_Tests` (24)
ran **45 passed, 0 failed**; the full `dotnet test Seqeron.sln` reported 0 failed (Genomics 18787
passed). Every hand-computed number reproduced against the live library.

## Canonical methods & source under test

Three public methods in `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`:

- `FindRareCodons` (per-codon, `:663`; code path `:609–629`) — a codon is rare when its per-family
  relative fraction is **strictly below** a threshold (default **0.15**). Reports `(nt position = i×3,
  codon, translated AA, frequency)` per rare codon.
- `CalculateMinMaxProfile` (%MinMax, `:720`; formula `:771–799`) — Clarke & Clark sliding-window score.
- `FindRareCodonClusters` (Sherlocc RCC, `:825`; `:825–901`) — 7-codon window, ≥ 4 slow-of-7 → cluster.

Shared helpers `SplitIntoCodons` (`:687`, `i+2 < length` → drops the trailing partial codon) and
`TranslateCodon` (unknown → "X", fraction 0).

## Stage A — description (algorithm faithfulness)

Re-validated from **authoritative external first sources retrieved this session**, ignoring the repo's
own TestSpec/Evidence:

- **Metric** is the per-amino-acid relative synonymous **fraction** f/Σf_family (Kazusa "fraction"
  column; each synonymous family sums to ~1.0) — a **relative-within-family** criterion, matching the
  spec's "< 15% relative frequency within their synonymous codon family," **not** an absolute per-1000
  frequency. Distinct from the relative-adaptiveness `w = f/f_max` used for CAI in `CalculateCAI`.
- **Threshold** default **0.15**, comparison strict `<` (a codon exactly at threshold is not flagged).
- **Rare E. coli set** independently corroborated against **Kazusa** (species 316407 / 83333) and the
  rare-codon literature — the canonical RIL set **AGA, AGG (Arg), AUA (Ile), CUA (Leu)** plus CGA (Arg),
  i.e. the codons supplemented by tRNAs in BL21-CodonPlus(RIL) strains. Sources: Wikipedia codon usage
  bias; **Shu et al. 2006** (CUA Leu ~0.04, consecutive rare Leu → ~3-fold translation inhibition);
  **Sharp & Li 1987** (rare criterion applied per synonymous family).
- **Cluster methods confirmed verbatim** against the two primary papers: **Clarke & Clark 2008** ("Rare
  Codons Cluster", PLoS ONE 3(10):e3412) — 18-codon window, Xij / Xmax,i / Xmin,i / Xavg,i definitions,
  the `(ΣXij−ΣXavg)/(Xmax−Xavg or Xavg−Xmin)×100` metric, sign convention (%Min reported negative),
  range −100…+100; and **Chartier, Gaudreault & Najmanovich 2012** (Bioinformatics 28(11):1438–1445,
  doi:10.1093/bioinformatics/bts149) — seven-position window, ≥ 4 pause positions of 7.
- **Edge-case semantics** all defined & sourced: empty/null → empty; non-÷3 → trailing partial codon
  dropped; sequence shorter than window → empty; window<1 / minRare<1 → throws; DNA T→U internal;
  case-insensitive; unknown codon → fraction 0 (always flagged); single-codon AAs (Met AUG=1.00,
  Trp UGG=1.00) never flagged unless threshold > 1.0; organism-specific (AGA rare in E. coli, common
  0.48 in yeast).

## Stage B — implementation

Code faithfully realises all three formulas. `FindRareCodons` = `freq < threshold` on per-family
fractions with strict `<`, position `i×3`. `CalculateMinMaxProfile` (`:771–799`) builds
`sumXij / sumXavg / sumMaxDelta / sumMinDelta` and branches `%Max` vs `%Min` exactly per Clarke & Clark;
`Xavg` is the per-family arithmetic mean; single-codon / unknown / stop codons contribute 0 to both
numerator and denominator (**no NaN**); window slides one codon → `codonCount − windowSize + 1` windows.
`FindRareCodonClusters` (`:825–901`) applies the verbatim window=7, ≥ 4-of-7 rule; a codon is rare by
the **same** `freq < 0.15` criterion as `FindRareCodons`; overlapping/adjacent qualifying windows merge
into maximal clusters.

**Numeric cross-check (hand-computed, then confirmed vs live library).** E. coli K12 Arg family (Kazusa
316407): Xavg = 0.16667, Xmax = 0.40, Xmin = 0.02. `3× CGC` (w 3) → **+100.0000**; `3× AGA` →
**−86.36364**; `CUG(0.50)+AGA` (w 2) → **+36.47059** — all reproduced exactly. Cluster hand-check:
`10× CGC + 7× AGA` → `FindRareCodons` flags all 7 at nt 30/33/36/39/42/45/48; `FindRareCodonClusters`
returns **one** cluster `codons 7..16, RareCount=7` (a 7-wide window at start 7 already holds 4 rare);
isolated rare codons yield **no** cluster while per-codon detection still flags them — the exact
capability the campaign added.

**Test-quality audit (HARD gate — PASS).** 45 tests assert exact sourced values (positions, codons, AAs,
frequencies to 1e-10/1e-12, cluster boundaries 0..6 / 0..9 / 7..16, %MinMax −86.363636…/+100/+36.47…,
bounds invariant [−100,100], determinism) traced to Clarke & Clark / Chartier and the hand-math — not
code echoes; no no-throw tautologies. Every public method/overload and every Stage-A path covered.

## Findings

- **Description-overreach note (not a defect).** `FindRareCodons` performs **per-codon** detection only;
  the run/cluster capability lives in the opt-in `CalculateMinMaxProfile` / `FindRareCodonClusters`.
  Cluster detection being out of scope for the base method is its documented contract, not a bug.
- **Sherlocc per-position criterion simplified (documented).** The paper derives the slow/pause cutoff
  statistically (extreme-value fit over a 7-codon usage average, thresholds 13–18 tested); the
  implementation uses the simpler per-codon fraction threshold (default 0.15) — the same criterion as
  `FindRareCodons`. The published **window = 7, ≥ 4-of-7** rule is reproduced verbatim and the
  per-position criterion is fully parameterizable (`rareThreshold`). An intentionally lighter-weight,
  table-driven variant; not a defect.
- The %MinMax metric uses the per-family relative **fraction** (Kazusa column), consistent with Clarke &
  Clark's `Xij` and with `FindRareCodons`; `w = f/f_max` is used separately for CAI. Agree.
- Minor table-vs-live-Kazusa rounding (e.g. CGA 0.06 vs ~0.07) changes no rare/common classification at
  tested thresholds.

No contradictions with the evidence artifact [[codon-rare-001-evidence]] or the concept
[[rare-codon-analysis]]; this report is the independent two-stage re-validation verdict for the same unit.
