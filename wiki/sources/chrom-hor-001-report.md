---
type: source
title: "Validation report: CHROM-HOR-001 (higher-order repeat (HOR) detection — ChromosomeAnalyzer.DetectHigherOrderRepeat)"
tags: [validation, chromosome, governance]
doc_path: docs/Validation/reports/CHROM-HOR-001.md
sources:
  - docs/Validation/reports/CHROM-HOR-001.md
source_commit: 26fb94a83697567ac130e7f9454ea3429b8ccf73
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CHROM-HOR-001

The two-stage **validation write-up** for test unit **CHROM-HOR-001** — **higher-order repeat (HOR)
detection**, the periodic multi-monomer arrays of alpha-satellite centromeric DNA — validated
2026-06-25. This is the *report* artifact feeding one row of the [[validation-ledger]] (one of that
board's "24 new campaign units," listed there as *α-satellite/HOR*, pending first validation); the
two-stage methodology is the [[validation-protocol]] under [[validation-and-testing]], and
[[test-unit-registry]] defines the unit. The HOR algorithm, its constants and the surrounding centromere
/ alpha-satellite / suprachromosomal-family detectors are synthesized in [[centromere-analysis]] (the
chromosome centromere/satellite anchor).

**Scope vs. sibling chromosome units.** [[chrom-alphasat-001-report]] / CHROM-ALPHASAT-001 validates the
narrow **monomer-detection** slice (171-bp periodicity + AT-richness + CENP-B box) and explicitly left
`DetectHigherOrderRepeat` **out of scope**; this unit is exactly that HOR slice — the *multi-monomer*
periodicity that sits one level above the monomer. [[chrom-cent-001-report]] / CHROM-CENT-001 covers the
whole `AnalyzeCentromere` surface (Levan arm-ratio + all four additive detectors) and only re-confirmed
HOR in passing; this report is the dedicated Stage-A/B validation of the HOR method itself.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN.** No code defect. One Stage-B **test-coverage
gap** closed this session (non-ACGT edge case). Suprachromosomal-family assignment remains a documented
**data-blocked boundary** (needs curated T2T/CHM13 HOR reference libraries — not attempted by this
method). Full **unfiltered** `dotnet test Seqeron.sln -c Debug`: **Seqeron.Genomics.Tests 18780 passed /
0 failed** (incl. the new test), 0 warnings on the changed test project.

## Canonical method & result

`DetectHigherOrderRepeat(string sequence, int monomerLength = 171)` →
`HorResult(HasHigherOrderStructure, MonomersPerUnit, HorUnitLengthBp, HorCopyNumber, MonomerCount,
MeanInterHorIdentity, MeanIntraHorIdentity)`, in
`src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs` (method `:751–835`, record
`:157`). Thresholds `InterHorMinIdentityPercent = 95.0` (`:728`) and `HorPeriodConsistencyFraction = 0.90`
(`:734`); per-pair identity via `SequenceAligner.GlobalAlign` + `CalculateStatistics` (`:773`), mean via
`MeanPairwiseIdentity` (`:841`). Tests: `ChromosomeAnalyzer_HigherOrderRepeat_Tests.cs`.

## Stage A — description (algorithm faithfulness)

Confirmed against primary literature opened live (values trace to the sources, not to code):

- **Monomer length 171 bp; HOR period = k monomers ⇒ HOR unit length = k × 171 bp.** McNulty & Sullivan
  2018 (PMC6121732, "fundamental 171 bp monomeric repeat units"; "HOR unit length is determined by where
  the next monomer shows nearly total identity to the first"); Rosandić et al. 2024 (PMC11050224, "a peak
  of period *n* corresponds to *n* × 171 bp"). ✓ maps to the code's `unitLengthBp = period × monomerLength`.
- **The defining identity ordering: inter-HOR ≥ intra-HOR.** Intra-HOR monomers are **50–70% identical**
  (distinct monomer types within one unit); inter-HOR copies "differ by only a few percent" i.e.
  **≥ ~95%** (Rosandić: "minimal divergence … typically less than 5%"). ✓ `InterHorMinIdentityPercent =
  95.0` maps exactly to "<5% divergence"; the invariant inter ≥ intra is the biological signature.
- **Copy number = ⌊monomer count / k⌋** — standard tiling count, consistent with "HOR units repeated …
  hundreds to thousands of times." ✓
- **Period-consistency 0.90** — period must hold across ≥90% of k-spaced pairs, encoding "largely
  uninterrupted … homogeneous array." Source-consistent.
- Corroboration: Willard 1985 / Alkan et al. 2007 (171-bp monomer, <5% inter-HOR divergence).

**Edge-case semantics (defined & sourced):** < 2 full monomers (incl. empty/null/single) → no periodicity
→ no HOR (period 1, NaN identities); a **homogeneous 1-monomer** array → smallest high-identity period is
k=1, which is **not** a multi-monomer HOR ⇒ `HasHigherOrderStructure = false` (the literature distinction
between a monomeric/homogeneous array and a multimeric HOR unit); trailing partial monomer → dropped by
the floor split; non-ACGT (e.g. N) → the aligner treats it as a non-matching residue, no special meaning.

**Stage A findings:** none. The stub TestSpec/report contained no incorrect claim; the description is
biologically correct and now backed by verbatim quotes.

## Stage B — implementation

- Splits into ⌊len/monomerLength⌋ full monomers, trailing partial dropped (`:762`). ✓
- Searches the **smallest** period k∈[1, monomerCount/2] where ≥90% of k-spaced monomer pairs are ≥95%
  identical (`:790–805`). ✓
- `copyNumber = monomerCount/period`, `unitLengthBp = period·monomerLength` (`:815–816`); inter-HOR
  identity = mean of monomers period apart (`:820–826`); intra-HOR = mean pairwise identity within the
  first unit, only for k≥2 (`:829–831`). ✓
- `HasHigherOrderStructure = period ≥ 2` — correctly excludes the homogeneous k=1 case (`:833`). ✓
- Numerical robustness: NaN only where defined (no inter pairs / k<2); no div-by-zero (guards `:763/826/851`);
  `monomerLength < 1` throws `ArgumentOutOfRangeException` (`:753`); mixed case folded via `ToUpperInvariant`
  (`:759`). ✓

**Independent cross-check (hand-computed, not repo-fixture echoes).** A k=4 unit repeated m=7× over a
disjoint-substitution background (symmetric difference 60 ⇒ Hamming intra = (171−60)/171 = **64.91%**,
inside the 50–70% band), exact copies ⇒ inter = 100%:

| Quantity | Hand-computed | Code | Match |
|---|---|---|---|
| HOR period (monomers/unit) | 4 | 4 | ✓ |
| HOR unit length (bp) | 4×171 = 684 | 684 | ✓ |
| Copy number | ⌊28/4⌋ = 7 | 7 | ✓ |
| Monomer count | 28 | 28 | ✓ |
| Mean inter-HOR identity | 100.0% | 100.0% | ✓ |
| Mean intra-HOR identity | 64.91% (Hamming) | 65.5% (aligner) | band ✓ |
| inter ≥ intra invariant | true | true | ✓ |
| Non-HOR control (10 divergent monomers) | no HOR | HasHOR=false, period 1, inter=NaN | ✓ |

The repo fixture (k=3/m=5) reproduces exactly: period 3, copy 5, unit 513 bp, 15 monomers, inter 100%,
intra 57.8947% = exact gapless Hamming (0 gaps). **Note on the intra value:** intra-HOR identity is a
*continuous aligner* measure — when the background is truly non-periodic and substitutions are tightly
banded the optimal alignment is gapless and aligner identity = Hamming exactly; residual short-period
autocorrelation lets the aligner insert a few opportunistic gaps and report a slightly higher value. That
is a **fixture property, not a code defect** — both land inside the sourced 50–70% band and every
structural quantity matches exactly. The custom `monomerLength=10` overload was also exercised (period 2 /
copy 5 / unit 20).

**Test-quality audit (HARD gate PASS):** the fixture builds monomers from a shared background with point
substitutions and **hand-derives** expected identity as (171−|A△B|)/171 — values trace to construction,
not code echoes. Covered: 3-mer HOR, dimeric HOR, monomeric-only → no HOR (self-validating: asserts all
pairs <95% first), homogeneous 1-mer → not HOR, inter ≥ intra invariant, trailing partial, empty/null/
single, invalid `monomerLength` throws, mixed case.

**Gap found & fixed this session (Stage-B test gap, not a code defect):** no explicit **non-ACGT** test (a
Stage-A required edge case). Added `DetectHigherOrderRepeat_NonAcgtTrailingPartialMonomer_IsIgnored`
(80-bp N tail dropped as a partial monomer ⇒ period/copy/inter unchanged; hand-derived). Logged in
FINDINGS_REGISTER.

## Findings & boundaries

- **No production-code defect.** The period search, ≥95%/≥90% thresholds, copy-number/unit-length
  arithmetic, the inter ≥ intra invariant and the k=1-not-a-HOR rule are all correct and faithfully
  implemented.
- **Documented boundary (not a defect):** suprachromosomal-family / chromosome-specific HOR **assignment**
  is not attempted — it needs curated HOR reference libraries (e.g. T2T-CHM13 HOR annotations).
  `DetectHigherOrderRepeat` correctly limits itself to period / copy number / unit length / inter-vs-intra
  identity. See [[centromere-analysis]] for the full HOR/SF writeup and the non-redistributable-library
  limitation.
