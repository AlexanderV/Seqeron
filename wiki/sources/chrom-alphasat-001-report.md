---
type: source
title: "Validation report: CHROM-ALPHASAT-001 (alpha-satellite monomer detection — ChromosomeAnalyzer.DetectAlphaSatellite / FindCenpBBoxes)"
tags: [validation, chromosome, governance]
doc_path: docs/Validation/reports/CHROM-ALPHASAT-001.md
sources:
  - docs/Validation/reports/CHROM-ALPHASAT-001.md
source_commit: c2d37439e969b3159523581e9415ead8a622e666
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CHROM-ALPHASAT-001

The two-stage **validation write-up** for test unit **CHROM-ALPHASAT-001** (alpha-satellite / alphoid
centromeric monomer detection), validated 2026-06-25. This is the *report* artifact feeding one row of
the [[validation-ledger]] — one of that board's "24 new campaign units" (listed there as
*α-satellite/HOR*, pending first validation); the two-stage methodology is the [[validation-protocol]]
under [[validation-and-testing]], and [[test-unit-registry]] defines the unit. The two validated
methods, their derived constants and the HOR/suprachromosomal detectors that surround them are
synthesized in [[centromere-analysis]] (the chromosome centromere/alpha-satellite anchor).

**Scope vs. [[chrom-cent-001-evidence]] / CHROM-CENT-001.** CHROM-CENT-001 covers the whole
`AnalyzeCentromere` surface plus all four additive detectors. This unit is the **narrow slice** that
validates only **monomer-periodicity + AT-richness + CENP-B-box detection** (`DetectAlphaSatellite` /
`FindCenpBBoxes`). Higher-order-repeat structure (`DetectHigherOrderRepeat`) and
**suprachromosomal-family / chromosome-specific HOR-family naming** are explicitly **out of scope** here
— family assignment is a documented **data-blocked boundary** (needs curated T2T/CHM13 reference HOR
libraries not in the repo).

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN.** No code defect. One Stage-B **test gap**
closed (non-ACGT input, see below). Full **unfiltered** `dotnet test Seqeron.sln -c Debug`:
**0 failed** (Seqeron.Genomics.Tests 18779 passed).

## Canonical methods & constants validated

`src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`:

- `DetectAlphaSatellite(string)` (`:630`) — too-short guard (`< 171+5+1`); **AT content over ACGT
  bases only**; **periodicity** = best self-similarity over the period window **[166, 176]** (fraction
  of bases identical to the base `period` upstream); CENP-B box count;
  `IsAlphaSatellite = periodicity ≥ 0.50 AND AT > 0.50`.
- `FindCenpBBoxes(string)` (`:692`) and `CountCenpBBoxes` (`:855`) → `MatchesIupac` (`:872`): **Y→C/T,
  R→A/G**, all other positions exact match, 0-based ascending positions.
- Sourced constants: `AlphaSatelliteMonomerLength = 171` (`:32`), `CenpBBoxConsensus =
  YTTCGTTGGAARCGGGA` (`:41`).
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/ChromosomeAnalyzer_AlphaSatellite_Tests.cs`.

## Stage A — description (algorithm faithfulness)

Confirmed against primary literature opened live:

- **Monomer length = 171 bp** — Willard 1985; Waye & Willard 1987; Hartley & O'Neill 2019 (PMC6121732,
  *"fundamental 171bp monomeric repeat units"*). ✓ matches `AlphaSatelliteMonomerLength`.
- **CENP-B box = 17-bp consensus `YTTCGTTGGAARCGGGA`** (Y=C/T, R=A/G) — Masumoto et al. 1989
  (*J Cell Biol* 109(4):1963, origin of the 17-bp box) and Kugou et al. (PMC4843215), which quotes the
  consensus verbatim (core = `TTCG…CGGG`). ✓ matches `CenpBBoxConsensus`.
- **AT-richness** — the alphoid monomer is widely described as AT-rich; the `AT > 0.5` above-balance
  gate is a reasonable threshold.

**Divergence noted (does NOT affect verdict):** PMC6121732 renders the box as
`5'-T/CTCGTTGGAAA/GCGGGA-3'`, which expands to a **16-bp** `YTCGTTGGAARCGGGA` — it drops one `T` after
the leading Y. This is a **typo in that review**; the canonical Masumoto 1989 / PMC4843215 form is the
17-bp `YTTCGTTGGAARCGGGA`, which is exactly what the code implements. No action needed.

Independent hand-/Python cross-check (no dependence on the C# code): consensus length 17; all four IUPAC
corners match at index 0 (`CTTCG…`, `TTTCG…`, `CTTCG…GCGGGA`, `TTTCG…GCGGGA`); `Y` rejects a leading
`A`; a fixed position-2 `T→A` breaks the match; a box after a 50-bp flank is reported at 0-based offset
50.

## Stage B — implementation

`DetectAlphaSatellite` realises the description exactly: periodicity is "fraction of bases identical to
the base `period` upstream"; AT content **excludes non-ACGT from the denominator**; CENP-B matching
implements Y/R degeneracy with exact match at fixed positions. `CountCenpBBoxes` and `FindCenpBBoxes`
**share `MatchesIupac`**, so `count == positions.Count` (verified: box-count test = 10 == FindCenpBBoxes
hits). Mixed-case input matches uppercase.

Cross-verification recomputed vs code (values trace to the literature consensus, not code echoes):

- **171-bp tandem array** (60 A + 40 T + 36 C + 35 G, 20 copies): periodicity@171 = **1.0**; AT content
  = **100/171 = 0.5847953…**. ✓
- **Box-carrying monomer** (77 A + 17-bp box + 77 T, 10 copies): **10** CENP-B boxes. ✓
- Negatives (random, AT-rich-non-repetitive, GC-rich 16-bp-period mismatch), empty / null / too-short,
  and the IUPAC wrong-base / fixed-position-violation cases all match.

**Gap found & fixed (Stage-B test gap, not a code defect):** the edge-case list included "non-ACGT" but
the fixture had no test exercising it. Added
`DetectAlphaSatellite_NonAcgtBases_ExcludedFromAtContentDenominator` — monomer 60 A + 40 T + 36 C +
30 G + **5 N** → AT content = **100/166** (N excluded from the ACGT denominator, hand-derived),
periodicity 1.0, still called alpha-satellite. This locks the defined non-ACGT semantics.

## Findings

- **No code defect.** The CENP-B consensus, monomer length, AT-richness gate, periodicity, and IUPAC
  degeneracy are all correct and faithfully implemented.
- **Test-quality gate PASS (HARD gate):** every public surface and Stage-A edge case is covered; the one
  non-ACGT gap was closed with a hand-derived exact-value test. Expected numbers are literature- or
  hand-derived, not code echoes.
- **Data-blocked boundary (unchanged, out of scope):** suprachromosomal-family / chromosome-specific
  HOR-family assignment requires curated T2T/CHM13 reference HOR libraries not present in the repo — see
  [[centromere-analysis]] for the full HOR/SF writeup and the non-redistributable-library limitation.
