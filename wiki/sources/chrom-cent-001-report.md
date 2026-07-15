---
type: source
title: "Validation report: CHROM-CENT-001 (centromere classification + alpha-satellite suprachromosomal-family assignment)"
tags: [validation, chromosome, governance]
doc_path: docs/Validation/reports/CHROM-CENT-001.md
sources:
  - docs/Validation/reports/CHROM-CENT-001.md
source_commit: d0034a86c805e107791d627b1505f78d3bef288d
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CHROM-CENT-001

The two-stage **validation write-up** (validated 2026-06-26) for test unit **CHROM-CENT-001** —
centromere classification plus **α-satellite suprachromosomal-family (SF) assignment**. This is the
*report* artifact feeding one row of the [[validation-ledger]]; the two-stage methodology is the
[[validation-protocol]] under [[validation-and-testing]], and [[test-unit-registry]] defines the unit.
The Levan arm-ratio classification, the alpha-satellite / CENP-B / HOR / SF detection layers and their
derived constants and limits are synthesized in [[centromere-analysis]] (the chromosome
centromere/satellite anchor). It is a **distinct artifact** from the evidence file
[[chrom-cent-001-evidence]] (that page catalogs sources/oracles from `docs/Evidence`; this page records
the Stage A/B pass over `docs/Validation/reports`).

**Re-validation trigger.** Limitation-fix commit `887a9945` **ADDED** SF assignment
(`AssignSuprachromosomalFamily` + `LoadBundledAlphaSatelliteReference` + new types). The prior validation
was **SUPERSEDED**; this session re-validated the unit fresh, with **primary focus on the new SF surface**
and a confirmation pass over the pre-existing canonical surface (Levan, α-satellite/CENP-B, HOR).

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN.** No defect found; **no code or test changed this
session**. Full unfiltered `dotnet test Seqeron.sln -c Debug`: **Failed: 0** (Seqeron.Genomics.Tests
**18860** passed). 275 `~ChromosomeAnalyzer` tests green (13 new SF tests + centromere + mutation-killer
suites).

## Canonical surface validated

`src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`:

- **NEW (primary focus):** `AssignSuprachromosomalFamily(sequence, reference=null)` (`:1090`),
  `LoadBundledAlphaSatelliteReference()`, `ClassifyFamily` (`:1187`); types `SuprachromosomalFamily`,
  `AlphaSatelliteReferenceMonomer`, `SuprachromosomalFamilyResult`, `AlphaSatelliteBoxType`; bundled CC0
  Dfam reference `Resources/AlphaSatelliteReference.fasta` (ALR / ALRa / ALRb).
- **Pre-existing (confirmation):** Levan classification (`CalculateArmRatio`,
  `ClassifyChromosomeByArmRatio`, `DetermineCentromereType`), α-satellite / CENP-B detection
  (`DetectAlphaSatellite`, `FindCenpBBoxes`), HOR detection (`DetectHigherOrderRepeat`).
- Test files: `ChromosomeAnalyzer_SuprachromosomalFamily_Tests.cs` (13, NEW),
  `ChromosomeAnalyzer_Centromere_Tests.cs`, `ChromosomeAnalyzer_MutationKillers_Tests.cs`,
  `ChromosomeAnalyzerTests.cs`.

## Stage A — description (algorithm faithfulness)

Confirmed against authoritative first-sources retrieved live:

- **Dfam REST** — `DF000000029` **ALR** (len 171), `DF000000014` **ALRa** (172, ≈83.7% to ALR),
  `DF000000015` **ALRb** (169); all classified `…Satellite;Centromeric`. All three consensus strings
  match the bundled FASTA **byte-for-byte**. Dfam data is **CC0** (Storer, Hubley, Rosen et al. 2021,
  *Mobile DNA* 12:2) → genuinely redistributable.
- **McNulty & Sullivan (2018)** (PMC6121732) — SF taxonomy verbatim: SF1 dimeric (J1·J2), SF2 dimeric
  (D1·D2), SF3 pentameric (W1–W5), SF4 monomeric (M1, acrocentric p-arms + Y), SF5 irregular (R1·R2);
  and the A/B-box rule verbatim ("B-type monomers contain CENP-B boxes, while A-type … contain a binding
  site for pJα"; A-type = J1,D2,W4,W5,M1,R2; B-type = J2,D1,W1–W3,R1).
- **Shepelev et al. (2009)** *PLoS Genet* 5:e1000641 — SF taxonomy origin (twelve monomer types in five
  families) confirmed cross-source.
- **Masumoto et al. (1989)** — 17-bp CENP-B box `YTTCGTTGGAARCGGGA` (Y=C/T, R=A/G).
- **Levan, Fredga & Sandberg (1964)** *Hereditas* 52(2):201–220 — arm-ratio r = L/S, centromeric index
  ci = 100·p/(p+q); cut-points 1.7 / 3.0 / 7.0 over m/sm/st/a/T (re-confirmed, unchanged).

**A/B typing hand-checked:** scanning the Dfam consensus for `[CT]TTCGTTGGAA[AG]CGGGA` → **ALRb** matches
at 0-based **126** (`CTTCGTTGGAAACGGGA`) → **B-type**; **ALR / ALRa** no match → **A-type**. Reproduces
the sourced rule (B carries the box, A does not) and the bundled `boxtype` headers.

**SF-assignment rule vs. published families.** Only two reproducible signals are available from the CC0
reference — the **HOR period** and the **A/B-box composition** of one HOR unit. The mapping is a faithful
reduction of the Shepelev/McNulty taxonomy:

| Signal | Code rule | Published basis |
|---|---|---|
| period multiple of 5 | → SF3 | SF3 pentameric (W1–W5) |
| period 1, all A-type | → SF4 | SF4 = M1 monomeric A-type |
| period 2 (A+B dimer) | → {SF1, SF2} | SF1/SF2 both dimeric A→B (cannot separate from CC0) |
| irregular A/B, both box types | → SF5 | SF5 = R1·R2 irregular |
| otherwise | → Unknown | — |

Edge cases (sourced/defined): empty/null/shorter than one 171-bp monomer → not alpha-satellite,
`Unknown`, period 0. A monomer below the **≥60% identity** gate to the closest reference → not
alpha-satellite (gate justified by alphoid divergence — ~16% from consensus, 50–70% between monomer
classes; Waye & Willard 1987 / PMC6121732; random DNA measured ~52%, real monomers ~85–98%).

**Stage A PASS** — every SF fact, the A/B rule, the CENP-B box motif/position, the three Dfam sequences,
the CC0 licence and the Levan thresholds are confirmed against external first-sources. The SF1-vs-SF2 and
non-period-5-SF3 residual is an **honest open boundary, not a hidden defect**.

## Stage B — implementation

- **Reference integrity:** `AlphaSatelliteReference.fasta` ALR/ALRa/ALRb lines are byte-identical to the
  Dfam REST `consensus_sequence` fields (lengths 171/172/169, `N` positions preserved); `boxtype` headers
  (A/A/B) match the sequence-derived CENP-B typing. `LoadBundledAlphaSatelliteReference` parses the
  embedded FASTA, upper-cases, returns 3 records with correct Name/Accession/Sequence/BoxType.
- **Code path:** `AssignSuprachromosomalFamily` guards empty/null/<171; splits into `len/171` whole
  monomers (trailing partial ignored); best-matches each monomer to every reference via
  `SequenceAligner.GlobalAlign` → `CalculateStatistics().Identity`; accepts ≥60% identity; period from
  `DetectHigherOrderRepeat`; first-unit A/B pattern → `ClassifyFamily`. `ClassifyFamily` matches the
  Stage-A table exactly.
- **Independent probe (own inputs, not test echoes):** ALRa×8 → **Sf4** (period 1, [A], 95.5%);
  (ALRa+ALRb)×6 → **Sf1OrSf2Dimeric** (period 2, [A,B], 97.5%); pentamer 3B+2A ×6 → **Sf3** (period 5,
  84.7%); irregular array → **Sf5**; **random DNA (seed 1234, 400 bp) → false / Unknown / 51.9%** (genuine
  negative below the 60% gate); ALRb×8 (B-only) → **Sf5**. `FindCenpBBoxes`: ALRb→`[126]`, ALRa→`[]`,
  ALR→`[]`. All match sourced expectations and test assertions.
- **Pre-existing surface (confirmation):** Levan `ClassifyChromosomeByArmRatio` normalises to
  r = max(p/q, q/p), ≤1.7 m / ≤3.0 sm / <7.0 st / else a, p=0 → Telocentric (correct per Levan 1964 — the
  defect fixed the prior round remains correct). The α-satellite/CENP-B/HOR detectors are unchanged by
  `887a9945` (additive only); the **M-SF-8** test asserts `DetectAlphaSatellite`/`DetectHigherOrderRepeat`
  are byte-identical before/after an SF call — the additive contract holds.
- **Test-quality audit (HARD gate PASS):** all public surface covered (M-SF-1..8, S-SF-1/2, C-SF-1, edge
  cases incl. empty-reference→throws); assertions are **sourced, not code echoes** (A/B typing → PMC6121732
  + box@126; SF families → the dimeric/pentameric/monomeric/irregular taxonomy; identity → measured
  Dfam-derived values; negative → the 60% gate); deterministic (local RNG seeds, no shared-RNG hazard).

**Stage B PASS** — the code faithfully realises the description; reference byte-identical to CC0 Dfam;
A/B typing correct; SF rule reproduces published signals on worked examples plus a real negative case.

## Findings & boundaries

- **No code defect; no code or test changed this session.** State ✅ CLEAN.
- **Documented open boundary (acceptable, not a defect):** **SF1 vs SF2** (both dimeric A→B) and SF3
  arrays whose period is not a multiple of 5 (e.g. dodecameric DXZ1) are **not resolved** from the CC0
  reference — they need the SF-resolved consensus monomer library (J1/J2/D1/D2/W1–W5/M1/R1/R2), which is
  **not CC0/redistributable**. Recorded in `Resources/README.md`, `LIMITATIONS.md` §2, and the TestSpec
  residual. See [[centromere-analysis]] for the full HOR/SF writeup and the non-redistributable-library
  limitation.
- **Runtime enforcement (LimitationPolicy):** the guarded branch — the unresolved `Sf1OrSf2Dimeric`
  (SF1-vs-SF2) call — has minimum access mode **`Permissive`** (`LimitationCatalog`). Under the default
  `LimitationPolicy.DefaultMode = Moderate` it throws `SeqeronLimitationException` (this guarded branch is
  allowed only under `Permissive`; see `LIMITATIONS.md` › Runtime enforcement). This confirms
  CHROM-CENT-001 is one of the units **named in the operating-envelope doc** — its LIMITED end-state
  (dimeric arrays resolvable only to `{SF1, SF2}`) is guarded at runtime, an additive policy layer that
  leaves the validated contract and CLEAN verdict unchanged.
- **Historical note:** the prior round (2026-06-25, Levan focus) found and fully fixed a real defect in
  `ClassifyChromosomeByArmRatio` (diverged from Levan 1964; omitted Subtelocentric) and de-green-washed
  its tests; that fix remains in place and is re-confirmed here. The present round re-validates the unit
  after the additive SF-assignment addition, which does not touch the Levan surface.
