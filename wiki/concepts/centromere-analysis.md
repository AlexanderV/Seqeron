---
type: concept
title: "Centromere analysis (arm-ratio class + alpha-satellite / HOR / SF)"
tags: [chromosome, algorithm]
mcp_tools:
  - analyze_centromere
  - arm_ratio
  - classify_chromosome_by_arm_ratio
sources:
  - docs/Evidence/CHROM-CENT-001-Evidence.md
  - docs/algorithms/Chromosome_Analysis/Centromere_Analysis.md
  - docs/algorithms/Chromosome_Analysis/Higher_Order_Repeat_Detection.md
  - docs/Validation/reports/CHROM-ALPHASAT-001.md
  - docs/Validation/reports/CHROM-CENT-001.md
  - docs/Validation/reports/CHROM-HOR-001.md
source_commit: b8c0053177f07ad0218835374e920df51029953a
created: 2026-07-09
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-cent-001-evidence
      evidence: "Test Unit ID: CHROM-CENT-001 ... Area: Chromosome ... Canonical Method: ChromosomeAnalyzer.AnalyzeCentromere"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:aneuploidy-detection
      source: chrom-cent-001-evidence
      evidence: "Both are Chromosome-analysis family units; Wikipedia Karyotype groups human chromosomes by centromere position, the same karyotype layer aneuploidy scores"
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-alphasat-001-report
      evidence: "CHROM-ALPHASAT-001 validation report; Area: Chromosome; Canonical methods ChromosomeAnalyzer.DetectAlphaSatellite / FindCenpBBoxes — the monomer-detection slice of this concept, validated as its own test unit"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-cent-001-report
      evidence: "CHROM-CENT-001 validation report (2026-06-26); Area: Chromosome; re-validated after additive AssignSuprachromosomalFamily; Stage A/B PASS, CLEAN, 18860 passed / 0 failed, no code changed"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: chrom-hor-001-report
      evidence: "CHROM-HOR-001 validation report (2026-06-25); Area: Chromosome; Canonical method ChromosomeAnalyzer.DetectHigherOrderRepeat — the multi-monomer HOR-periodicity slice of this concept, validated as its own test unit; Stage A/B PASS, CLEAN, 18780 passed / 0 failed, one non-ACGT test gap closed, no code defect"
      confidence: high
      status: current
---

# Centromere analysis (arm-ratio class + alpha-satellite / HOR / SF)

Analysing the **centromere** of a chromosome — the constriction that links sister chromatids and
splits the chromosome into a short **p** arm ("petit") and long **q** arm. This is the second ingested
unit of the **Chromosome-analysis** family and its centromere/satellite anchor, a sibling of
[[aneuploidy-detection]] (which scores whole-chromosome copy number). Future chromosome-family units
(telomere, arm-ratio as a standalone, synteny, GC-skew — the compositional-skew statistic
now synthesized on [[nucleotide-composition-skew]]) get their own concepts. Validated under test
unit **CHROM-CENT-001**; the validation record is [[chrom-cent-001-evidence]], [[test-unit-registry]]
tracks the unit, and [[algorithm-validation-evidence]] describes the artifact pattern. The two-stage
validation verdict is [[chrom-cent-001-report]] — Stage A/B PASS, ✅ CLEAN (18860 passed / 0 failed, no
code changed), re-validated 2026-06-26 after the additive `AssignSuprachromosomalFamily` (SF) surface,
whose unresolved `Sf1OrSf2Dimeric` branch is runtime-guarded to `Permissive` mode.

The canonical method is `ChromosomeAnalyzer.AnalyzeCentromere(...)`; four **opt-in, additive** detectors
were layered onto it over successive sessions without changing its defaults. The base method is a
**generic tandem-repeat-density heuristic** — its `AlphaSatelliteContent` is a repeat score, *not* an
alpha-satellite-specific measurement — using a sliding window with k-mer frequency, low GC-variability
as the discriminating feature, and repeat-content estimation at k=15.

## Arm-ratio classification (Levan 1964)

Centromere **position** is classified from the arm ratio **q/p** (long arm / short arm), following the
Levan, Fredga & Sandberg (1964) nomenclature exactly:

| Arm ratio (q/p) | Levan sign | Class |
|-----------------|-----------|-------|
| ≤ 1.7 | M, m | Metacentric |
| (1.7, 3.0] | sm | Submetacentric |
| (3.0, 7.0) | st | Subtelocentric |
| ≥ 7.0 | t | Acrocentric |
| p = 0 (ratio ∞) | T | Telocentric |

The boundary values 1.7 / 3.0 / 7.0 come directly from Levan table entries. **Caveat:** practical human
karyotype classes come from cytogenetic (microscopic) observation, not genomic coordinate ratios, so
Levan thresholds applied to sequence-derived positions may reclassify borderline chromosomes. Telocentric
chromosomes are not present in normal humans.

## Alpha-satellite DNA — the specific detectors

Human centromeric heterochromatin is built from **α-satellite (alphoid) DNA**: a ~**171 bp** monomer
tandemly repeated and organised into higher-order arrays. Four additive detectors resolve it:

- **`DetectAlphaSatellite` / `FindCenpBBoxes`** — alpha-satellite-specific signal from three defining
  molecular signatures: tandem **171-bp periodicity** (±5 bp tolerance; monomers diverge 10–40%),
  **≥0.50 self-similarity** (lower bound of 50–70% intra-array monomer identity), **AT-richness > 0.50**,
  and the **17-bp CENP-B box** (IUPAC consensus `YTTCGTTGGAARCGGGA`, Y=C/T R=A/G; present in only a
  subset of monomers). **No consensus monomer string is embedded** — detection is period + AT + box, so
  no alphoid sequence is invented. This monomer-detection slice (`DetectAlphaSatellite` /
  `FindCenpBBoxes` + the 171-bp and 17-bp-box constants) was independently re-validated — Stage A/B PASS,
  CLEAN, no code defect — as its own test unit **CHROM-ALPHASAT-001**; see the report
  [[chrom-alphasat-001-report]] (periodicity search window [166,176], AT over ACGT bases only,
  `IsAlphaSatellite = periodicity ≥ 0.50 AND AT > 0.50`; note the PMC6121732 16-bp CENP-B rendering is a
  typo, the code uses the correct 17-bp form).

- **`DetectHigherOrderRepeat`** — the **HOR** organisation: split the array into consecutive 171-bp
  monomers, align monomer-vs-monomer (`SequenceAligner.GlobalAlign` + `CalculateStatistics`, EMBOSS-style
  identity), and report the HOR **period** = the smallest block size *k* for which monomers *k* apart are
  **≥95% identical** (inter-HOR divergence <5%) across **≥90%** of the array. Period 1 = homogeneous
  single-monomer (1-mer) array, *not* a multi-monomer HOR. Reports period, unit length (k×171 bp), copy
  number (⌊monomers/k⌋), and mean inter- vs intra-HOR identity. The key contrast: **intra-HOR** monomers
  are only 50–70% identical (distinct monomer types) while **inter-HOR** copies are 97–100% identical.
  This HOR method was independently validated — Stage A/B PASS, ✅ CLEAN, no code defect (one non-ACGT
  test gap closed; 18780 passed / 0 failed) — as its own test unit **CHROM-HOR-001**; see the report
  [[chrom-hor-001-report]] (`InterHorMinIdentityPercent = 95.0` = "<5% divergence", period-consistency
  0.90, `HasHigherOrderStructure = period ≥ 2`; independent k=4/m=7 cross-check: period 4, unit 684 bp,
  copy 7, inter 100% / intra 64.91%). Distinct from the monomer-detection slice
  [[chrom-alphasat-001-report]] (which left HOR out of scope) and the whole-centromere unit
  [[chrom-cent-001-report]]. The result is a read-only `HorResult` record struct (the 7 fields above);
  `monomerLength` defaults to 171 and must be ≥ 1 (else `ArgumentOutOfRangeException`); the search
  returns the **smallest** qualifying period (a dimeric array is period 2, not 4 or 6) and a trailing
  partial monomer is silently dropped. Cost is `O(M² · L²)` worst case (`M` monomers, `L`≈171 bp) —
  up to `O(M²)` Needleman-Wunsch monomer pairs at `O(L²)` each — but every ordered pair `(min, max)` is
  aligned at most once via a memoisation dictionary, so multi-kb arrays run in milliseconds. The
  detector recovers HOR *structure* only; it does **not** decompose **cascading / nested HOR-of-HORs**
  (Rosandić 2024) — that and SF-family labelling are left to dedicated centromere tooling (HORmon,
  alpha-CENTAURI) with curated T2T/CHM13 reference libraries.

- **`AssignSuprachromosomalFamily`** — SF1–SF5 assignment backed by a **bundled CC0 Dfam reference**
  (ALR/ALRa = A-type, ALRb carries the CENP-B box = B-type). Best-match each monomer (≥60% identity
  gate) to confirm identity and A/B-type, take the HOR period, and assign the SF from period + A/B
  composition: **SF3** pentameric (period%5==0), **SF4** monomeric period-1 A-type, **{SF1, SF2}**
  dimeric A→B, **SF5** irregular A/B mix.

## Derived parameters (all source-backed unless flagged)

| Parameter | Value | Basis |
|-----------|-------|-------|
| Monomer length / tandem period | 171 bp | Willard 1985; Alkan 2007 (167–171 bp) |
| Period search tolerance | ±5 bp | 10–40% monomer divergence allowance |
| Min self-similarity to call a tandem array | 0.50 | 50–70% intra-array identity |
| Min AT content | > 0.50 | AT-rich alphoid monomer |
| CENP-B box | 17 bp, `YTTCGTTGGAARCGGGA` | Masumoto 1989 |
| Inter-HOR min identity | ≥ 95% (<5% divergence) | Rosandić 2024; Alkan 2007; McNulty 2018 |
| HOR-period consistency across array | 0.90 | arrays are homogeneous, repeated hundreds–thousands× |
| Min identity to call alpha-satellite | ≥ 60% | **ASSUMPTION** — conservative gate |
| SF3 ⇔ period a multiple of 5 | period%5==0 | **ASSUMPTION** — pentameric-ancestry proxy |

## Limitations

- The base `AnalyzeCentromere` `AlphaSatelliteContent` is a **repeat-density heuristic**, not
  alpha-satellite-specific — the specific signal lives in the additive detectors. A
  [[research-grade-limitations|research-grade]] simplification.
- **SF1 vs SF2 not separated** (both dimeric with the identical A→B box pattern), and diverged-pentamer
  SF3 arrays whose period isn't a multiple of 5 (e.g. DXZ1, period 12) are not tagged. Both need an
  **SF-resolved consensus monomer library** (J1/J2/D1/D2/W1–W5/M1/R1/R2) that is published only in
  non-machine-retrievable supplements or unlicensed HMM repositories — **not redistributable**. Callers
  holding such a set can pass it to `AssignSuprachromosomalFamily(sequence, reference)`.
