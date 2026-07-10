---
type: source
title: "Validation report: CHROM-TELO-001 (telomere analysis — TTAGGG/CCCTAA terminal-repeat detection + T/S-ratio length)"
tags: [validation, chromosome, governance]
doc_path: docs/Validation/reports/CHROM-TELO-001.md
sources:
  - docs/Validation/reports/CHROM-TELO-001.md
source_commit: 9dfe8fee4470a739dd91e9192efd5d7319ec5c50
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CHROM-TELO-001

The two-stage **validation write-up** (validated 2026-06-24) for test unit **CHROM-TELO-001** —
telomere analysis: detecting the conserved vertebrate telomere repeat at each chromosome end and
estimating telomere length directly from the detected run and from a qPCR **T/S ratio**. This is the
*report* artifact feeding one row of the [[validation-ledger]]; the two-stage methodology is the
[[validation-protocol]] under [[validation-and-testing]], and [[test-unit-registry]] defines the unit.
The repeat-detection/orientation rule, the 70 % per-window purity model, T/S linearity and the
invariants are synthesized in [[telomere-analysis]] (the chromosome telomere anchor). It is a
**distinct artifact** from the evidence file [[chrom-telo-001-evidence]] (that page catalogs the
literature sources/oracles from `docs/Evidence`; this page records the Stage A/B pass over
`docs/Validation/reports`).

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN.** No defect found; **no code changed this
session**. Build succeeded (0 warnings / 0 errors); the 33 Telomere tests
(`ChromosomeAnalyzer_Telomere_Tests.cs`) pass.

## Canonical surface validated

`src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`:

- `AnalyzeTelomeres(chromosomeName, sequence, telomereRepeat="TTAGGG", searchLength=10000, minTelomereLength=500, criticalLength=3000)`
  (`:250–288`) — empty guard, upper-casing, RC derivation via `DnaSequence.GetReverseComplementString`,
  5′/3′ windowing, threshold/critical flags.
- `MeasureTelomereLength` (`:293–340`) — tandem-window scan, per-window similarity ≥ 0.7, purity
  accumulation.
- `EstimateTelomereLengthFromTSRatio(tsRatio, referenceRatio=1.0, referenceLength=7000)` (`:345–352`).
- Constant `HumanTelomereRepeat = "TTAGGG"` (`:18`); pinned by a constant test.

## Stage A — description (algorithm faithfulness)

Confirmed against authoritative first-sources retrieved live:

- **Moyzis et al. (1988)** *PNAS* 85:6622 and **Meyne, Ratliff & Moyzis (1989)** *PNAS* 86:7049
  (PMC297991, PMID 2780561) — biotinylated `(TTAGGG)n` probes hybridise to the telomeres of all
  chromosomes in **91 vertebrate species** (bony fish → mammals; common ancestor > 400 Mya). `TTAGGG`
  is *the* 6-bp vertebrate telomere repeat.
- **Reverse-complement strand** — the sourced head-to-head arrangement `5'(TTAGGG)n-(CCCTAA)m 3'`
  confirms **CCCTAA** is the C-rich complementary strand. Hand check: complement(TTAGGG)=AATCCC,
  reversed = **CCCTAA**. ✓
- **Cawthon (2002)** *NAR* 30:e47 (PMID 12000852) — the T/S ratio (telomere signal T vs single-copy
  gene S, relative to a reference) is **proportional** to average telomere length; the linear
  `length = referenceLength × (tsRatio / referenceRatio)` faithfully realises "proportional".
- **Species table** (Wikipedia / Evidence) — vertebrates `TTAGGG`, Arabidopsis `TTTAGGG`, Tetrahymena
  `TTGGGG`, Bombyx `TTAGG` — consistent with the configurable `telomereRepeat` parameter.

**Orientation & tandem semantics (in-spec):** 5′ terminus modelled with the RC `CCCTAA` scanned forward
from index 0; 3′ with forward `TTAGGG` scanned backward from the terminus (correct for a single-strand
5′→3′ reference where the true 5′ end reads CCCTAA). Search is restricted to a `searchLength` window from
each end; counting is contiguous, non-overlapping windows that **stop at the first window below the
similarity threshold**, so only tandem copies contiguous to the terminal run are counted (an internal
non-terminal motif is not). Length = accepted windows × repeatLen; purity = matchingBases / totalBases.

**Edge cases (all defined & sourced):** empty → no telomere, `IsCriticallyShort = true`; no
repeats/homopolymer → length 0, not detected; below `minTelomereLength` → `HasTelomere = false` (length
still measured); region shorter than the repeat → (0,0), no crash; divergent repeats → 70 % per-window
similarity (≥ 5/6 for a 6-mer), purity reflects the divergence.

**Independent hand cross-check (reproduced spec values):** `[1000×A]+(TTAGGG)×200` (len 2200) → 3′ length
1200, purity 1.0; `(CCCTAA)×200+[1000×A]` → 5′ length 1200, purity 1.0; `searchLength 600` on the 3′ case
→ 100 windows → 600 (window truncation); divergent `TTAGGA×200` → 5/6 per window → length 1200, purity
5/6; T/S {1.5, 0.5, 2.0}@ref 7000 → {10500, 3500, 14000}, referenceRatio 2.0 → 3500, T/S 0.0 → 0. ✓

**Stage A PASS** — repeat, RC strand and T/S proportionality confirmed against the cited primary sources;
worked examples reproduce the spec values.

## Stage B — implementation

- **Formula realised:** default motif `"TTAGGG"`, RC → `CCCTAA` for the 5′ end; 5′ searched with RC
  forward from start, 3′ with forward motif backward from terminus; tandem counting contiguous /
  non-overlapping / accept-while-similarity-≥0.7-else-break; end-region search
  `5'=sequence[..min(searchLength,len)]`, `3'=sequence[max(0,len-searchLength)..]`; length = windows ×
  repeatLen; purity = matchingBases/totalBases; case-insensitive (`ToUpperInvariant` on both); custom
  repeat length supported (`repeatLen = repeatUnit.Length`, e.g. 7-mer `TTTAGGG` → 5/7 ≥ 0.7); T/S =
  `referenceLength * tsRatio / referenceRatio`. ✓
- **Cross-verification table recomputed vs code (via tests):** (TTAGGG)×200@3′ → 1200/1.0; (CCCTAA)×200@5′
  → 1200/1.0; both ends ×150 → 900 each; A×1000 → 0/not detected; empty → not detected/critically short;
  50 reps (300) with min 500 → Has=false, len 300; divergent TTAGGA×200 → 1200/(5/6); ×2000 searchLen
  15000 → 12000; searchLen 600 → 600; custom TTTAGGG×150 → 1050/1.0; T/S {1.5,0.5,2.0,1@ref2,0} →
  {10500,3500,14000,3500,0}. All ✓.
- **Variant/delegate consistency:** two public methods only; both behave per spec; no `*Fast` variants.
- **Test-quality audit (PASS):** `ChromosomeAnalyzer_Telomere_Tests.cs` — **33 tests** (incl.
  parameterised) assert **exact** lengths, exact purity (1.0, 5/6), exact T/S values, name preservation,
  threshold/critical logic, case-insensitivity, custom 7-mer repeat, search-window truncation and
  invariants — not tautologies; covers every Stage-A edge case; a constant test pins
  `HumanTelomereRepeat == "TTAGGG"`.

**Stage B PASS** — the implementation faithfully realises the validated description; correct motif,
correct RC strand on 5′, strictly tandem terminal counting; all cross-check values reproduced.

## Findings & boundaries

- **No code defect; no code changed this session.** State ✅ CLEAN.
- **In-spec assumption (noted, not a defect):** the backward 3′ scan is **phase-anchored to the
  terminus**, so detection assumes the terminal repeat run ends exactly at the sequence end — the
  biological norm and the only case the spec/tests exercise.
- **Documented configurability (not biological constants):** `criticalLength` (3000), `minTelomereLength`
  (500), `searchLength` (10000) and `referenceLength` (7000) are engineering defaults (Cawthon's clinical
  values differ by assay). See [[telomere-analysis]] for the full purity/T-S/orientation writeup and the
  research-grade (repeat-run, not shelterin/t-loop) scope.
