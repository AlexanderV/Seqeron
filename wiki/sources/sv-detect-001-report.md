---
type: source
title: "Validation report: SV-DETECT-001 (structural-variant detection from paired-end mapping signatures)"
tags: [validation, structural-variant, governance]
doc_path: docs/Validation/reports/SV-DETECT-001.md
sources:
  - docs/Validation/reports/SV-DETECT-001.md
source_commit: 1848b38435fea02da3a3b741832a07b43dedbb42
ingested: 2026-07-11
created: 2026-07-11
updated: 2026-07-11
---

# Validation report: SV-DETECT-001

The two-stage **validation write-up** for test unit **SV-DETECT-001** (Structural Variant
Detection from Paired-End Mapping (PEM) signatures — span + orientation → SV type),
validated 2026-06-15 (Area: StructuralVar). This is the *report* artifact that feeds one
row of the [[validation-ledger]]; it records the validator's **verdict** on both the
algorithm description and the shipped code. The two-stage methodology is the
[[validation-protocol]]; the discordant-pair algorithm itself is summarized in
[[discordant-pair-sv-detection]]. Distinct from the pre-implementation
[[sv-detect-001-evidence]] artifact.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS-WITH-NOTES · State: ✅ CLEAN.** One real
algorithm defect — RF (everted) pairs mis-classed as concordant, with no Duplication
branch — found and **completely fixed** this session across code, tests, description,
Evidence and TestSpec (logged FINDINGS_REGISTER A14). Full unfiltered suite
**6485 passed / 0 failed** (baseline 6484 + 1 new test M9; the SV-DETECT-001 fixture = 15
passed), `dotnet build` 0 errors, no new warnings in the changed files.

Canonical methods: `StructuralVariantAnalyzer.ClassifySV`, `DetectSVs`,
`FindDiscordantPairs` (+ helpers `IsConcordantOrientation`, `ClusterDiscordantPairs`).

## Stage A — description (algorithm faithfulness)

- Sources retrieved this session: **Medvedev, Stanciu & Brudno 2009** (Nat Methods 6(11s),
  author PDF, text extracted with pypdf) — the PEM signature catalogue verbatim: deletion =
  span **greater** than insert size, insertion = span **smaller**, inversion = same-strand
  (flipped) pair, cross-chromosome linking = translocation, and the corner case that an
  insertion larger than the fragment is invisible to span and never recovers the inserted
  sequence; **BreakDancer README** — `-c` default **3** with bounds `mean ± c·std`, `-r`
  minimum support default **2**, SV codes DEL/INS/INV/ITX/CTX; **cureffi.org / SAM FLAG
  0x02** — FR is the proper orientation, "RF, FF or RR" are abnormal; **DELLY (Rausch 2012)**
  and **SVXplorer (Kumar 2020)** — FR cluster ⇒ deletion, **RF (everted) cluster ⇒
  tandem-duplication**, FF/RR ⇒ inversion.
- Formula confirmed: span cutoff `μ ± c·σ` (default c=3, BreakDancer); DEL/INS/INV/CTX
  signature mapping (Medvedev 2009). Inclusive bounds (discordant iff strictly outside) is a
  sourced convention.
- **DEFECT corrected in the description (Stage A note):** the original description
  (SV_Detection.md), TestSpec and Evidence claimed **RF is concordant** ("FR/RF both proper,
  point inward"). This is **wrong** for the short-insert FR library the unit models — RF is
  the basic **tandem-duplication** signature and is discordant (confirmed by DELLY, LUMPY/Manta,
  SVXplorer, and the unit's own cited cureffi/BWA source). The description, TestSpec and
  Evidence were corrected this session; DELLY + SVXplorer added as sources 5–6.
- **Stage A: PASS-WITH-NOTES** (one sourced RF-concordance correction; otherwise correct).

## Stage B — implementation (code review + cross-check)

- Code path: `src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/StructuralVariantAnalyzer.cs`
  — `FindDiscordantPairs` (L154), `IsConcordantOrientation` (L204), `ClassifySV` (L225),
  `DetectSVs` (L265), `ClusterDiscordantPairs` (L286). Span bounds `μ ± c·σ` with strict
  `<`/`>` (exactly-at-bound concordant), the `ClassifySV` order (inter-chr→INV→RF-DUP→DEL→INS→
  Complex), the min-support gate and the cluster sweep all realised correctly.
- **Defect found & fixed (the RF orientation bug):** before the fix,
  `IsConcordantOrientation` returned true for **both FR and RF**, so an RF same-chromosome
  pair within span bounds was never flagged discordant, and `ClassifySV` had no Duplication
  branch — `SVType.Duplication` existed in the enum but was **unreachable** from the PEM path
  (an everted pair fell through to ComplexRearrangement). **Fix:** `IsConcordantOrientation`
  → concordant iff **FR only** (`strand1=='+' && strand2=='-'`); added an RF branch to
  `ClassifySV` (`strand1=='-' && strand2=='+'` → `SVType.Duplication`). No other caller relied
  on the old behaviour.
- Cross-verification (μ=400, σ=50, c=3 ⇒ bounds [250,550]) recomputed vs the fixed code:
  DEL (span 5000), INS (span 100), INV (same strand), **DUP (RF/everted — discordant; was
  wrong before fix)**, CTX (inter-chr, and inter-chr precedence over orientation), FR
  concordant, both boundaries inclusive, and the min-support cases (3× DEL → 1 Deletion with
  SupportingReads=3; 1× DEL below minSupport 2 → empty) — all match sources.
- Test-quality audit (HARD gate) PASS: a **green-washed test** (S4
  `FindDiscordantPairs_RfOrientationWithinBounds_NotDiscordant`, which locked in the wrong
  behaviour) was rewritten to `…_IsDiscordantDuplication` (RF **is** discordant → Duplication,
  citing DELLY/SVXplorer). Added **M9** `ClassifySV_RfEvertedOrientationSameChr_ReturnsDuplication`
  covering the previously-unreachable branch. All tests assert exact `SVType` / counts /
  `SupportingReads`; no weakened assertions, widened tolerances, or skips.
- **Stage B: PASS-WITH-NOTES** — one real defect fixed; tests strengthened.

## Findings

- **State ✅ CLEAN.** One real algorithm defect (RF mis-classed as concordant; no Duplication
  branch) **completely fixed** across code, tests, description, Evidence and TestSpec —
  logged FINDINGS_REGISTER **A14**.
- Out-of-scope by design (documented, not defects): split-read insertions > fragment;
  copy-paste/interspersed duplications and overlapping FR+RF complex clusters fold to
  `ComplexRearrangement`; linear cluster sweep rather than BreakDancer connection scoring.

See the report at `docs/Validation/reports/SV-DETECT-001.md`.
