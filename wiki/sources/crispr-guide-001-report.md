---
type: source
title: "Validation report: CRISPR-GUIDE-001 (CRISPR guide RNA on-target efficacy scoring — Doench-2014 Rule Set 1 + Doench-2016 Rule Set 2 / Azimuth, CrisprDesigner.CalculateOnTargetDoench2014 / CalculateOnTargetRuleSet2)"
tags: [validation, primer, governance]
doc_path: docs/Validation/reports/CRISPR-GUIDE-001.md
sources:
  - docs/Validation/reports/CRISPR-GUIDE-001.md
source_commit: 8ab783ae77cc9e5a6a05c3daefb18a3a4ad4a52d
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CRISPR-GUIDE-001

The two-stage **validation write-up** for test unit **CRISPR-GUIDE-001** (CRISPR guide RNA design —
extracting candidate gRNAs at PAM sites and scoring their **on-target efficacy**), area **MolTools**,
re-validated 2026-06-24. This is the *report* artifact that feeds one row of the [[validation-ledger]];
it records the validator's independent **verdict** on both the algorithm description (Stage A) and the
shipped code (Stage B), inside the wider [[validation-and-testing]] campaign. The two learned scoring
models, their provenance, invariants, and edge cases are synthesized in the concept
[[crispr-guide-rna-design]] (the MolTools reagent-design anchor for CRISPR guides, sibling to the
primer/probe units); [[test-unit-registry]] defines the unit. This session is an independent
re-confirmation of two prior fixes — the Doench-2014 "Rule Set 1" grounding (commit `129c2ca`) and the
Doench-2016 Rule Set 2 / Azimuth addition (commit `57730b9`).

## Verdict

**Stage A: PASS · Stage B: PASS · End state: ✅ CLEAN.** No defect found; no code or test change
required. The three CRISPR test classes (`CrisprDesigner_Doench2014_Tests`, `CrisprDesigner_RuleSet2_Tests`,
`CrisprDesigner_GuideRNA_Tests`) ran **54 passed, 0 failed** (build 0 warnings / 0 errors). This is an
independent re-validation from a fresh context: the **CRISPOR reference `doenchScore.py` was re-downloaded**
this session and the model re-grounded from scratch rather than trusting the repo's copied table.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CrisprDesigner.cs`:

- `CalculateOnTargetDoench2014(string context30Mer)` (`:521–568`; constants/table `:473–502`) — **Rule
  Set 1**, a logistic-regression linear model over a fixed **30-nt context** `[4 nt 5′]+[20 nt
  protospacer]+[3 nt PAM]+[3 nt 3′]`. `score = intercept`; a **GC term** over `seq[4..24)` adds
  `abs(10−gcCount)·gcWeight` with `gcWeight = gcLow if gcCount ≤ 10 else gcHigh`; a feature loop adds a
  `weight` when `CompareOrdinal(seq,pos,modelSeq,0,len)==0` (≡ `seq[pos:pos+len]==modelSeq`); output
  `1/(1+e^−score)` → (0,1), scaled ×100.
- `CalculateOnTargetRuleSet2(...)` — **Rule Set 2 / Azimuth**, a trained scikit-learn **GBRT** (not a
  coefficient table). Engine `AzimuthRuleSet2.cs` (sklearn-free reader + featurizer + tree traversal);
  embedded pickles `Resources/azimuth_rs2_{nopos,full}.bin`; extractor
  `scripts/azimuth/extract_azimuth_model.py`; externally-derived oracle CSVs
  `scripts/azimuth/oracle/{nopos,full}_oracle.csv` (947 rows).
- Tests: `CrisprDesigner_Doench2014_Tests` (9), `CrisprDesigner_RuleSet2_Tests` (13),
  `CrisprDesigner_GuideRNA_Tests` (heuristic evaluator) → 54 total.

## Stage A — description (algorithm faithfulness)

Grounded against **primary sources opened independently this session**: Doench et al. 2014 (*Nat
Biotechnol* 32:1262, PMID 25184501) for Rule Set 1; the CRISPOR reference `doenchScore.py` (Haeussler et
al. 2016, *Genome Biol* 17:148), **re-downloaded raw** and read for the constants/`params` directly; and
Microsoft Research **Azimuth** (BSD-3-Clause) for Rule Set 2 provenance.

- **Constants confirmed byte-identical to the reference:** intercept **0.59763615** (ref L30), `gcLow
  −0.2026259` / `gcHigh −0.1665878` (ref L31–32); GC term `abs(10−gc)·gcWeight` with boundary
  `gcCount ≤ 10 → gcLow` (ref L38–43); output `1/(1+exp(−score))`.
- **Coefficient table — full independent diff:** the reference `params` (via `exec` of the downloaded
  file) vs the repo's `DoenchParams` (via regex) → **ref 70 / repo 70, ordered exact match, symmetric
  difference ∅**. Including the reference's intentional quirks (`(24,'AG'/'CG'/'TG')` reuse the
  `(24,'A'/'C'/'T')` weights; `(26,'GT')=(27,'T')=0.11787758`).
- **Worked examples reproduced from scratch** (reference formula in Python over the repo's parsed table):
  `TATAGCTGCGATCTGAGGTAGGGAGGGACC` → 0.7130893 (ref 0.713089368, Δ 2.1e-08); `TCCGCACC…GGCGC` →
  0.0189838 (Δ 3.2e-09); M-003 `AAAAA…GGAAA` ×100 → 4.4338168085 (< 1e-9). The ~1e-8 deltas are
  float-print precision in the reference's quoted literals; the model reproduces exactly.
- **Edge-case semantics** (documented, stricter than the guard-free reference): wrong length / null /
  empty / non-ACGT throw; lowercase is upper-cased → identical value; an **NGG PAM guard** (offsets
  25–26 == `GG`) enforces SpCas9 specificity — an **input guard, not a scoring-model change**.
- **Rule Set 2 / Azimuth** re-confirmed at the level of provenance and oracle externality: both `.bin`
  model files + extractor present; the oracle CSVs carry both `ref_score` (verified reference prediction)
  **and** `upstream` (Microsoft Azimuth fixture) with an `agrees` flag — proving the test oracles are
  **externally derived, not read off the C# output** — plus the documented ~38 % upstream-fixture drift.
  The detailed sklearn-free recovery (1e-13 featurizer match, CPython-2.7 column order, bit-identical
  tree traversal) is accepted from the prior report; model files and tests unchanged.

## Stage B — implementation

- **Rule Set 1 code path** (`:521–568`) is line-for-line faithful to the reference (score init, GC term,
  feature loop, sigmoid ×100).
- **Rule Set 2 code path** is additive (reader + featurizer + GBRT traversal, thin wrappers); no existing
  signature or test changed; nopos/full wrappers delegate to one engine.
- **Cross-verification vs the actual C#:** 54 passed / 0 failed. M-001/M-002 lock the reference
  worked-example scores ×100 (71.3089368437 / 1.89838463593, tol 1e-4); M-003 locks the
  independently-recomputed 4.4338168085 — all three reproduced by hand this session.
- **Test-quality audit (HARD gate): PASS, no green-washing.** M-001/M-002/M-003 are exact
  externally-sourced cross-checks; a wrong intercept (0.59763615 → 0.69763615) shifts the score ~2 points
  ≫ tol and fails them. Edge tests cover wrong length (29/31/short), null/empty, non-ACGT (`N`),
  lowercase equivalence, range/ordering, and non-NGG PAM rejection. Rule Set 2 tests assert against the
  verified reference and the upstream-agreeing subset (counts locked), never against C# output.
- **Numerical robustness:** sigmoid bounded; GBRT is a finite-tree sum; no overflow/precision concern.

## Findings

- **No code defect, no test change (State CLEAN).** Rule Set 1 layout (4+20+3+3), intercept, gcLow/gcHigh,
  GC-term boundary, output transform, and the full 70-entry coefficient table are byte-identical to the
  re-downloaded CRISPOR reference; three worked-example scores reproduced to ≤ 2e-8. Rule Set 2 / Azimuth
  provenance and oracle externality re-confirmed.
- **No follow-ups.** The full unfiltered suite was not re-run because no code changed.
