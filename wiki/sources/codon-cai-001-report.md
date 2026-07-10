---
type: source
title: "Validation report: CODON-CAI-001 (Codon Adaptation Index — CAI, CodonOptimizer.CalculateCAI)"
tags: [validation, annotation]
doc_path: docs/Validation/reports/CODON-CAI-001.md
sources:
  - docs/Validation/reports/CODON-CAI-001.md
source_commit: 01b6d4e55e6b0e58c9b6d0b08a8c4cf532287e5d
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CODON-CAI-001

The two-stage **validation write-up** for test unit **CODON-CAI-001** (Codon Adaptation Index — the
geometric-mean relative-adaptedness score of a coding sequence vs a reference codon-usage table),
validated 2026-06-24 with a **2026-06-25 re-validation** in a fresh context. This is the *report*
artifact that feeds one row of the [[validation-ledger]]; it records the validator's independent
**verdict** on both the algorithm description (Stage A) and the shipped code (Stage B), and the wider
campaign is [[validation-and-testing]]. The algorithm, its invariants, oracles and edge cases are
synthesized in the concept [[codon-adaptation-index]] (the CAI anchor in the codon-usage-bias family,
built on the [[relative-synonymous-codon-usage]] synonymous-family normalization); [[test-unit-registry]]
defines the unit. Distinct from [[codon-cai-001-evidence]] — the pre-implementation evidence artifact
sourced from `docs/Evidence/` — this page is the independent two-stage re-validation verdict.

## Verdict

**Stage A: PASS · Stage B: PASS · End state: CLEAN.** No code defect; no production code touched. The
CAI fixture (`CodonOptimizer_CAI_Tests.cs`) ran **34 passed, 0 failed**; the full
`Seqeron.Genomics.Tests` suite **18787 passed, 0 failed**. Every cross-check value was reproduced by
independent hand computation (Python) to **≤ 1e-10**. This is an independent re-validation from a fresh
context: the unit had been reset to ⬜ in the 2026-06-25 re-reset because the opt-in
`excludeSingleCodonAminoAcids` parameter was added during the limitation-elimination campaign; it is now
re-validated against externally re-fetched first sources.

## Canonical method & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`:

- `CalculateCAI(string codingSequence, CodonUsageTable table, bool excludeSingleCodonAminoAcids = false)`
  (`:473–504`) — `w = Math.Max(codonFreq / maxFreq, 1e-6)`; `logSum += Math.Log(w);
  return Math.Exp(logSum / count)` — the exact geometric mean `exp((1/L)·Σ ln w_i)`. Stop codons
  (`"*"`) are `continue`-skipped (excluded from `L`); no-data amino acids give `w = NaN` and are
  skipped; empty / null / all-stop → returns 0. `T→U` and `ToUpperInvariant()` normalization applied.
- Helper `CalculateRelativeAdaptiveness` (`:506–522`) and the **derived** `SingleCodonAminoAcids` set,
  built in the static ctor (`:131–144`) from `AminoAcidToCodons` groups of size 1 excluding `"*"` →
  `{M, W}` for the standard code (not hard-coded). When `excludeSingleCodonAminoAcids: true`, Met/Trp
  are `continue`-skipped before scoring (`:494`), exactly the Sharp & Li / Jansen exclusion.
- `CalculateCAI` is the **sole canonical entry** — `OptimizeSequence` calls it for original/optimized
  CAI (`:264`, `:307`); no divergent re-implementation.
- Tests: `tests/Seqeron/Seqeron.Genomics.Tests/CodonOptimizer_CAI_Tests.cs` (34 tests).

## Stage A — description (algorithm faithfulness)

Confirmed against **Wikipedia "Codon Adaptation Index"** (re-fetched 2026-06-25, verbatim
`w_i = f_i/max(f_j)` and `CAI = (∏ w_i)^{1/L}` = geometric mean), **Sharp & Li (1987)** *Nucleic Acids
Res.* 15(3):1281–1295 (PMID 3547335, the original CAI paper), and **Jansen et al. (2003)** (PMC2684136),
which quotes Sharp & Li verbatim that single-codon families (AUG, UGG) "should be excluded" because
"their corresponding w value will always be 1 regardless of codon usage bias of the gene." Formula,
geometric-mean form, stop exclusion, and the opt-in single-codon exclusion all check out; log base is
irrelevant (cancels in the exp/log pair; code uses natural log consistently).

The **former divergence D-A1 is resolved**: strict Sharp & Li/Jansen single-codon-AA exclusion is now
selectable via `excludeSingleCodonAminoAcids: true`, and the historical inclusive (`w=1`) behaviour is
the documented default — so Stage A is **PASS**, not PASS-WITH-NOTES. The two remaining deviations are
documented, bounded fallbacks for *partial custom tables* (a case Sharp & Li did not face) and are
mathematically benign: the **`1e-6` clamp** for a codon absent from the table whose AA still has present
codons (`f_max > 0`, avoids `ln 0 = -∞`), and the **`w = NaN` skip** for a no-data amino acid
(`f_max ≤ 0`), which drops the codon from `L`.

**Independent hand cross-checks** (Python, exact) reproduced against the in-code Kazusa E. coli K12
table: `CUAACU` → 0.17056…, `AGAAGG` → 0.07071…, `CUGCUA` → 0.28284…, 4×CUG+1×CUA → 0.60342…. The
single-codon-exclusion table: `AUGUGG` (Met+Trp only) → 1.0 inclusive / **0.0** exclusive (L=0);
`AUGCUACUA` → 0.18566… / **0.08**; `AUGUGGCUA` → 0.43089… / **0.08**; `CUGCUA` (no Met/Trp) → 0.28284…
under either flag (exclusion no-op). Zero-frequency clamp against a custom `Leu={CUG:1.0}` table:
`CUACUG` → w=[1e-6, 1.0] → **0.001**; `CUA` → **1e-6**; `UUUCUG` (Phe no-data, skipped) → **1.0**. All
match fixture assertions to ≤ 1e-10.

## Stage B — implementation

Every worked example was re-traced against the code and reproduced exactly. Stage-B cross-verification
table (hand value vs test assertion): `CUAACU`/E.coli 0.17056, `CUGCCGACC`/E.coli 1.0,
`CUGCCGACC`/Yeast 0.4084, `AGAAGG`/E.coli 0.07071, `CUGCUA`/E.coli 0.28284, `AUGCUGCCGACC`/Human 0.7657,
4×CUG+1×CUA/E.coli 0.6034 — all ✓. The formula is realised correctly (`w = f/f_max` with the `1e-6`
clamp; `exp(logSum/count)` geometric mean), stop and optional single-codon-AA exclusion both fire, and
the clamp/skip fallbacks behave as documented.

**Test-quality audit:** 34 tests (30 pre-existing + **4 added this session**), all asserting **exact
sourced numeric values** (Sharp & Li / Jansen / Kazusa hand computation), not no-throw / tautology /
code-echo; deterministic. Coverage spans empty/null, single Met/Trp/both, all-optimal, rare, range
invariant, organism specificity (E. coli / Yeast / Human), DNA/lowercase, stop exclusion (incl.
mid-sequence), geometric-mean sensitivity & monotonicity, hand-calculated, performance, incomplete-codon,
and the full exclusion mode. The **coverage gap closed this session**: the previous suite had no test
for the documented `1e-6` zero-frequency clamp nor the `NaN` no-data-AA skip. Four tests added using a
public partial table (`CreateCodonTableFromSequence("CUG", …)`):
`AbsentCodonWithPresentSynonym_ClampsWeightToEpsilon` (CUACUG→0.001),
`AllCodonsAbsentFromFamily_ClampsToEpsilon` (CUA→1e-6),
`AminoAcidWithNoFrequencyData_IsSkipped` (UUUCUG→1.0),
`AllCodonsHaveNoFrequencyData_ReturnsZero` (UUUUUU→0).

## Findings

- **No code defect, State CLEAN.** The only change this session was additive (4 edge-case tests); no
  production code was touched. Every cross-check value reproduced by independent hand computation
  (≤ 1e-10); 34 fixture tests + 18787 full-suite tests pass.
- **No open follow-ups.** The former D-A1 single-codon-AA divergence is resolved (opt-in exclusion, with
  the sourced Jansen 2003 quote confirming the rule) and the coverage gap on the clamp/skip fallbacks is
  closed. See [[codon-adaptation-index]] for the standing cross-page nuance (RSCU's "0.5 pseudocount"
  reference-table convention vs this implementation's `1e-6` score-time clamp — different value, different
  stage, both guarding `log(0)`; not a source contradiction).
