---
type: source
title: "Validation report: CODON-RSCU-001 (Relative Synonymous Codon Usage — RSCU + codon counting, CodonUsageAnalyzer.CalculateRscu / CountCodons)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/CODON-RSCU-001.md
sources:
  - docs/Validation/reports/CODON-RSCU-001.md
source_commit: e3c96b23abab4c437f04a5ba414ba4d07d96f11b
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CODON-RSCU-001

The two-stage **validation write-up** for test unit **CODON-RSCU-001** (Relative Synonymous Codon
Usage — RSCU — plus the supporting codon-counting operation), validated 2026-06-15. This is the
*report* artifact that feeds one row of the [[validation-ledger]]; it records the validator's
independent **verdict** on both the algorithm description (Stage A) and the shipped code (Stage B).
The wider campaign is [[validation-and-testing]]; the measure itself is synthesized in
[[relative-synonymous-codon-usage]] and [[test-unit-registry]] defines the unit. Distinct from the
pre-implementation [[codon-rscu-001-evidence]] artifact (sourced from `docs/Evidence/`), and distinct
from the sibling report [[annot-codonusage-001-report]] — which validates the *same measure* on a
**different method** (`GenomeAnnotator.GetCodonUsage`); this unit validates the
`CodonUsageAnalyzer` surface.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · End state: ✅ CLEAN.** **No algorithm or code defect** —
the implementation faithfully realises the validated RSCU formula. The only issue was a
**test-coverage gap**: two documented Stage-A branches (the absent-family `expected>0 ? … : 0` guard,
and the stop-codon 3-fold family) were not exercised by the original 16-test fixture. All three were
closed in-session with sourced, hand-computed exact-value tests (test-only, **zero code change**).
Fixture grew 16 → 19 tests. Full unfiltered suite **6526 passed / 0 failed** (1 pre-existing
unrelated skipped benchmark), `dotnet build` 0 errors, the one changed file builds warning-free.

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonUsageAnalyzer.cs`:

- `CalculateRscu(DnaSequence)` / `(string)`, core `CalculateRscuCore` (`:88–110`) — groups
  `CodonToAminoAcid` by amino acid, computes `expected = totalCount / numSynonymous`, then
  `rscu[codon] = expected > 0 ? observed/expected : 0`. This is exactly `n_i·x_j/Σx` (since
  `observed/(Σx/n) = n·observed/Σx`), with the absent-family guard returning 0.
- `CountCodons(DnaSequence)` / `(string)`, core `CountCodonsCore` (`:37–52`) — non-overlapping
  triplets from offset 0 (`i += 3`), trailing 1–2 bases ignored (`i+3 <= length`), `IsValidCodon`
  excludes any non-ACGT triplet.
- String overloads (`:29–35`, `:80–86`): null/empty → empty dict; `ToUpperInvariant()` before
  counting (case-insensitive). `DnaSequence` overloads use `ArgumentNullException.ThrowIfNull`; both
  overloads delegate to the same `*Core`.

## Stage A — description (algorithm faithfulness)

Formula `RSCU_{i,j} = n_i · x_{i,j} / Σ_k x_{i,k}` confirmed **verbatim** against the LIRMM / Rivals
RSCU methods page (retrieved this session), and corroborated by **GenomicSig** (CRAN) and **seqinr
`uco`** — all giving the same "observed / expected-under-uniform" quantity, no-bias value **1.00**,
bounds **[0, n_i]**, primary attribution **Sharp, Tuohy & Mosurski (1986)**. Edge semantics all
sourced: single-codon families (Met, Trp) ⇒ RSCU = 1.0 (degenerate n=1 case); absent family (0/0)
canonically undefined ⇒ repository convention returns **0** (cubar uses a pseudocount instead —
documented divergence that only affects families that never occur); stop codons grouped as one 3-fold
synonymous family (reference tools commonly *exclude* stops, but this cannot change any amino-acid
codon's RSCU — documented convention). **Stage A: PASS**, no material divergence.

**Independent hand cross-check** (from the sourced formula, not code output):
Leu `CTGCTGCTGCTA` (CTG×3, CTA×1) ⇒ RSCU(CTG)=6·3/4=**4.5**, RSCU(CTA)=**1.5**, Σ=6=n_i;
Phe `TTTTTTTTC` (TTT×2, TTC×1) ⇒ **4/3**, **2/3**; Phe equal `TTTTTC` ⇒ **1.0** each;
Met `ATGATG` ⇒ **1.0**; stop `TAATAGTGA` ⇒ **1.0** each; stop biased `TAATAATGA` (TAA×2, TGA×1) ⇒
TAA=**2.0**, TGA=**1.0**, TAG=**0.0**, Σ=3=n_i.

## Stage B — implementation

The group-by-amino-acid + per-family ratio computes the sourced formula exactly; the single-codon
degenerate case falls out naturally; stops form a 3-fold family; the absent-family branch returns 0.
All six hand-computed cases were verified by adding tests and running them — the code reproduces the
externally-derived values to **`.Within(1e-10)`**. String and `DnaSequence` overloads delegate to the
same `*Core` (S5/S6 confirm exact agreement; case-insensitivity confirmed).

**Test-quality audit (HARD gate) — PASS after fix.** Every MUST/COULD value (4.5, 1.5, 4/3, 2/3, 1.0,
family sum = 6, bounds [0,6]) traces to the independently-retrieved LIRMM/GenomicSig formula, not to
code output — a deliberately-wrong implementation (per-thousand normalisation, or per-codon families)
would fail M1/M2/M5. Exact assertions with `Within(1e-10)`; no weakened Greater/Contains/range
substitutes, no skips, no widened tolerances.

## Findings (test-only, closed in-session)

- **Absent-family branch** untested (Assumption #1) — added
  `CalculateRscu_AbsentFamily_ReturnsZeroForEveryCodon` (`ATGATG` ⇒ all Leu/Phe = 0).
- **Stop-codon 3-fold family** untested (Assumption #2) — added
  `CalculateRscu_StopCodonFamily_TreatedAsThreeFoldFamily` (`TAATAGTGA` ⇒ 1.0 each) and
  `CalculateRscu_StopFamilyBiased_ComputesFamilyRatioAndSumsToDegeneracy` (`TAATAATGA` ⇒
  2.0/1.0/0.0, Σ=3).

No code defect; the two coverage gaps are the entirety of the PASS-WITH-NOTES qualifier. No follow-ups.
