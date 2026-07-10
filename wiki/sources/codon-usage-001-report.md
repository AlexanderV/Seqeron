---
type: source
title: "Validation report: CODON-USAGE-001 (raw codon-usage table + TVD comparison, CodonOptimizer.CalculateCodonUsage / CompareCodonUsage)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/CODON-USAGE-001.md
sources:
  - docs/Validation/reports/CODON-USAGE-001.md
source_commit: b0db43a8692338b09003a14ca3e87a130a9b5a63
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: CODON-USAGE-001

The two-stage **validation write-up** for test unit **CODON-USAGE-001** (Codon Usage Analysis —
the raw per-codon **count** table plus a distribution-comparison similarity), validated
2026-06-24, Area *Codon*. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's independent **verdict** on both the algorithm
description (Stage A) and the shipped code (Stage B); the wider campaign is
[[validation-and-testing]] and [[test-unit-registry]] defines the unit. The measure itself (raw
`Dictionary<codon,int>` counts + the Total-Variation-Distance similarity) is synthesized in the
concept [[codon-usage-comparison]]. Distinct from the pre-implementation evidence artifact
[[codon-usage-001-evidence]] (sourced from `docs/Evidence/`) — this page is the independent
two-stage re-validation verdict, sourced from `docs/Validation/reports/CODON-USAGE-001.md`.

## Verdict

**Stage A: PASS-WITH-NOTES · Stage B: PASS · End state: CLEAN.** No code defect; no code change
required. Build green. The unit filter ran **22/22 pass**; the full `Seqeron.Genomics.Tests`
suite reported **18208 passed, 0 failed**. The single Stage-A note is a **scope mismatch only**,
not a formula error (see below).

## Canonical methods & source under test

In `src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/CodonOptimizer.cs`:

- `CalculateCodonUsage(seq)` (`:634–652`) — canonical. `string.IsNullOrEmpty` → empty dict;
  `ToUpperInvariant()` then `Replace('T','U')`; `SplitIntoCodons` → in-frame triplets; tallies into
  a `Dictionary<string,int>` of **raw counts**. `Count(c) = |{ i : seq[i:i+3] = c }|`.
- `CompareCodonUsage(seq1, seq2)` (`:657–681`) — comparison. Per-codon frequency
  `f_i(c) = count_i(c) / total_i`; returns `1 − Σ_c|f₁(c) − f₂(c)| / 2` over the union of observed
  codons; returns `0` when the union is empty or either total is `0`. This is exactly
  **Total-Variation-Distance similarity** (TVD = ½·L¹ between two discrete distributions), value in
  `[0, 1]`.
- `SplitIntoCodons` (`:687–695`) — loop `for (i = 0; i + 2 < len; i += 3)` ≡ `i ≤ len−3`: complete
  in-frame codons from **frame 0 only**; trailing 1–2 nt dropped. Verified for len 3/5/9/192.
- Tests: `tests/.../CodonOptimizer_CodonUsage_Tests.cs` (22 instances).

## Scope note (the only Stage-A note)

The session prompt framed the unit around the **full Kazusa / EMBOSS-`cusp`** output set (per-codon
count, **frequency-per-1000**, **fraction within the synonymous family**, and **RSCU**). The
*authoritative* unit definition (per `tests/TestSpecs/CODON-USAGE-001.md`, the Evidence doc,
`ALGORITHMS_CHECKLIST_V2.md:1291`, and the named test file) is the **narrower pair**: raw counts +
TVD similarity. This unit does **not** compute per-1000 frequency, per-family fraction, or RSCU —
those richer columns belong to adjacent units:

- **RSCU** → CODON-RSCU-001 / `CodonUsageAnalyzer.CalculateRscu` ([[relative-synonymous-codon-usage]]).
- **Per-codon frequency / fraction tables** → CODON-STATS-001 / SEQ-CODON-FREQ-001
  ([[codon-stats-001-report|CODON-STATS-001]] / [[seq-codon-freq-001-evidence|SEQ-CODON-FREQ-001]]).

So Stage A is PASS-WITH-NOTES **purely** for this scope framing; the unit's own counts + TVD
definitions are mathematically correct.

## Stage A — description (algorithm faithfulness)

Confirmed against the **EMBOSS `cusp` manual** (Number = observed raw count; Frequency = per-1000;
Fraction = proportion within the synonymous set), the **Kazusa Codon Usage Database** row format
`[triplet] [aa] [fraction] [freq/1000] ([number])` (e.g. E. coli `UUU F 0.57 19.7 (101)`), and
**Wikipedia "Codon usage bias"**. The TVD similarity's four invariants (identity = 1, symmetry,
range [0,1], disjoint → 0) are genuine mathematical facts of TVD. Edge semantics all match standard
practice: empty/null → empty dict; incomplete trailing 1–2 nt dropped (Kazusa/cusp convention);
frame 0 only (input assumed in-frame CDS); DNA `T` treated as RNA `U`; TVD empty/one-empty → 0.

**Independent hand cross-checks (reproduced):**

- **Counts** `ATGGCTGCTTAA` → (`AUGGCUGCUUAA`, frame 0) codons AUG, GCU, GCU, UAA ⇒
  `{AUG:1, GCU:2, UAA:1}`, total 4.
- **TVD M9** `AUGAUGCCCUUU` vs `AUGUUUUUUCCC` → Σ|Δ| = 0.5 ⇒ sim **0.75**. ✓
- **TVD M7** `CUGCUGCUGCUA` vs `CUACUACUACUG` → Σ|Δ| = 1 ⇒ sim **0.5**. ✓
- **TVD S6** `AUGGCUAUG` vs `AUGUUUAUG` → Σ|Δ| = 2/3 ⇒ sim **2/3**. ✓

## Stage B — implementation

Formula realised correctly. `CalculateCodonUsage` emits raw integer counts with **no**
normalisation (correctly matching the unit's stated output, not the full Kazusa table);
`CompareCodonUsage` normalises each sequence by its **own** total, computes ½·L¹, returns `1 − TVD`.
`CompareCodonUsage` delegates to `CalculateCodonUsage` for both inputs, so frequencies derive from
the same counting path — consistent, no `*Fast`/delegate divergence. The recomputed
cross-verification table (Counts `AUGGCUGCU`→{AUG:1,GCU:2}; incomplete `AUGGC`→{AUG:1}; all-64→64
keys ×1; DNA→RNA `ATGGCTTAA`→no `T`; sum invariant total = 6; TVD 0.5 / 0.75 / 0.0 disjoint / 2/3;
empty & one-empty → 0) reproduced exactly, all confirmed by passing tests.

**Test-quality audit:** the 22 tests assert **exact** values (`Is.EqualTo(...).Within(1e-10)` for
TVD; exact ints for counts) across M1–M10 / S1–S6 / edge cases (null, too-short, single-codon,
mixed-case, T→U, sum invariant, symmetry, disjoint→0, empty/one-empty). Deterministic; no
tautological no-throw assertions.

## Findings

- **No code defect and no test change (End state CLEAN).** Every worked example reproduced exactly;
  the code faithfully realises the validated counts + TVD formulas; frame handling (frame 0,
  trailing-codon drop) matches the Kazusa/cusp convention.
- **Documentation-only follow-up (non-blocking):** the prompt's "frequency-per-1000 / per-family
  fraction / RSCU" emphasis should be re-scoped to the adjacent units (CODON-RSCU-001,
  CODON-STATS-001, SEQ-CODON-FREQ-001), not this unit.

## Sources & cross-checks

External oracles: EMBOSS `cusp` manual, Kazusa Codon Usage Database (format + readme), Wikipedia
"Codon usage bias", and standard Total-Variation-Distance (½·L¹) probability theory — full URLs in
`docs/Validation/reports/CODON-USAGE-001.md`. No contradictions surfaced.
