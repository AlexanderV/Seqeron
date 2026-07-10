---
type: source
title: "Validation report: ANNOT-CODONUSAGE-001 (Relative Synonymous Codon Usage — RSCU)"
tags: [validation, annotation, governance]
doc_path: docs/Validation/reports/ANNOT-CODONUSAGE-001.md
sources:
  - docs/Validation/reports/ANNOT-CODONUSAGE-001.md
source_commit: 987ea6c1cf04c61c6257f0034ea4d51e00e0fffc
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation report: ANNOT-CODONUSAGE-001

The two-stage **validation write-up** for test unit **ANNOT-CODONUSAGE-001** (Relative Synonymous
Codon Usage — RSCU), validated 2026-06-15. This is the *report* artifact that feeds one row of the
[[validation-ledger]]; it records the validator's **verdict** on both the algorithm description and
the shipped code. The measure itself is summarized in [[relative-synonymous-codon-usage]]; the
two-stage methodology is the [[validation-protocol]]. Distinct from the pre-implementation
[[annot-codonusage-001-evidence]] artifact.

## Verdict

**Stage A: PASS · Stage B: PASS-WITH-NOTES · End-state: ✅ CLEAN.** **No algorithm or code defect** —
the implementation faithfully realises the validated RSCU formula. The only issues were three
Stage-A invariants / edge branches not directly exercised by the original tests; all three were
closed in-session with sourced assertions (test-only, **zero code change**). Full unfiltered suite
**6568 passed / 0 failed** (1 pre-existing unrelated skipped benchmark), `dotnet build` 0 warnings /
0 errors. Fixture grew 11 → 13 tests. No follow-ups.

## Canonical method

`GenomeAnnotator.GetCodonUsage(IEnumerable<string>)` and its `GetCodonUsage(IEnumerable<string>,
GeneticCode)` overload, at `Seqeron.Genomics.Annotation/GenomeAnnotator.cs:922-992`. The standard
overload delegates to the `GeneticCode` overload. Codon counts are pooled across all input sequences
via case-insensitive (`ToUpperInvariant`) in-frame triplet stepping (`i += 3`, stopping at
`Length-3` so a partial trailing codon is dropped); only A/C/G/T codons are counted. Families are
built from `code.CodonTable` (RNA-keyed, `U→T` conversion), stop codons (`'*'`) skipped, giving
exactly the 61 sense codons. RSCU = `(double)nI * x / familyTotal`, with `familyTotal == 0 ⇒ 0.0`.
`ArgumentNullException` on null `codingSequences` / `code`. The codon table
(`GeneticCode.cs:230+`) matches NCBI table 1.

## Stage A — description (algorithm faithfulness)

- Formula `RSCU = n_i · x_{i,j} / Σ_j x_{i,j}` confirmed **verbatim** against the LIRMM / Rivals RSCU
  methods page (symbols, normalisation by the family total, multiplication by family size `n_i`, and
  the [0, n_i] range). The "observed / expected-under-uniform" phrasing (PMC2528880, cubar
  `est_rscu`) is the same quantity: expected = Σx/n_i, so observed/expected = n_i·x/Σx. Two
  independent reference implementations (cubar `est_rscu`, and the definition attributed to Sharp
  et al.) corroborate.
- Synonymous families and stop codons confirmed against **NCBI genetic-code table 1**: Leu / Ser /
  Arg are the degeneracy-6 families, Met = ATG and Trp = TGG are single-codon, stops =
  TAA / TAG / TGA are excluded (RSCU defined over the 61 sense codons).
- Edge semantics all sourced: single-codon amino acids ⇒ RSCU = 1.0 (n_i=1); range [0, n_i];
  **whole-family zero count ⇒ 0.0 per member** (the base RSCU denominator is undefined; the unit
  deliberately does **not** apply the CAI 0.5 pseudocount, since that is a CAI convention — Sharp &
  Li 1987 — not part of plain RSCU). This is a documented, sourced design choice, not a defect.
- Independent hand cross-check from the formula + NCBI families reproduced: Leu `CTTCTTCTGTTA` ⇒
  CTT = 3.0, CTG = 1.5, TTA = 1.5, others 0, Σ over family = 6.0 = n_i; Phe `TTTTTC` ⇒ 1.0 / 1.0;
  Met `ATGATG` ⇒ 1.0; Trp `TGGTGGTGG` ⇒ 1.0.

## Stage B — implementation (code review + cross-check)

- The code computes `(double)nI * x / familyTotal` per family — the exact `n_i·x/Σx`. Running the
  tests confirmed `Count == 61`, Trp / Met = 1.0, Leu values 3.0 / 1.5 / 1.5. A 8-row
  cross-verification table (M1–M6 + Trp + the 61-codon set) all matched to `.Within(1e-10)`, including
  pooled inputs (M4 `["CTTCTT","CTGTTA"]` = M1) and `ATGTAA` (ATG = 1.0; no stop key emitted).
- Variant/delegate consistency: the `GeneticCode.Standard` overload equals the default overload
  (delegation). The legacy raw-count `GetCodonUsage(string)` is a different method, out of scope,
  left unchanged.
- **Test-quality audit (HARD gate) PASS:** all expected values (3.0, 1.5, 1.0, 6.0) trace to the
  LIRMM formula + NCBI families and were recomputed by hand — no tautologies / code-echoes; exact
  `.EqualTo(...).Within(1e-10)` assertions, no ranges / AtLeast / Contains where exact values are
  known, no skips or weakened tolerances.

## Findings (all test-only, closed in-session)

- **Trp single-codon branch** untested (INV-02 covered Met only) — added
  `GetCodonUsage_TryptophanSingleCodon_ReturnsOnePointZero` (TGG = 1.0).
- **61-codon cardinality / stop exclusion** untested (INV-04) — added
  `GetCodonUsage_OutputContainsExactlyThe61SenseCodons_NoStops` (Count == 61, no TAA/TAG/TGA).
- **Empty-enumerable branch** — C2 tested only `new[]{""}`; extended to also assert
  `Array.Empty<string>()` ⇒ 61-codon all-zero map.

No code defect; the three coverage gaps are the entirety of the PASS-WITH-NOTES qualifier.
