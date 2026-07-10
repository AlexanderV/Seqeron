---
type: source
title: "Validation findings register — full disposition of every campaign note"
tags: [validation, testing, governance]
doc_path: docs/Validation/FINDINGS_REGISTER.md
source_commit: 9710ae416a27c89fd4c0699c0a7631c0df408224
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation findings register — full disposition of every campaign note

The governance ledger that collects **every** note, limitation, and follow-up raised across
all 86 per-unit validation reports in `docs/Validation/reports/` and assigns each one **exactly
one disposition**. It is the "what did we decide to do about each finding" side of the
per-unit [[validation-and-testing|validation campaign]] — the campaign produces findings; this
register triages them. The triage mechanism itself is written up as
[[validation-findings-disposition]].

## The four disposition categories

Each finding lands in exactly one bucket (dated 2026-06-12, Seqeron.Genomics):

- **A — FIXED-NOW.** Non-radical: doc/comment/spec-text corrections or small safe code changes,
  applied in the pass. Dozens of rows (A1-A65). Many are **test-quality** fixes with *no algorithm
  defect* — e.g. replacing tautological/code-echo assertions with externally sourced exact values.
  A minority are real code defects caught during re-validation (e.g. A11 `LogRatioToCopyNumber`
  round-half-to-even to match CNVkit; A63 `ClassifyChromosomeByArmRatio` rewritten to Levan 1964
  thresholds; A64 MobiDB-lite `f+` including histidine; A65 miRNA `CountSeedSitesInUtr` greedy
  non-overlapping scan per Garcia 2011).
- **B — FEASIBLE -> IMPLEMENT.** Needs a real code change but is safely implementable with strict
  tests in its own context (B1-B5): PROSITE reject unsupported metachars, Fst locus-count throw,
  EMBL `""`-unescape, splice-predict invariant consistency, RNA-PARTITION McCaskill outside
  recursion. All marked done.
- **C — NOT-POSSIBLE (radical).** Originally "requires a redesign / public-API change / new model,
  documented not changed." Notably, **most C items were subsequently IMPLEMENTED as approved
  breaking/additive changes** (C1 Kraken taxonomy-tree + k-mer-LCA + RTL classifier; C2 N-ary
  multifurcating phylo trees; C3 unrooted-bipartition Robinson-Foulds; C4 progressive/guide-tree MSA;
  C5 circular restriction digest; C6 reverse-strand RBS; C7 all four published CRISPR scoring models
  — Doench 2014, MIT/Hsu, CFD, Rule Set 2/Azimuth; C8 EMBL remote refs / single-dot ranges). Each
  records the retrieved primary source, an independent cross-check, and mutation-verified tests.
- **D — BY-DESIGN.** The note documents correct, sourced intended behaviour; "fixing" it would be
  wrong. No change (e.g. oncology LOH merge semantics, HLA homozygous-loss label, RNA-MFE
  NNDB-vs-ViennaRNA parameter rounding).

## Recurring theme: green-washing detection

A large fraction of category-A rows are **not** algorithm defects but **weak tests that "green-wash"
the code** — tautological assertions (`Significant == (PValue < 0.05)`), `Contains`/range checks that a
deliberately-wrong impl would still pass, or code-echo tests that recompute the expected value from the
implementation's own constants. The fix pattern is uniform: re-derive the oracle **independently** from
an externally retrieved primary source (paper, reference implementation, verbatim data file), lock the
exact literal, and **mutation-check** that a seeded bug now fails. This is the same "tests must encode
the algorithm/business spec, not a buggy impl" discipline the [[build-quality-gate]] applied to the
S1244/S125 Sonar giants.

## Status and the re-validation reset

The document opens with a **superseding banner**: on **2026-06-24** a full re-validation reset moved all
units back to pending in `docs/Validation/VALIDATION_LEDGER.md` and `docs/checklists/` for a fresh
end-to-end re-verification (triggered by the extensive code churn of the limitations-elimination campaign).
The dispositions below the banner are therefore **retained as historical evidence but SUPERSEDED** until
each unit is re-validated. Read this register as a snapshot of the 2026-06-12 pass plus its per-row
completion notes, not as live status — the ledger is ground truth for current state.

## Where this fits

- [[validation-and-testing]] — the correctness/validation strategy this register operationalizes finding-by-finding.
- [[validation-findings-disposition]] — the triage governance process (the A/B/C/D mechanism) abstracted from this document.
- [[build-quality-gate]] — shares the "review, don't green-wash" / spec-not-impl standing rule.
- [[test-unit-registry]] / [[definition-of-done]] — the per-unit tracking and bar that produce these findings.
