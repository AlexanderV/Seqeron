---
type: source
title: "Validation protocol — the two-stage, one-session-per-unit methodology"
tags: [validation, testing, governance, methodology]
doc_path: docs/Validation/VALIDATION_PROTOCOL.md
source_commit: 4e36db1f9a54daf2800a5d6647b30778b5187a7b
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Validation protocol — the two-stage, one-session-per-unit methodology

The **methodology document** the whole validation campaign runs on: how each Seqeron test unit is
independently re-validated. Where the [[validation-ledger]] records *where each unit stands* and the
[[test-unit-registry]] defines *what a unit is*, this protocol (`docs/Validation/VALIDATION_PROTOCOL.md`)
defines **how a validation session is conducted**. Its premise: existing TestSpecs, Evidence docs, and
unit tests were *authored alongside the code* and may share its blind spots, so each algorithm is
re-validated against authoritative **external** sources in a **fresh context**, deliberately separate from
the implementer's — a second, adversarial pair of eyes rather than the same context marking its own work.

## Two ordered stages (A before B)

Validation of a unit proceeds in two ordered stages, and **B is only meaningful once A passes**:

- **Stage A — validate the description.** Confirm the biology/maths is correct *in the abstract*: the
  cited sources are authoritative and actually say what's claimed (open them; don't trust the label),
  every formula matches the source exactly (symbols, normalisation, units, log base, edge conventions),
  definitions/conventions are explicit and standard (coordinate base 0/1, inclusive/exclusive ends,
  strand, percentage vs fraction, ambiguity codes), edge cases have a *defined, sourced* expected
  behaviour (never "implementation-defined"), an **independent cross-check** exists (reference-tool
  output, a paper's worked example, a Rosalind dataset, or a hand computation, with exact numbers), and
  the stated invariants are genuinely true.
- **Stage B — validate the implementation.** Confirm the code faithfully realises the *validated*
  description: it computes the real formula (not an approximation), handles each Stage-A edge case,
  reproduces every cross-verification value against the actual code, keeps variants/`*Fast`/delegates
  consistent with the canonical, is numerically robust (no precision loss/overflow/div-by-zero on stated
  ranges), and its **tests are real** — asserting exact sourced literals, deterministic, covering the
  edge cases, not tautologies or "no-throw."

The cardinal rule: **never silently change code to match a description Stage A proved wrong.** If the
description is wrong, fix the description too — fixing the code to a wrong spec is itself a defect.

## One algorithm = one session (fresh context)

Each unit is validated in its **own session** so no reasoning from a previous unit contaminates the next.
The protocol ships a **per-session prompt template** (validate `{UNIT-ID}`; do Stage A then Stage B; use
peer-reviewed papers, reference textbooks, Wikipedia *with its cited primary refs*, reference
implementations such as Biopython/EMBOSS/samtools, and public datasets/Rosalind; cross-check ≥1
independent reference implementation or hand-computed dataset). This external-primary-source, independent
cross-check discipline is the same one behind the [[validation-findings-disposition|green-washing remedy]]:
re-derive the oracle independently from a retrieved primary source, lock the exact literal, and
mutation-check — tests must encode the algorithm's spec, not the implementation's own constants.

## Two completion end-states

Every finished session ends in **exactly one** of two states — no half-fixes:

| End-state | Meaning |
|-----------|---------|
| **✅ FIXED / CLEAN** | No defect found, or every defect *completely* fixed in-session (code corrected, tests added to lock the sourced values, `dotnet build` + the unit's tests pass). Algorithm fully functional. |
| **🔧 LIMITED** | A defect/gap could not be fully fixed here. The report must record precisely **why** (root cause, correct behaviour per source, why the fix is out of reach — needs a dataset, model, upstream API, or larger redesign) and **what** is missing to make it functional. |

Per-stage verdicts are **✅ PASS · 🟡 PASS-WITH-NOTES** (minor, documented divergence — a by-design scope
boundary, not a defect) **· ❌ FAIL · ⬜** not yet validated. A Stage-A FAIL stops the session: log it, do
not proceed to B until resolved.

## Report template and environment

Each session writes a report at `docs/Validation/reports/{UnitID}.md` from a fixed template (Stage-A
sources/formula/edge-case/cross-check/findings; Stage-B code-path/formula-evidence/recomputed-table/
variant-consistency/test-audit/findings; final verdict + follow-ups). The repo targets **net10.0** (SDK
at `~/.dotnet`, no system PATH); any code-touching session must leave `Seqeron.Genomics.Tests` building
and green — the green baseline was **4484 passed / 0 failed (2026-06-12)**. Phase-1 scope was the **86
implemented (☑) units**. The 86 original per-unit reports were consolidated into the ledger +
[[findings-register]] and archived in git (`git show cb113ce:docs/Validation/reports/{UnitID}.md`).

## Where this fits

- [[validation-ledger]] — the live pass/fail board this protocol produces; every row is two per-stage marks + an end-state, driven by these two-context sessions.
- [[test-unit-registry]] — defines the `{UNIT-ID}` units each session validates.
- [[validation-findings-disposition]] — triages the notes/defects a session surfaces; shares the independent-oracle-derivation + mutation-check remedy.
- [[validation-and-testing]] — the overall correctness strategy this per-unit campaign operationalizes.
