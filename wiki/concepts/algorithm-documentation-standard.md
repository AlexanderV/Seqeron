---
type: concept
title: "Algorithm documentation standard (spec-doc template)"
tags: [documentation, validation, methodology]
sources:
  - docs/templates/algorithm-doc-template.md
source_commit: fee738e3f9d25fb455df96903ed5c5e8e2080cc2
created: 2026-07-18
updated: 2026-07-18
graph:
  relationships:
    - predicate: relates_to
      object: concept:algorithm-validation-evidence
      source: algorithm-doc-template
      evidence: "7.3 Related Tests, Evidence, or Documents ... Evidence: [<TEST-UNIT-ID>-Evidence.md](../../../docs/Evidence/<TEST-UNIT-ID>-Evidence.md); the spec doc's Test Unit ID header field and INV/ASM IDs are referenced from the Evidence artifact"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:definition-of-done
      source: algorithm-doc-template
      evidence: "Definition of Done requires a TestSpec and Evidence per unit; the algorithm-doc template is the third leg (the spec) whose Test Unit ID header field ties the three together"
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:validation-and-testing
      source: algorithm-doc-template
      evidence: "5.3 Conformance to Theory / Spec ... structurally enforced ... Implemented / Intentionally simplified / Not implemented — the doc records how the implementation was validated against the cited theory"
      confidence: medium
      status: current
---

# Algorithm documentation standard (spec-doc template)

The canonical structure every algorithm specification under `docs/algorithms/<group>/<file>.md` must follow, defined by `docs/templates/algorithm-doc-template.md`. It is the **spec** leg of the project's three-part documentation methodology — **spec → evidence → testspec** — captured by three sibling templates in `docs/templates/`: this `algorithm-doc-template.md`, the `evidence-template.md` (per-unit literature-traced source record, see [[algorithm-validation-evidence]]), and the `testspec-template.md`. The header's `Test Unit ID` field is the key that ties a spec doc to its Evidence artifact and TestSpec (see [[definition-of-done]]).

## Header metadata

Every spec doc opens with a five-row table:

| Field | Meaning |
|-------|---------|
| Algorithm Group | Area / subdomain |
| Test Unit ID | The `<UnitID>` (or `N/A`) tying to Evidence/TestSpec |
| Related Projects | Project names / `N/A` |
| Implementation Status | Exactly one of the five status values below |
| Last Reviewed | `YYYY-MM-DD` |

### Implementation Status vocabulary

Exactly one value, with a fixed meaning:

- **Production** — complete; matches the cited theory/spec for in-scope inputs; no known correctness gaps; tests cover the contract.
- **Simplified** — implements the core mechanism but omits parameters/refinements the cited source defines; gaps listed under 5.3 "Intentionally simplified".
- **Reference** — educational/illustrative; correctness prioritized over performance, scale, or full parameter coverage.
- **Framework** — provides API and extension points; users must supply domain data (gene sets, signature matrices, scoring tables) for production use.
- **Experimental** — under active development; contract/output may change without backwards compatibility.

## Required section order

Eight top-level (`##`) sections, kept in order; optional subsections that do not apply are removed rather than left empty:

1. **Overview** — 3–6 sentences: what it does, what problem it solves, when to use it, and whether it is exact, heuristic, probabilistic, or specification-driven.
2. **Scientific / Formal Basis** — theory only, no repository behavior. Subsections: `2.1 Domain Context`, `2.2 Core Model` (formula/recurrence/state machine/parsing model, cited), `2.3 Modeling Assumptions` (optional; `ASM-NN` IDs), `2.4 Properties and Invariants` (`INV-NN` IDs, stated so tests can check them), `2.5 Comparison with Related Methods` (optional).
3. **Contract** — `3.1 Inputs and Parameters`, `3.2 Output / Return Value`, `3.3 Preconditions and Validation` (indexing 0/1-based, coordinate inclusivity, alphabet, case, normalization, exception types).
4. **Algorithm** — `4.1 High-Level Steps`, `4.2 Decision Rules / Reference Tables / Data Structures` (optional; cite the origin of every numeric parameter table), `4.3 Complexity`.
5. **Implementation Notes** — repository code only, must not restate Section 2: `5.1 Location and Entry Points`, `5.2 Current Behavior`, `5.3 Conformance to Theory / Spec` (structurally enforced, below), `5.4 Deviations and Assumptions` (optional).
6. **Edge Cases and Limitations** — `6.1 Edge Cases` table, `6.2 Limitations`.
7. **Examples and Related Material** (optional) — `7.1 Worked Example`, `7.2 Applications`, `7.3 Related Tests/Evidence/Documents`, `7.4 Change History`.
8. **References** — full bibliographic citation + DOI/stable URL each; primary literature or official specs preferred.

## Enforced conventions

- **Structurally-enforced 5.3** — `5.3 Conformance to Theory / Spec` MUST contain three labeled bold blocks, each a bullet list using `(none)` rather than deletion: **Implemented (verbatim from the cited theory/spec)**, **Intentionally simplified** (feature → approximation → user-observable consequence), **Not implemented** (out-of-scope item → what users should rely on instead). A free-text paragraph here is non-conforming.
- **Citation discipline** — every formula, threshold, invariant, biological claim, or file-format rule carries an inline citation; the same fact is not repeated across Overview, Section 2, and Implementation Notes.
- **Identifier policy** — modeling assumptions use zero-padded `ASM-NN`, invariants use `INV-NN`; both are referenceable from tests, Evidence docs, and the Notes column of 5.4.
- **File-relative links** — Markdown resolves links relative to the current file; all in-repo links use file-relative paths (e.g. `../../../src/Foo/Bar.cs`), never repo-rooted. For docs at `docs/algorithms/<group>/<file>.md` the prefix to repo root is `../../../`.
- **References only** — use `## References`; do not create a separate `Sources` section.

## Multi-algorithm files

For files documenting peer methods together (UPGMA + NJ, Bray-Curtis + Jaccard, ESTIMATE + CIBERSORT, or triples), the hierarchy is fixed: Sections 2 and 4 stay single `##` umbrellas with per-algorithm `###` children (`2.A`, `2.B`, `2.C` …) using H4 **named, unnumbered** subsections (Domain Context, Core Model, Modeling Assumptions, Properties and Invariants); Section `5.3` stays a single `###` with `####` per-algorithm children each carrying the three required blocks. The letter→algorithm mapping is declared once as a quoted line under the Section 2 heading (e.g. `> A = UPGMA, B = Neighbor-Joining`). ID prefixes in combined files use the **algorithm token**, not the letter, so identifiers stay unique and readable from tests/evidence: `INV-UPGMA-01`, `ASM-NJ-02`. Sections 1, 3, 5.1, 5.2, 5.4, 6, 7, 8 are shared; a shared-table row is split only when genuinely algorithm-specific, prefixed with the algorithm token in brackets (e.g. `[UPGMA] negative branch lengths impossible`).

## Relationship to the rest of the methodology

This spec-doc template is one of three sibling documentation templates in `docs/templates/`. Its companions — `evidence-template.md` (the per-unit source record, whose shape [[algorithm-validation-evidence]] already describes) and `testspec-template.md` (the TestSpec required by the [[definition-of-done]]) — are expected to reconcile onto this page as they are ingested, so it becomes the shared home for the project's **documentation standard**. Together they operationalize [[scientific-rigor]] at documentation time and feed the [[validation-and-testing]] campaign.
