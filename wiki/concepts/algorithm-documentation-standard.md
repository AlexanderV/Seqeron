---
type: concept
title: "Algorithm documentation standard (spec-doc template)"
tags: [documentation, validation, methodology]
sources:
  - docs/templates/algorithm-doc-template.md
  - docs/templates/evidence-template.md
  - docs/templates/testspec-template.md
source_commit: f6d643e426d704f1898703b290b6ac2fd0e46060
created: 2026-07-18
updated: 2026-07-18
graph:
  relationships:
    - predicate: relates_to
      object: concept:algorithm-validation-evidence
      source: algorithm-documentation-standard
      evidence: "7.3 Related Tests, Evidence, or Documents ... Evidence: [<TEST-UNIT-ID>-Evidence.md](../../../docs/Evidence/<TEST-UNIT-ID>-Evidence.md); the spec doc's Test Unit ID header field and INV/ASM IDs are referenced from the Evidence artifact"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:definition-of-done
      source: algorithm-documentation-standard
      evidence: "Definition of Done requires a TestSpec and Evidence per unit; the algorithm-doc template is the third leg (the spec) whose Test Unit ID header field ties the three together"
      confidence: medium
      status: current
    - predicate: relates_to
      object: concept:validation-and-testing
      source: algorithm-documentation-standard
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

## Evidence-artifact template

The **evidence** leg — `docs/templates/evidence-template.md` — defines the shape of the per-unit Evidence artifact stored at `docs/Evidence/<TEST_UNIT_ID>-Evidence.md`. It is the literature-traced source record behind a unit's tests; the [[algorithm-validation-evidence]] hub owns the artifact pattern and the per-file index, while this page records the template's mandated structure alongside its spec-doc sibling. The file is keyed by the same `Test Unit ID` as the spec doc's header (plus the algorithm name and a **Date Collected**), so spec ↔ evidence ↔ testspec share one identifier.

Required sections, in order:

1. **Online Sources** — one `###` subsection per source, each carrying `URL`, `Accessed` date, and an **Authority rank (1–5, per the Authoritative sources policy)**, followed by **numbered** Key Extracted Points (definitions, recurrences, complexity, worked examples).
2. **Documented Corner Cases and Failure Modes** — grouped by source; only cases actually documented in authoritative sources are admitted (no invented edge cases).
3. **Test Datasets** — one table per dataset, each with its own `Source:` citation; the canonical worked example(s) with exact parameters and expected outputs serve as the oracle / differential-test fixture. Large inputs (sequences, matrices) are stored in a separate file and referenced, not inlined.
4. **Assumptions** — behaviour **not** confirmed by authoritative sources; every entry MUST carry a bold `ASSUMPTION:` prefix with its justification, keeping unverified claims visibly separated from sourced facts.
5. **Recommendations for Test Coverage** — tiered as **MUST / SHOULD / COULD Test**; each MUST item cites its backing `Evidence: <source>`, each SHOULD/COULD item gives a `Rationale`.
6. **References** — every reference MUST include a verifiable link (DOI or stable URL); a reference without a link is treated as incomplete.
7. **Change History** — dated entries, the first being `Initial documentation.`

### Conventions distinct from the spec-doc template

- **Sourced-vs-assumed split** — the Evidence artifact structurally separates what authoritative sources confirm (Online Sources, Corner Cases, References) from what they do not (the `ASSUMPTION:`-prefixed Assumptions), a discipline the spec doc does not itself impose.
- **Authority ranking** — each source is rank-scored 1–5, so downstream tests can weight evidence by source authority.
- **Test-coverage tiering** — MUST/SHOULD/COULD is the Evidence artifact's own vocabulary for prioritising which behaviours the unit's tests must exercise, each tied back to a source or rationale, feeding the [[validation-and-testing]] campaign and the [[definition-of-done]]'s "Evidence documented" gate.

## TestSpec template

The **testspec** leg — `docs/templates/testspec-template.md` — defines the shape of the per-unit TestSpec stored at `TestSpecs/<TestUnitID>.md`, the first of the [[definition-of-done]]'s six required criteria ("TestSpec created"). Its header carries the same `Test Unit ID` (plus **Area**, **Algorithm**, a **Status** checkbox — `☐ In Progress` → `☑ Complete`, a fixed **Owner** of *Algorithm QA Architect*, and **Last Updated**), so spec ↔ evidence ↔ testspec share one identifier. Where the spec doc states the contract and the Evidence artifact traces the literature, the TestSpec is the **executable-test plan and audit trail**: it turns the Evidence artifact's tiered recommendations into concrete, ID'd test cases and proves that existing tests actually cover them.

Seven required sections, in order:

1. **Evidence Summary** — a condensed carry-over from the Evidence artifact: `1.1 Authoritative Sources` (a table of `# / Source / Authority Rank / DOI or URL / Accessed`, every source requiring a verifiable link, DOI preferred over PubMed/PMC), `1.2 Key Evidence Points` (each cited to a source), `1.3 Documented Corner Cases` (from evidence only; if none, the explicit sentence "No authoritative sources explicitly specify corner cases for <X>." rather than invention), and `1.4 Known Failure Modes / Pitfalls` (each source-cited).
2. **Canonical Methods Under Test** — a `Method / Class / Type / Notes` table where **Type** is one of three fixed values setting test depth: **Canonical** (deep evidence-based testing), **Delegate** (smoke verification only — 1–2 tests proving delegation), **Internal** (tested indirectly via canonical methods).
3. **Invariants** — `INV-N` rows, each with a statement, a **Verifiable** Yes/No flag, and an **Evidence** cell that is either a source or the literal **ASSUMPTION**; the same `INV-NN` identifiers the spec doc's §2.4 declares, now made checkable.
4. **Test Cases** — the tiered case list, mapped one-to-one onto the Evidence artifact's **MUST / SHOULD / COULD** coverage tiers: `4.1 MUST Tests (Required — every row needs Evidence)` with columns `ID / Test Case / Description / Expected Outcome / Evidence` and `M1…` IDs whose Evidence cell is a source or **ASSUMPTION**; `4.2 SHOULD Tests (Important edge cases)` and `4.3 COULD Tests (Nice to have)`, both replacing the Evidence column with **Notes**. Complex inputs (sequences, matrices) are described briefly here with the full data kept in the test code or a referenced data file.
5. **Audit of Existing Tests** — the traceability and consolidation record, the section that distinguishes a TestSpec from a plain plan: `5.1 Discovery Summary` (where tests were found, file paths); `5.2 Coverage Classification` (a `Area/Test Case ID / Status / Notes` table with a fixed status vocabulary — **✅ Covered, ⚠ Weak, ❌ Missing, 🔁 Duplicate**); `5.3 Consolidation Plan` (the **Canonical file** that tests converge into, plus what to **Remove** and why); `5.4 Final State After Consolidation` (`File / Role / Test Count`); `5.5 Phase 7 Work Queue` (every ❌ and ⚠ from §5.2 listed and implemented, each marked **✅ Done** or **⛔ Blocked**, with a tally whose **Remaining must be 0**); and `5.6 Post-Implementation Coverage` (a re-audit where **ALL** rows must be **✅** — if any remain ❌ or ⚠, the Test Unit **cannot** be marked ☑ Complete).
6. **Assumption Register** — a total count plus a table of every **ASSUMPTION** referenced anywhere in the spec and where it is used; zero assumptions is the stated goal and each one is treated as a risk.
7. **Open Questions / Decisions** — unresolved items, with "None" permitted only if genuinely nothing is open.

### Conventions distinct from the two sibling templates

- **Pass/fail is structural, not prose** — the unit's completion gate is mechanical: §5.5 Remaining = 0 and §5.6 all-✅. This is the concrete enforcement behind the [[definition-of-done]]'s "Edge cases covered" and "Tests pass" criteria — the TestSpec cannot honestly close while an ❌/⚠ survives.
- **Evidence → test-case traceability** — every MUST test carries an Evidence cell (source or ASSUMPTION) and the tiers mirror the Evidence artifact's MUST/SHOULD/COULD recommendations, so each test traces up to a sourced behaviour and any unsourced test is visibly flagged. The Assumption Register (§6) aggregates those flags for review.
- **Existing-test audit and consolidation** — unlike the spec and evidence templates (which describe intended behaviour), the TestSpec inventories and rationalises the *actual* test files (discovery → classify → consolidate → work-queue → re-audit), making it the artifact that ties written tests back to the specified contract.

## Relationship to the rest of the methodology

This page is the shared home for the project's **documentation standard**, and it now captures all **three** sibling templates in `docs/templates/`: the spec-doc template (`algorithm-doc-template.md`, above), the `evidence-template.md` (the per-unit source record, under "Evidence-artifact template"; its runtime artifact pattern lives on [[algorithm-validation-evidence]]), and the `testspec-template.md` (the per-unit executable-test plan and audit trail, under "TestSpec template"). The documentation standard is complete: **spec → evidence → testspec**, keyed end-to-end by a single `Test Unit ID` and gated by the [[definition-of-done]]. Together they operationalize [[scientific-rigor]] at documentation time and feed the [[validation-and-testing]] campaign.
