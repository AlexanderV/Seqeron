---
type: index
title: "Ingestion backlog — docs/algorithms reconciliation + queued sources"
tags: [meta, coverage]
created: 2026-07-09
updated: 2026-07-20
---

# Ingestion backlog

Coverage reconciliation for `docs/algorithms/**` (kept **in scope** by the
[coverage exclude policy](SCHEMA.md#coverage-exclude-policy)). Generated during the
2026-07-09 lint pass; regenerate when concept pages are added or algorithm docs
change. As of the 2026-07-20 lint, both the per-domain algorithm-doc backlog and
the queued source batches are **fully closed** — see the sections below.

The **pending** rows are a real coverage gap that folds into the main per-algorithm
ingest campaign (the same campaign advancing the `docs/Evidence/**` files) — not a
separate effort. A pending algorithm doc is resolved when a concept page lists it in
`sources:`; at that point it moves to the covered table.

Status at generation: **227** algorithm docs covered-via-concept, **0** pending — the per-domain

## Covered via concept (done) — moved

The full chronological resolution-note history lives in **[[backlog-covered]]**, and
the **Covered via concept (done)** table (227 algorithm docs → concept pages) lives
in **[[backlog-covered-table]]** (history split out 2026-07-18, table further split
out 2026-07-20, to keep pages under the size cap). Status: **227** algorithm docs
covered-via-concept, **0** pending.

## Pending (fold into the ingest campaign)

The per-domain pending tables live in **[[backlog-pending]]** to keep this hub under the page-size cap. A pending row is resolved when a concept page lists the algorithm doc in `sources:`, at which point it moves to the *Covered via concept* table above. **The per-domain algorithm-doc pending backlog is now fully closed — all domain tables in [[backlog-pending]] are empty.**

## Queued source batches (approved 2026-07-09) — CLOSED

The batches approved in the 2026-07-09 lint triage are all resolved (verified by the
2026-07-20 coverage check — every doc below is referenced by a wiki page):

- **Testing methodology checklists (10)** — `docs/checklists/01…10_*.md` ingested;
  each has a concept page (e.g. [[property-based-testing]]) plus a `wiki/sources/` summary.
- **Validation governance ledgers (4)** — `docs/Validation/FINDINGS_REGISTER.md`,
  `LIMITATIONS.md`, `VALIDATION_LEDGER.md`, `VALIDATION_PROTOCOL.md` ingested (covered by
  [[validation-and-testing]], [[test-unit-registry]], [[validation-ledger]],
  [[validation-protocol]], [[validation-findings-disposition]]).
- **MCP top-level docs** — `docs/mcp/README.md` ingested as [[mcp-tool-catalog]].
  `docs/mcp/MCP_STATUS.md` and `docs/mcp/traceability.md` are **coverage exclusions**
  per [SCHEMA](SCHEMA.md#coverage-exclude-policy) (live status ledger / generated
  traceability matrix), not ingestion targets.

## Notes

- `docs/algorithms/README.md` and `docs/algorithms/CANONICAL_MAP.md` are index/map
  docs, not algorithm units. `CANONICAL_MAP.md` is ingested as the source page
  [[canonical-algorithm-map]] (canonical-identity map: alias→canonical IDs, folder
  buckets, legacy baselines) — the identity counterpart to this coverage ledger.
  `README.md` is **resolved index-only (2026-07-16)** — a coverage exclusion, not a
  wiki page. It is purely navigational (section headers linking to the per-algorithm
  folders, every one of which is already synthesized by a concept in the *Covered via
  concept* table) plus a Canonicalization section whose content is already owned by
  [[canonical-algorithm-map]]. No distinct synthesis to capture, so no dedicated page.
- The `docs/Evidence/**` campaign (175 of 213 remaining) is the primary driver: each
  Evidence ingest typically creates or extends the concept that also covers the
  matching algorithm doc, clearing a pending row here as a side effect.
