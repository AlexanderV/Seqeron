---
type: index
title: "Ingestion backlog — docs/algorithms reconciliation + queued sources"
tags: [meta, coverage]
created: 2026-07-09
updated: 2026-07-17
---

# Ingestion backlog

Coverage reconciliation for `docs/algorithms/**` (kept **in scope** by the
[coverage exclude policy](SCHEMA.md#coverage-exclude-policy)) plus source batches
queued for ingestion. Generated during the 2026-07-09 lint pass; regenerate when
concept pages are added or algorithm docs change.

The **pending** rows are a real coverage gap that folds into the main per-algorithm
ingest campaign (the same campaign advancing the `docs/Evidence/**` files) — not a
separate effort. A pending algorithm doc is resolved when a concept page lists it in
`sources:`; at that point it moves to the covered table.

Status at generation: **227** algorithm docs covered-via-concept, **0** pending — the per-domain

## Covered via concept (done) — moved

The **Covered via concept (done)** table and the full chronological resolution-note
history now live in **[[backlog-covered]]** (split out 2026-07-18 to keep this hub
under the page-size cap). Status at that split: **227** algorithm docs covered-via-concept,
**0** pending.

## Pending (fold into the ingest campaign)

The per-domain pending tables live in **[[backlog-pending]]** to keep this hub under the page-size cap. A pending row is resolved when a concept page lists the algorithm doc in `sources:`, at which point it moves to the *Covered via concept* table above. **The per-domain algorithm-doc pending backlog is now fully closed — all domain tables in [[backlog-pending]] are empty.**

## Queued source batches (approved 2026-07-09)

Approved for ingestion in the 2026-07-09 lint triage; pending `/wiki:ingest`.

### Testing methodology checklists (10) — `docs/checklists/`

- `docs/checklists/01_PROPERTY_BASED_TESTING.md`
- `docs/checklists/02_METAMORPHIC_TESTING.md`
- `docs/checklists/03_FUZZING.md`
- `docs/checklists/04_MUTATION_TESTING.md`
- `docs/checklists/05_SNAPSHOT_TESTING.md`
- `docs/checklists/06_ALGEBRAIC_TESTING.md`
- `docs/checklists/07_ARCHITECTURE_TESTING.md`
- `docs/checklists/08_DIFFERENTIAL_TESTING.md`
- `docs/checklists/09_COMBINATORIAL_TESTING.md`
- `docs/checklists/10_CHARACTERIZATION_TESTING.md`

### Validation governance ledgers (4) — `docs/Validation/`

- `docs/Validation/FINDINGS_REGISTER.md`
- `docs/Validation/LIMITATIONS.md`
- `docs/Validation/VALIDATION_LEDGER.md`
- `docs/Validation/VALIDATION_PROTOCOL.md`

### MCP top-level docs (3) — `docs/mcp/`

- `docs/mcp/MCP_STATUS.md`
- `docs/mcp/README.md`
- `docs/mcp/traceability.md`

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
