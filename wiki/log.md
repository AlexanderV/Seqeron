# Wiki Log

Append-only chronological record of operations on the wiki. Each entry begins with `## [YYYY-MM-DD] <op> | <description>` so it's parseable with `grep "^## \[" log.md | tail -N`.

Operations:
- `ingest` — a source was processed into the wiki.
- `query` — a question was answered against the wiki (typically only logged when the answer was filed back as synthesis).
- `lint` — a health check was run.
- `schema` — the schema was modified.
- `shard` — an index was sharded.

---

## [2026-07-09] ingest | README.md → readme (source) + 5 concepts + 1 gotcha
   Established hub pages from the project front page: three-front-doors, skill-layer,
   layered-architecture, scientific-rigor, validation-and-testing, research-grade-limitations.
   graph: +7 nodes, +4 typed edges

## [2026-07-09] ingest | ALGORITHMS_CHECKLIST_V2.md → algorithms-checklist-v2 (source) + 2 concepts
   Created test-unit-registry, definition-of-done. Updated validation-and-testing and
   research-grade-limitations (campaign-added pending re-validation, unverified complexity claims).
   graph: +3 nodes, +2 typed edges
