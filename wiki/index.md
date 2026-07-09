# Wiki Index

The catalog of all pages in this wiki. Each entry: a wikilink to the page and a one-line summary. The LLM reads this first when answering queries to identify candidate pages.

Keep summaries tight — one line each. The index is engineered to be cheap to read; a fat index defeats its purpose.

When this file exceeds ~300 lines or the wiki passes ~150 pages, shard into `wiki/indexes/<type>.md` and replace this file with a directory of shards. See the `scaling-playbook.md` reference in the `llm-wiki` skill for the migration procedure.

---

## Sources

- [[readme]] — the project's front-page `README.md`: what Seqeron is, its three entry points, and headline facts.
- [[algorithms-checklist-v2]] — the test-unit validation registry: 364 units (255 done / 109 proposed), DoD, evidence, coverage.
- [[advanced-testing-checklist]] — 2026-03-19 effectiveness/gap analysis rating the ten testing techniques by applicability, coverage, effort, and P0–P3 priority.
- [[mcp-methods-audit]] — a 2026-01-23 census of the C# API surface: 277 public static methods across 54 classes, with counting rules and the tool-count caveat.
- [[mcp-checklist]] — the SUPERSEDED (2026-07-01) MCP build tracker: a 12-server/241-tool plan never built; DoD gates, error-code catalog, two-tests-per-tool. Live status now in `docs/mcp/MCP_STATUS.md`.
- [[mcp-plan]] — the SUPERSEDED (2026-07-01) MCP design doc (sibling of the checklist): v3→v4 rationale, full 241-tool inventory across 12 servers, and the still-valid standards (naming, JSON-Schema 2020-12, error ranges, SemVer, 2-tests-per-tool, doc contract).
- [[mcp-prompt]] — the CURRENT one-tool-per-session MCP-completion subagent prompt: names the real 11 servers (server→project→file), the 3-part tool DoD (binding+tests+docs), and reverses the old no-business-asserts test policy. Ground truth: `docs/mcp/MCP_STATUS.md`.
- [[sonar-gate-plan]] — the Sonar-ratchet tracker: 66/66 SonarAnalyzer rules moved from report-only to blocking (or silenced-with-justification); Groups A–E, notable behaviour changes, and review-not-blind-fix on the S1244/S125 giants.

## Concepts

- [[three-front-doors]] — one algorithm engine exposed through skills, the C# API, and MCP, with identical results.
- [[skill-layer]] — the Agent-Skill routing + discipline layer that keeps 427 tool schemas out of the model's context.
- [[layered-architecture]] — the strict up-only dependency layering (Levels 0–4) enforced by architecture tests.
- [[scientific-rigor]] — runtime honesty: `LimitationPolicy`, tool-only computation, provenance on every result.
- [[validation-and-testing]] — 22k+ tests across ten methodologies plus the per-unit validation campaign.
- [[build-quality-gate]] — the build-time static-analysis gate: SonarAnalyzer ratcheted to blocking under `TreatWarningsAsErrors`, fix-or-silence, review-not-blind-fix.
- [[test-unit-registry]] — the area-prefixed Test Unit ID scheme and per-unit record behind the validation effort.
- [[definition-of-done]] — the six-criterion acceptance bar (TestSpec, tests, ≥80% coverage, edge cases, CI, evidence) each unit must clear.

## Gotchas

- [[research-grade-limitations]] — beta, not for clinical use; simplified-subset implementations; internal-only validation.

## Synthesis

(populated as query answers are filed back)
