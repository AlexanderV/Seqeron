# Project and governance sources

Source pages for project-level documentation, governance, validation methodology, MCP, and skills.

Return to the [[index|Wiki Index]].

- [[readme]] — the project's front-page `README.md`: what Seqeron is, its three entry points, and headline facts.
- [[algorithms-checklist-v2]] — the test-unit validation registry: 364 units (255 done / 109 proposed), DoD, evidence, coverage.
- [[canonical-algorithm-map]] — the `docs/algorithms/` canonical-identity map: alias→canonical test-unit IDs (SEQ-COMPOSITION→SEQ-STATS, SEQ-TM→SEQ-THERMO, GENOMIC-TANDEM→REP-TANDEM), folder-bucket normalization (MolTools/Population_Genetics/K-mer/RnaStructure, four merges complete), one-canonical-doc rules, and the kept legacy/baseline methods (UPGMA, JC/K2P, chi-sq HWE, Nussinov, OLC). Complements [[backlog]] (identity vs coverage).
- [[advanced-testing-checklist]] — 2026-03-19 effectiveness/gap analysis rating the ten testing techniques by applicability, coverage, effort, and P0–P3 priority.
- [[mcp-methods-audit]] — a 2026-01-23 census of the C# API surface: 277 public static methods across 54 classes, with counting rules and the tool-count caveat.
- [[mcp-checklist]] — the SUPERSEDED (2026-07-01) MCP build tracker: a 12-server/241-tool plan never built; DoD gates, error-code catalog, two-tests-per-tool. Live status now in `docs/mcp/MCP_STATUS.md`.
- [[mcp-plan]] — the SUPERSEDED (2026-07-01) MCP design doc (sibling of the checklist): v3→v4 rationale, full 241-tool inventory across 12 servers, and the still-valid standards (naming, JSON-Schema 2020-12, error ranges, SemVer, 2-tests-per-tool, doc contract).
- [[mcp-prompt]] — the CURRENT one-tool-per-session MCP-completion subagent prompt: names the real 11 servers (server→project→file), the 3-part tool DoD (binding+tests+docs), and reverses the old no-business-asserts test policy. Ground truth: `docs/mcp/MCP_STATUS.md`.
- [[mcp-readme]] — the MCP front-door guide (`docs/mcp/README.md`): why MCP (real computation, structured I/O, local/auditable), the current 427-tools/11-servers table, stdio/Codex connection, and two worked tool-only workflows; supersedes the mcp-plan/mcp-checklist tool counts.
- [[skills-strategy]] — the Claude Code skills plan of record (`docs/skills/STRATEGY.md`): many independent dual-mode skills over the source of truth, point-don't-duplicate, the generated-catalog anti-drift mechanism, and the (complete) phase plan.
- [[golden-skills-regression]] — the golden skills regression set (`docs/skills/golden/README.md`): 12 hard dual-mode tasks with expected routing/pipeline/rigor/shape as a manual eval aid; the Phase-3 acceptance instrument for the skills strategy.
- [[sonar-gate-plan]] — the Sonar-ratchet tracker: 66/66 SonarAnalyzer rules moved from report-only to blocking (or silenced-with-justification); Groups A–E, notable behaviour changes, and review-not-blind-fix on the S1244/S125 giants.
- [[findings-register]] — the validation-campaign disposition ledger: every note across all 86 per-unit reports triaged into fixed-now / feasible / not-possible / by-design; green-washing detection; a 2026-06-12 snapshot superseded by the 2026-06-24 re-validation reset.
- [[validation-ledger]] — the live per-unit validation status board (VALIDATION_LEDGER.md): Stage A / Stage B pass-fail matrix across three phases (86 implemented + 24 new campaign units, 148 Phase-2, 12 enhanced), refreshed by reset banners after code churn; ground truth for *where things stand*, superseding the findings-register snapshot.
- [[validation-protocol]] — the validation *methodology* (VALIDATION_PROTOCOL.md): one fresh session per unit in a context separate from the implementer's, Stage A (validate the description vs external primary sources) before Stage B (validate the code), ending in exactly one of ✅ CLEAN / 🔧 LIMITED; the per-session prompt, checklists, verdict legend, and report template the ledger operationalizes.
- [[limitations]] — the validated operating-envelope document (LIMITATIONS.md): the human-readable catalog of what the library does NOT do, every row BY-DESIGN + `✅ CLEAN`; three kinds (irreducible / data-blocked / scope) across ~13 units; the research-vs-clinical disclaimer.
