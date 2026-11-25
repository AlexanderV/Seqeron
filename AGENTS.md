## LLM Wiki

This project maintains an LLM-curated wiki at `wiki/` following Andrej Karpathy's "LLM Wiki" pattern (https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f).

Before answering questions that rely on knowledge accumulated in this project, read `wiki/index.md` (or the relevant shard under `wiki/indexes/` if the wiki has been sharded) and use its one-line summaries to find the pages you need. Cite with `[[wikilinks]]`. If the index does not surface good candidates, fall back to `wiki_search.py` from the `llm-wiki` skill for BM25-ranked retrieval.

Sources are the project's own documentation — everything under `docs/**` **plus** the top-level markdown files in the repo root (`README.md`, `ALGORITHMS_CHECKLIST_V2.md`, `ALGORITHMS_ROADMAP.md`, and similar root `*.md`). There is **no `raw/` directory**, and the wiki never edits these files. To capture a source, follow the `llm-wiki` skill's ingest workflow (`/wiki:ingest <repo-relative-path>`, e.g. `/wiki:ingest docs/api/parsers.md` or `/wiki:ingest README.md`): decide placement by page type (`source`, `module`, `api`, `concept`, `decision`, `gotcha`, `synthesis`); reference docs by repo-relative path and record `sources:` + `source_commit:` in frontmatter rather than copying prose; identify touched pages and make surgical `str_replace` updates rather than rewrites; update the index; append a one-line entry to `wiki/log.md`.

A page is **stale** when any path in its `sources:` has commits after its `source_commit` (`git log <source_commit>..HEAD -- <path>`); `/wiki:lint` flags these.

Scaling discipline: atomic pages (400-line soft cap, 800-line hard cap), sharded indexes past ~150 pages or 300 index lines, required YAML frontmatter on every page, `[[wikilinks]]` for every cross-reference.

Full conventions live in `wiki/SCHEMA.md`. Treat it as authoritative when it disagrees with this summary.

## Licences and assemblies

The repository's main assemblies are MIT (`LICENSE`). Ports of GPL / LGPL programs (made by reading their source) live in **separate assemblies named by subject**
(`Seqeron.Genomics.Kalign`, `...PartialOrderAlignment`, `...ClustalW`, `...TCoffee`, ... and the MCP servers `Seqeron.Mcp.MultipleAlignment`, ...), each with its own `LICENSE`,
an SPDX header (`// SPDX-License-Identifier: <id>` + origin note) in EVERY `.cs` file (also of its test project), `<PackageLicenseExpression>` in its csproj, its own test
project(s) and tool docs under `docs/mcp/tools/<server>/`. Rules for agents:

* The MIT assemblies must NEVER reference a copyleft assembly (project reference, `InternalsVisibleTo`, assembly reference, reflection string, MCP tool). The copyleft assemblies may
  reference the MIT ones; what they need from an MIT assembly is made `public` there (documented), never reached through `InternalsVisibleTo`.
* A new port of a GPL / LGPL / AGPL program goes straight into the matching copyleft assembly (new assembly named by meaning if the subject is new - never "Gpl" / "Copyleft"),
  with the header, LICENSE, tests, MCP tools (server of that licence) and an entry in `docs/Validation/review-2026-09/LICENSES.md`. A permissive port stays in the MIT assemblies and
  gets its notice in `THIRD-PARTY-NOTICES.md`. A custom / non-commercial / unknown licence is NOT decided by an agent: it goes into the "NEEDS OWNER DECISION" section of LICENSES.md.
* `tests/Seqeron/Seqeron.Genomics.Tests/Architecture/LicenceSeparationTests.cs` enforces the separation (and the SPDX headers); `python3 scripts/skills/gen-catalog.py --check`,
  `python3 scripts/docs/sync_mcp_tool_docs.py --check` and `python3 scripts/skills/sync-github-skills.py --check` know the extra MCP servers.
* Inventory, black-box oracles and open owner decisions: [docs/Validation/review-2026-09/LICENSES.md](docs/Validation/review-2026-09/LICENSES.md); notices of the permissive programs:
  [THIRD-PARTY-NOTICES.md](THIRD-PARTY-NOTICES.md).
