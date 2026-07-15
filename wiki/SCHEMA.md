# Wiki Schema

This file is the configuration for this wiki. It documents the conventions, page types, tag taxonomy, and any workflow customizations. The LLM reads this first when entering the wiki, and its conventions override the defaults documented in the `llm-wiki` skill.

This file is **co-evolved with the user**. When the LLM notices a recurring pattern in your edits or feedback that isn't here, it will propose adding it. When something here stops fitting, prune it.

## Wiki location

- Wiki root: `wiki/`
- Raw sources: **the project's own documentation** — there is NO `raw/` directory. Two locations count as sources:
  1. everything under `docs/**`, and
  2. top-level markdown files in the repo root (e.g. `README.md`, `ALGORITHMS_CHECKLIST_V2.md`, `ALGORITHMS_ROADMAP.md`).
  Sources are already immutable in git. The LLM reads them but NEVER edits them.
- Asset/image storage: `wiki/assets/`

## Source-of-truth rules (project-specific)

1. Sources are the files under `docs/**` **plus** the top-level markdown files in
   the repo root (`README.md`, `ALGORITHMS_CHECKLIST_V2.md`, `ALGORITHMS_ROADMAP.md`,
   and any similar root `*.md`). All are read-only: the LLM reads and cites them
   but NEVER edits them.
2. Wiki pages **reference** sources by repo-relative path — they must not copy large
   passages verbatim. Synthesize, link, and cite.
3. Every page derived from a source records in frontmatter:
   - `sources:` — list of repo-relative paths (e.g. `docs/api/parsers.md`, `README.md`)
   - `source_commit:` — the git commit hash of HEAD used as the synthesis baseline
4. **Staleness rule:** for each path in `sources:`, find its latest commit after
   `source_commit`. The page is stale unless its own Markdown file was refreshed in
   that same commit or a later descendant commit. This makes a source plus its
   ingest page safe to commit atomically without hiding a later source-only change.
   In the index, a staged source change likewise requires its affected page to be
   staged, so the pre-commit check catches drift before it enters history.
   `/wiki:lint` runs `scripts/wiki_stale.py` to apply the rule deterministically.
5. Ingest is triggered per changed source: `/wiki:ingest <repo-relative-path>`
   (e.g. `/wiki:ingest docs/api/parsers.md` or `/wiki:ingest README.md`).

## Coverage exclude policy

`scripts/wiki_coverage.py` treats every markdown file under `docs/**` as an
expected wiki source. Several subtrees are **generated or reference-only** and are
NOT ingestion targets — they would otherwise swamp the coverage signal (they
accounted for ~693 of 1135 uncovered docs at the 2026-07-09 lint). Always pass
these `--exclude` globs when running the coverage check:

```
--exclude 'docs/mcp/tools/**'          # 427 generated per-tool reference docs
--exclude 'docs/Validation/reports/**' # generated per-run validation reports
--exclude 'docs/refactoring/**'        # internal/historical refactoring logs
--exclude 'docs/skills/_generated/**'  # generated skill artifacts
--exclude 'docs/templates/**'          # doc templates
--exclude 'docs/mcp/MCP_STATUS.md'     # live campaign/status ledger, not a knowledge source
--exclude 'docs/mcp/traceability.md'   # generated tool→test traceability matrix
--exclude 'docs/skills/golden/tasks.md' # golden regression task list, reference-only
```

Everything else under `docs/**` remains in scope. `docs/algorithms/**` is
explicitly **in scope**: each algorithm doc is reconciled against its synthesizing
concept page (recorded in `wiki/backlog.md`), not excluded.

## Page types

This wiki uses these page types, each with a dedicated subdirectory:

- `source` (in `wiki/sources/`) — one summary page per ingested documentation file (maps 1:1 to a path under `docs/**`).
- `module` (in `wiki/modules/`) — pages about a module/package of the library.
- `api` (in `wiki/api/`) — pages about the public API surface (types, functions, contracts).
- `concept` (in `wiki/concepts/`) — cross-cutting ideas, patterns, abstractions.
- `decision` (in `wiki/decisions/`) — architectural decisions and their rationale.
- `gotcha` (in `wiki/gotchas/`) — non-obvious behavior, limitations, sharp edges.
- `synthesis` (in `wiki/synthesis/`) — cross-cutting analyses, comparisons, query answers filed back.

Add additional types here as the wiki evolves.

## Tag taxonomy

(Empty initially. Add tags here as you adopt them, with one-line descriptions. Keep this list small and disciplined — a wiki with 200 tags has effectively no tags.)

Example structure:
- `methodology` — pages about research or analytical methods.
- `open-question` — pages or sections that flag unresolved questions.
- `contested` — pages where sources contradict.

## Page sizing

- Soft cap: 400 lines / ~2,000 words. Consider splitting beyond this.
- Hard cap: 800 lines. Must split.

## Frontmatter requirements

Every page must have:
- `type`
- `title`
- `tags`
- `created`
- `updated`

Plus type-specific:
- `source` pages: `doc_path` (repo-relative path under `docs/**`), `source_commit`, `ingested`
- Non-source pages: `sources` (repo-relative `docs/**` paths) and `source_commit`

## Optional graph metadata

Pages may declare typed graph metadata under a top-level `graph:` key. This is the source of truth for the compiled knowledge graph under `wiki/graph/`. Markdown remains canonical; the graph is a regenerable index. Pages without `graph:` still appear as nodes (derived from `type`/`kind`) and still contribute `mentions` edges from body `[[wikilinks]]`.

```yaml
graph:
  node_id: api:needleman-wunsch-align   # optional; default <node_type>:<slug>
  node_type: api                         # optional; default mapped from type via ontology
  canonical: true                        # mark as canonical when multiple slugs alias the same element
  aliases: [NeedlemanWunsch, global-align]
  relationships:
    - predicate: implements
      object: concept:global-alignment
      source: needleman-wunsch-align     # source/derived page slug
      evidence: "Fills the DP matrix with affine gap penalties per docs/api/alignment.md"
      confidence: high               # high | medium | low
      status: current                # current | historical | proposed | disputed | superseded
      # optional:
      # valid_from: 2025-01-15
      # valid_to: 2026-03-01
      # notes: "..."
      # raw_ref: "docs/api/alignment.md#L42"
      # contradicts: edge-id-or-source-slug
      # supersedes: edge-id-or-source-slug
```

Required fields on every relationship: `predicate`, `object`, `source`, `evidence`, `confidence`, `status`. Predicates and the subject/object types they accept are declared in `wiki/graph/ontology.yaml`. Typed semantic edges must be supported by an explicit source — never emit one inferred from training data alone.

## Index structure

Sharded. The top-level `wiki/index.md` is a small directory that routes readers to:

- `wiki/indexes/sources-project.md` — project-level and governance source pages.
- `wiki/indexes/sources-validation-a-m.md` — per-unit evidence artifacts and validation reports with slugs A–M.
- `wiki/indexes/sources-validation-n-z.md` — per-unit evidence artifacts and validation reports with slugs N–Z.
- `wiki/indexes/concepts.md` — concept pages.
- `wiki/indexes/gotchas.md` — gotcha pages.
- `wiki/indexes/meta.md` — coverage and maintenance indexes.
- `wiki/indexes/synthesis.md` — synthesis pages.

Ingests update the relevant shard rather than the top-level directory. The root routing index and
every shard share a ~300-line cap; split any index that exceeds it by a stable sub-category. For
fuzzy discovery at this wiki's scale,
`python .claude/skills/llm-wiki/scripts/wiki_search.py "<query>"` is the sanctioned
fallback after index-first navigation.

### Query retrieval policy

When index summaries do not surface a good candidate, normalize the question into compact
search terms in the dominant language of this wiki's corpus. Apply this based on mismatch
with the corpus language, never by singling out a particular input language. Preserve stable
identifiers, symbols, method names, and quoted strings.

Search `concept` pages first. They are the synthesized retrieval surface for explanatory
questions. If that pass is insufficient, search all page types with provenance deduplication
and prefer a `concept` representative:

```text
python .claude/skills/llm-wiki/scripts/wiki_search.py "<normalized terms>" --type concept --top 10
python .claude/skills/llm-wiki/scripts/wiki_search.py "<normalized terms>" --top 10 --dedup-provenance --prefer-type concept
```

Use `source` pages directly when the question asks for evidence, a validation report, or what
a specific document says. Query the graph only for relational questions. If a filtered pass
does not return a credible candidate, broaden the search rather than forcing an answer.

## Graph layer

The wiki has an optional compiled graph layer under `wiki/graph/`:

- `wiki/graph/ontology.yaml` — declares node types and predicates. **Tracked.** Edit this when you introduce new predicates or domain types.
- `wiki/graph/nodes.jsonl`, `wiki/graph/edges.jsonl` — generated. Track in git only if you want graph diffs in PRs.
- `wiki/graph/graph.sqlite` — generated. Gitignored by default.
- `wiki/graph/graph.graphml` — generated. Track only if you want to diff it.

Generation is reproducible from markdown via `scripts/wiki_graph_extract.py`. The graph can be deleted at any time and rebuilt without losing knowledge — markdown is canonical.

## Workflow customizations

(Empty initially. Document any deviations from the default ingest/query/lint workflows here.)

## User preferences

(Empty initially. As the user expresses style preferences — "always include a 'Why this matters' section on concept pages", "never use bullet lists in summaries", "prefer comparative tables for synthesis pages" — capture them here so they persist across sessions.)

## Lint cadence

- Structural lint: after every 5 ingests.
- Semantic lint: weekly or after every 20 ingests.
- Gap-finding: monthly.
- Graph lint + extract: after every ingest that adds typed `graph.relationships`.

Adjust based on the wiki's growth rate.
