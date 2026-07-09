## LLM Wiki

This project maintains an LLM-curated wiki at `wiki/` following Andrej Karpathy's "LLM Wiki" pattern (https://gist.github.com/karpathy/442a6bf555914893e9891c11519de94f).

Before answering questions that rely on knowledge accumulated in this project, read `wiki/index.md` (or the relevant shard under `wiki/indexes/` if the wiki has been sharded) and use its one-line summaries to find the pages you need. Cite with `[[wikilinks]]`. If the index does not surface good candidates, fall back to `wiki_search.py` from the `llm-wiki` skill for BM25-ranked retrieval.

Sources are the project's own documentation — everything under `docs/**` **plus** the top-level markdown files in the repo root (`README.md`, `ALGORITHMS_CHECKLIST_V2.md`, `ALGORITHMS_ROADMAP.md`, and similar root `*.md`). There is **no `raw/` directory**, and the wiki never edits these files. To capture a source, follow the `llm-wiki` skill's ingest workflow (`/wiki:ingest <repo-relative-path>`, e.g. `/wiki:ingest docs/api/parsers.md` or `/wiki:ingest README.md`): decide placement by page type (`source`, `module`, `api`, `concept`, `decision`, `gotcha`, `synthesis`); reference docs by repo-relative path and record `sources:` + `source_commit:` in frontmatter rather than copying prose; identify touched pages and make surgical `str_replace` updates rather than rewrites; update the index; append a one-line entry to `wiki/log.md`.

A page is **stale** when any path in its `sources:` has commits after its `source_commit` (`git log <source_commit>..HEAD -- <path>`); `/wiki:lint` flags these.

Scaling discipline: atomic pages (400-line soft cap, 800-line hard cap), sharded indexes past ~150 pages or 300 index lines, required YAML frontmatter on every page, `[[wikilinks]]` for every cross-reference.

Full conventions live in `wiki/SCHEMA.md`. Treat it as authoritative when it disagrees with this summary.
