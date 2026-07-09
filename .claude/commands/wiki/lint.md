---
description: Run a structural and semantic health check on the LLM Wiki.
argument-hint: "[--suggest-pages]"
---

Run a lint pass on the wiki using the `llm-wiki` skill's lint workflow.

Arguments (if any): $ARGUMENTS

1. Read `wiki/SCHEMA.md` first if you haven't this session — it may declare additional lint rules.
2. Run the structural lint script: `python .claude/skills/llm-wiki/scripts/wiki_lint.py wiki/`. Add `--suggest-pages` if I asked for it. This catches orphans, broken wikilinks, oversized pages, missing frontmatter, stale pages, and duplicate slugs.
2a. Run the staleness check: `python .claude/skills/llm-wiki/scripts/wiki_stale.py wiki/`. It flags pages whose `sources:` docs changed after `source_commit`. For each stale page, propose a re-ingest of the changed doc.
2a2. Run the coverage check: `python .claude/skills/llm-wiki/scripts/wiki_coverage.py wiki/ --docs-dir docs`. For each uncovered doc, ask me whether to ingest it or exclude it (trivial files can be excluded via --exclude).
2b. If `wiki/graph/ontology.yaml` exists, also run `python .claude/skills/llm-wiki/scripts/wiki_graph_lint.py wiki/`. Triage typed-edge findings (unknown predicates, missing evidence, broken object refs, alias collisions) the same way as structural ones — propose fixes, don't apply silently. After approved fixes, suggest running `python .claude/skills/llm-wiki/scripts/wiki_graph_extract.py wiki/`.
3. Walk through the structural findings with me and propose a fix for each — don't apply fixes silently. Wait for my approval per item or in batch.
4. Run the semantic pass: read the most recently updated pages and the most-linked-to hubs (use `python .claude/skills/llm-wiki/scripts/wiki_search.py "" --top-linked 10` to find them). Look for contradictions with older pages, unstated cross-references, and concepts mentioned in many pages but lacking their own page.
5. Surface gaps: topics I've expressed interest in but haven't ingested yet, claims that newer sources have superseded.
6. After applying approved fixes, update the index entries that changed and append one line to `log.md` summarizing the lint pass.
7. If the lint report is too long for me to triage comfortably, that's a signal the cadence is wrong (lint more often) or the schema needs revision — call it out.
