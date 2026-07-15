---
description: Query the LLM Wiki — answer a question from accumulated knowledge with citations.
argument-hint: "<question>"
---

Answer this question against the wiki using the `llm-wiki` skill's query workflow:

Question: $ARGUMENTS

Follow the full query procedure documented in the skill:

1. Read `wiki/SCHEMA.md` first if you haven't this session.
2. Read `wiki/index.md` (or the relevant `wiki/indexes/<type>.md` shard if the index has been sharded) to identify candidate pages.
3. Read the candidate pages. Follow `[[wikilinks]]` selectively — don't recursively chase every link.
4. If the index doesn't surface good candidates, normalize the question into compact search terms in the wiki corpus's dominant language. Do this for any query whose wording does not align with the corpus language; preserve identifiers, symbols, method names, and quoted strings. Search concepts first with `python .claude/skills/llm-wiki/scripts/wiki_search.py "<normalized query terms>" --type concept --top 10`. If that does not surface a good candidate, retry across all page types with `--dedup-provenance --prefer-type concept`. Use `--tag <tag>` only when it follows from the question.
5. If the question asks "what links to X", use `python .claude/skills/llm-wiki/scripts/wiki_search.py "" --backlinks <slug>` instead of grep.
5b. If the question is relational ("what's connected to X", "who proposed Y", "trace the path from A to B") and `wiki/graph/graph.sqlite` exists, run `python .claude/skills/llm-wiki/scripts/wiki_graph_query.py wiki/ neighbors --node <id>` (or `facts` / `edges` / `path`) to get structured neighbors before reading pages. Use the typed-edge results to choose which wiki pages to open — do not answer from graph rows alone for high-stakes claims.
6. Synthesize the answer with `[[wikilink]]` citations to the wiki pages used. Surface contradictions explicitly rather than picking a side.
7. If the wiki has no relevant content, say so plainly — do not confabulate. Suggest sources I might want to ingest to fill the gap.
8. If the answer represents new connection-making, offer to file it back into `wiki/synthesis/` so future queries benefit. Default to offering; let me decline for trivial answers.
