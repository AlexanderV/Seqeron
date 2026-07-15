# Vendored: llm-wiki-plugin

- Upstream: https://github.com/praneybehl/llm-wiki-plugin
- Upstream commit: 2d52db91187dc41e2370d064e6ac397c9fef1c71 (2026-05-29)
- License: MIT (see ./LICENSE), Copyright (c) 2026 Praney Behl
- Scope: this directory AND the slash commands in `.claude/commands/wiki/`
  are derived from the same upstream and covered by the same license.
- Pattern credit: Andrej Karpathy's "LLM Wiki" gist (April 2026) — concept only,
  no text copied.

## Local modifications (not upstream)

1. `assets/SCHEMA.md.template` — adapted for docs-as-source: sources are the
   repo's own `docs/**` (no `raw/`), library page types
   (module/api/concept/decision/gotcha/synthesis), `sources:`+`source_commit:`
   frontmatter, git-based staleness rule.
2. `scripts/init_wiki.py` — added `--no-raw` flag.
3. `scripts/wiki_stale.py` — NEW: flags pages whose sources changed after
   `source_commit` or no longer exist.
4. `scripts/wiki_coverage.py` — NEW: reports `docs/**` files not referenced
   by any wiki page.
5. `assets/hooks/pre-commit`, `assets/wiki-ci.yml.sample` — NEW: deterministic
   git/CI gates.
6. `.claude/commands/wiki/*.md` — script paths rewritten from plugin-relative
   to repo-relative; init/lint/ingest extended for the rules above.

To update from upstream: diff this directory against the upstream commit above,
re-apply modifications 1–6 (or cherry-pick), bump the commit hash here.
