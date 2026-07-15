# wiki-ingest-doc

**Steady-state, aspect-routed ingest of one `docs/` file into the Seqeron LLM Wiki (`wiki/`).**
One doc in → the right thing happens in the wiki, cheaply, without a full subagent unless the doc actually earns one.

This is the atomic, per-document counterpart to the bulk backfill campaign
(`/wiki:ingest` looping over `tools/wiki-ingest/next_pending.py`). Use the **campaign** to drain a
large `- [ ]` backlog; use **this skill** when one (or a few) docs land in the steady state.

---

## TL;DR

- **Trigger:** "ingest this doc into the wiki", "I added `<file>`, update the wiki", "reflect this validation report", "wire the new MCP tool". It auto-fires via its `SKILL.md` description.
- **It routes by aspect** — an algorithm spec, an Evidence file, a validation report, an MCP tool, a governance doc and a checklist are each ingested *differently*.
- **It is cheap by design** — two bundled scripts do the mechanical 90%; a subagent is spawned only for the ~10% of docs that need real synthesis or a concept correction.
- **It captures sharp edges** — non-obvious behavior becomes a `gotcha`, the highest-signal layer for an assistant.

---

## When it fires — and when it must not

| Fire this skill | Do NOT — use instead |
|---|---|
| One new/changed `docs/` file lands | A multi-hundred-file backlog → `/wiki:ingest` campaign |
| Adding an algorithm (its spec + Evidence + report + tool — run once per file) | Re-organising the whole wiki → `llm-wiki` skill + lint scripts |
| A validation report needs reflecting | Querying the wiki → `/wiki:query` |
| A new MCP tool needs wiring | Auditing health → `wiki_lint` / `wiki_stale` / `wiki_coverage` |

**Adding one algorithm usually drops several files of different aspects.** Run the skill once per file
and let each route by its group — they still converge on the same concept via reuse/cross-linking; only
the *treatment* differs.

---

## Aspect routing (the heart of the skill)

| Aspect (path) | Treatment |
|---|---|
| `docs/algorithms/**`, `docs/Evidence/**` | **Per-algorithm.** Enrich the existing concept with the canonical spec (algorithm doc = primary source, Evidence = validation artifact); new concept only if genuinely distinct. Check `wiki/backlog.md`. Hub `[[algorithm-validation-evidence]]`. |
| `docs/Validation/reports/**` | **Verdict-split** (see `report_verdict.py`). CLEAN → one registry row, **no page, no subagent**. DEFECT → full subagent: report page + concept correction + gotcha. |
| `docs/Validation/{LIMITATIONS,FINDINGS_REGISTER,VALIDATION_LEDGER,VALIDATION_PROTOCOL}.md` | **Governance.** Create/enrich the canonical governance concept (`[[operating-envelope-and-limitation-policy]]`, `[[validation-findings-disposition]]`, …), not a bare source page. Never force the algorithm hub. |
| `docs/mcp/tools/**` | **Not its own page.** `mcp_map.py` appends the tool to the concept's `mcp_tools:` frontmatter + one catalog line. |
| `docs/mcp/README.md` | Overview source page. |
| `docs/mcp/MCP_STATUS.md`, `traceability.md` | **Generated registries** → mark done, one log line, excluded. |
| `docs/checklists/**` | Testing methodology → `[[validation-and-testing]]`, `[[advanced-testing-checklist]]`. |
| everything else (top-level, `refactoring/`, `skills/`, `templates/`, `external/`) | Judge by content; no forced hub. |

---

## Bundled scripts

Both are stdlib-only, **dry-run by default**, UTF-8-safe, and live in `scripts/`. They exist so the
common cases never burn a subagent — the script does the deterministic part, the subagent only takes
the ambiguous tail.

### `scripts/mcp_map.py` — MCP tool → concept

Finds the concept that documents a tool's wrapped method and idempotently appends the tool to that
concept's `mcp_tools:` frontmatter.

```bash
# one new tool (the per-doc path)
python .claude/skills/wiki-ingest-doc/scripts/mcp_map.py --wiki wiki --tools docs/mcp/tools/alignment/global_align.md
python .claude/skills/wiki-ingest-doc/scripts/mcp_map.py --wiki wiki --tools <FILE> --apply --trust-proposed

# the whole layer (one-time batch)
python .claude/skills/wiki-ingest-doc/scripts/mcp_map.py --wiki wiki --all
```

Confidence tiers in the JSON report:

| Tier | Meaning | Risk |
|---|---|---|
| **confirmed** | already in `mcp_tools:` | zero — ground truth |
| **matched** | Method-ID **provenance bridge** (backlog map / concept `sources:` / literal) | trustworthy — written on `--apply` |
| **proposed** | slug/title lexical match + class-family veto | REVIEW — token collisions happen (`normalize_variant` → *tumor-normal*); needs `--trust-proposed` |
| **ambiguous** | several candidate concepts | subagent decides |
| **unmatched** | no concept (a helper, or a gap) | subagent maps or records in the catalog |

On the full 427 tools today: **≈211 confirmed / 6 matched / 7 proposed / 13 ambiguous / 190 unmatched.**
This script **owns the MCP layer end-to-end** — it folded in and retired the old one-shot
`tools/wiki-ingest/mcp_map.py` (dormant: not called by any live workflow). It carries the provenance
bridges + a curated `REJECT` set (drops cross-domain false matches the class-veto misses, e.g.
`normalize_variant`) + `wrapped_source` + `--catalog PATH` (regenerates the JSON behind
`wiki/concepts/mcp-tool-catalog.md`) + the per-file / `--check` / tiered interface. `--apply` writes
`confirmed`+`matched`; lexical `proposed` only with `--trust-proposed` after review.

### `scripts/report_verdict.py` — validation report → verdict registry

Classifies a report by its **Stage-B verdict field only** (not the prose, which says "not a defect"),
and registers the cheap ones.

```bash
python .claude/skills/wiki-ingest-doc/scripts/report_verdict.py --wiki wiki --reports docs/Validation/reports/ALIGN-SEMI-001.md
python .claude/skills/wiki-ingest-doc/scripts/report_verdict.py --wiki wiki --all --apply
```

| Class | Signal | Action |
|---|---|---|
| **clean** | plain PASS, or PASS-WITH-NOTES that only closed test-coverage / test-quality / code-echo gaps | one row in `wiki/sources/validation-verdicts.md` — **no page, no subagent** |
| **defect** | `FAIL`, off-by-one, rounding / probability / fidelity defect, code gap, a named `` `Method` `` | spawn the full subagent (report page + concept correction + gotcha) |

On the full 255 reports today: **239 clean / 16 defect** — ~94 % skip the subagent entirely. That is the
whole point of the recommendation "full subagent only when a report found a defect". Registry rows carry
a `[[concept]]` link (resolved via `mcp_map`'s matcher).

### `scripts/find_concept.py` — duplicate-concept guard

Before a subagent creates a new concept, ask "does one already exist?" — the single worst incremental-wiki
failure is a forked duplicate. Reuses `mcp_map`'s matcher.

```bash
python .claude/skills/wiki-ingest-doc/scripts/find_concept.py --wiki wiki --method SequenceAligner.GlobalAlign
python .claude/skills/wiki-ingest-doc/scripts/find_concept.py --wiki wiki --name "Codon Adaptation Index"
```

Read-only. Exit **0** = EXISTS (enrich it), **2** = NONE (create is OK), **1** = ambiguous (pick by hand).

### `scripts/batch.py` — one-time backlog driver

Ties a mapper to the checklist for the two big sweeps, and **prints** (never runs) the git commit line.

```bash
python .claude/skills/wiki-ingest-doc/scripts/batch.py reports --wiki wiki --apply   # register 239 CLEAN verdicts, mark done
python .claude/skills/wiki-ingest-doc/scripts/batch.py mcp     --wiki wiki --apply   # mark 211 CONFIRMED tools done
```

The risky tails (proposed / ambiguous / unmatched / defect) are left pending for a subagent.

### Shared conventions

All scripts are stdlib-only, **dry-run by default** (`--apply` to write), UTF-8-safe, and **fail loud**
(exit 2) if run outside the repo root. Both mappers expose a **`--check`** CI gate — exit 1 if the given
tools/reports are not yet reflected in the wiki (wire it into CI to fail when a new `docs/` file was added
but never ingested).

---

## Cost model

```
new doc ─► aspect router
             │
             ├─ MCP tool         ─► mcp_map.py        (no subagent)
             ├─ CLEAN report     ─► report_verdict.py (no subagent)  ← ~94% of reports
             ├─ generated ledger ─► mark done         (no subagent)
             │
             └─ everything else  ─► ONE general-purpose subagent ─► main agent commits
                (algorithm spec, Evidence, DEFECT report, governance doc, misc)
```

The scripts turn "427 tools + 255 reports" from *~680 subagents* into *2 script runs + ~16 defect
subagents + ~200 catalog/tail decisions in a single subagent*. That is the ≈55–60 M-token saving the
value audit identified.

---

## What it guarantees — and what it does NOT (read this)

**Guarantees (going forward):** every new doc is routed to the correct treatment; concepts are reused
not duplicated; CLEAN reports and single MCP tools never spawn a page or a subagent; every touched page
carries `sources:` + `source_commit:`; sharp edges become gotchas.

**Does NOT, by itself:**

- **Drain the existing backlog.** 856 pending are the campaign's job. This skill only keeps *new* docs correct. A "complete" wiki also needs: the campaign for the ~210 algorithm docs, `mcp_map.py --all` + `report_verdict.py --all` as one-time batches, and a one-time gotcha backfill of already-found defects.
- **Guarantee synthesis quality.** The concept pages are as good as the subagent's reading; the skill guides but does not grade.
- **Resolve the fuzzy tails.** ~6–40 % of `proposed` mappings need a human/subagent glance — hence `--trust-proposed` and the ambiguous/unmatched buckets.
- **Enforce global consistency.** Duplicate concepts, cross-page contradictions and staleness are caught by periodic `wiki_lint.py` / `wiki_stale.py` / `wiki_coverage.py` (from the `llm-wiki` skill), which this skill references but does not run.

> **A "perfect" wiki = this skill (steady state) + the campaign (backlog) + the two `--all` batches + a gotcha backfill + periodic lint/audit.** This skill is the central, highest-value piece — not the whole recipe.

---

## End-to-end examples

**1. New algorithm spec** (`docs/algorithms/Statistics/New_Thing.md`) → subagent enriches or creates the concept from the canonical spec; marks the backlog slug covered; commit.

**2. New Evidence** (`docs/Evidence/NEW-UNIT-001-Evidence.md`) → subagent writes `wiki/sources/new-unit-001-evidence.md` and enriches the concept; commit.

**3. CLEAN report** (`docs/Validation/reports/NEW-UNIT-001.md`, Stage B PASS) →
`report_verdict.py --reports … --apply` adds one registry row; main agent marks done + commits. **No subagent.**

**4. DEFECT report** (Stage B `FAIL → FIXED`) → `report_verdict.py` flags it; subagent writes the report page, corrects the concept if it described the old behavior, adds a gotcha; commit.

**5. New MCP tool** (`docs/mcp/tools/<server>/<tool>.md`) → `mcp_map.py --tools … ` shows the match;
`--apply --trust-proposed` writes the `mcp_tools:` entry; add one catalog line; commit. **No page.**

---

## Files

```
.claude/skills/wiki-ingest-doc/
├── SKILL.md                 # the skill itself (routing + procedure + guardrails)
├── README.md                # this file
└── scripts/
    ├── mcp_map.py           # MCP tool  → concept mcp_tools: frontmatter (class-family veto, --check)
    ├── report_verdict.py    # report    → CLEAN registry / DEFECT subagent split (--check)
    ├── find_concept.py      # guard     → does a concept already exist? (anti-duplicate)
    └── batch.py             # driver    → one-time reports/mcp sweep + mark-done + commit hint
```

## Extending

- New doc aspect? Add a row to the routing table in `SKILL.md` and, if it's mechanical and repetitive, a matching script in `scripts/`.
- Both scripts share the same shape (`--tools/--reports` · `--all` · `--pending`, dry-run default, JSON report + stderr summary) — copy one as a template.
- Keep the scripts **precision-first**: when a match is uncertain, hand it to the subagent rather than guessing. A wrong auto-edit is worse than an unmatched row.
