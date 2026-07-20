---
name: wiki-ingest-doc
description: |
  Ingest ONE newly-added or changed `docs/` file into the project LLM Wiki at `wiki/`, routing by ASPECT — each doc group is ingested differently. Use this whenever a single source doc is added/changed and its knowledge must land in the wiki, whatever its aspect: an algorithm spec (`docs/algorithms/**`), a validation artifact (`docs/Evidence/**`), a validation verdict (`docs/Validation/**`), a testing methodology (`docs/checklists/**`), an MCP-layer doc (`docs/mcp/**`), or a top-level/refactoring/skills/templates doc. Adding one algorithm typically drops SEVERAL such files (spec + Evidence + report + MCP tool) — run this per file; each is routed by its own group rules, not bundled blindly. Trigger phrases: "ingest this doc into the wiki", "I added <file>, update the wiki", "process the new Evidence/report/spec/tool", "reflect this validation report in the wiki", "wire the new MCP tool into the wiki". Do NOT use this for the multi-hundred-file backfill — that is the `/wiki:ingest` loop over `tools/wiki-ingest/next_pending.py`. This skill is the atomic, aspect-routed counterpart for the steady state: one doc = one economical, correctly-routed ingest.
---

# Ingest one added doc into the LLM Wiki (aspect-routed)

Steady-state maintenance skill. When a single `docs/` file is added or changed, this folds it into `wiki/` **using the routing rules for its aspect** — different groups ingest differently. It is the atomic counterpart to the bulk backfill campaign (`/wiki:ingest` + `tools/wiki-ingest/next_pending.py`); use the campaign to drain a large `- [ ]` backlog, use **this** when one (or a few) docs arrive.

Adding an algorithm usually creates several files of different aspects (spec, Evidence, validation report, MCP tool). Do **not** collapse them into one blind pass — run this once per file, and let each route by its group. The files still converge on the same concept via reuse/cross-linking; the *treatment* differs.

## Input

The path(s) of the added/changed doc(s). If the user names an algorithm or unit-id instead of paths, glob `docs/` for the related files first, then process each one through the routing below.

## Aspect routing — pick {CONTEXT} by the file's group

- **`docs/algorithms/**`, `docs/Evidence/**`** → per-algorithm. Hub `[[algorithm-validation-evidence]]` is appropriate; a dedicated concept ONLY if the algorithm is genuinely new, else REUSE the existing `wiki/concepts/` page. For `docs/algorithms/**` check `wiki/backlog.md`: a done-via-concept slug is already covered (skip); a pending one is ingested as **enrichment of the matching concept with the canonical primary spec** (algorithm doc = primary source; Evidence = validation artifact).
- **`docs/Validation/**`** → governance; never force the `[[algorithm-validation-evidence]]` hub. This group has THREE sub-shapes — route by which one the file is:
  - **Governance doc** (`LIMITATIONS.md`, `FINDINGS_REGISTER.md`, `VALIDATION_LEDGER.md`, `VALIDATION_PROTOCOL.md`, or similar policy/methodology text) → create or ENRICH a durable governance concept, not just a source page. Canonical homes: `LIMITATIONS.md` → `[[operating-envelope-and-limitation-policy]]`; `FINDINGS_REGISTER.md` → `[[validation-findings-disposition]]`; `VALIDATION_LEDGER.md` / `VALIDATION_PROTOCOL.md` → a source-summary tied to `[[validation-and-testing]]` + `[[test-unit-registry]]` (+ `[[validation-ledger]]` / `[[validation-protocol]]`). Wire into the governance cluster (`[[validation-and-testing]]`, `[[build-quality-gate]]`, `[[scientific-rigor]]`). If the canonical concept already exists, ENRICH it — do not fork.
  - **Per-unit report** (`docs/Validation/reports/<UNIT-ID>.md`) → classify FIRST with the bundled `scripts/report_verdict.py` (see "Bundled scripts"), then split:
    - **CLEAN** (no code defect — plain PASS, or PASS-WITH-NOTES that only closed test gaps): the cheap path, NO subagent and NO standalone page. Run `report_verdict.py --wiki wiki --reports {FILE} --apply` — it appends a one-line verdict row to the registry `wiki/sources/validation-verdicts.md`. The main agent then just marks done + commits. (~94% of reports are CLEAN.)
    - **DEFECT** (a real code defect was found/fixed — the concept may describe the old behavior): the expensive path — spawn the full subagent, which writes `wiki/sources/<unit>-report.md`, corrects the concept, and adds a gotcha (step 5). (~6% of reports.)
  - **Generated registry** (`MCP_STATUS.md`, `traceability.md`, status/matrix ledgers) → coverage-excluded: mark done, one log line, NO page.
- **`docs/checklists/**`** → testing methodology. Link to `[[validation-and-testing]]` and `[[advanced-testing-checklist]]`.
- **`docs/mcp/**`** → MCP layer. Link to `[[three-front-doors]]`, `[[skill-layer]]`, `[[mcp-plan]]`, `[[mcp-checklist]]`, `[[mcp-prompt]]`. **A single `docs/mcp/tools/**` tool doc is NOT its own page** — add the tool name to the `mcp_tools:` frontmatter of the concept it wraps, and one line to `wiki/concepts/mcp-tool-catalog.md`. Do the frontmatter append with the **bundled script** (see "Bundled scripts" below), not by hand:
  ```
  python .claude/skills/wiki-ingest-doc/scripts/mcp_map.py --wiki wiki --tools {FILE}          # dry-run: shows confirmed / proposed / unmatched
  python .claude/skills/wiki-ingest-doc/scripts/mcp_map.py --wiki wiki --tools {FILE} --apply --trust-proposed   # write, AFTER eyeballing the proposed concept
  ```
  A `proposed` match is lexical and can collide across domains (e.g. `normalize_variant` → *tumor-normal*), so glance at it before `--trust-proposed`; `confirmed` (already present) is zero-risk; `unmatched` means no concept yet. **Fallback — no catalog yet:** create it minimally (title + the tool's server section + the MCP-hub cross-links above), seeded with this tool. If the tool wraps an algorithm with no concept, record it under the catalog's "Unmapped / other" section (create a concept only if the algorithm is genuinely new). `README.md` → an overview source page; `MCP_STATUS.md` / `traceability.md` → generated registries, excluded (mark done, no page).
- **everything else** (top-level, `docs/refactoring/**`, `docs/skills/**`, `docs/templates/**`, `docs/external/**`) → judge by content; no forced hub.

## Procedure

**Cheap fast-paths first — skip the subagent when a script can do it deterministically:**
- **Validation report** → run `scripts/report_verdict.py --wiki wiki --reports {FILE} --apply`. If it classifies the report **clean**, the registry row is already written — go straight to the checklist + commit below, NO subagent. Only if it reports **defect** do you continue to the subagent path.
- **Single MCP tool doc** → run `scripts/mcp_map.py` (see the `docs/mcp/**` routing) to append `mcp_tools:`; the subagent, if any, only writes the catalog line.
- **Generated registry** (`MCP_STATUS.md`, `traceability.md`) → mark done, one log line, no subagent, no page.

**Verify after every write** (cheap paths and subagent alike): once anything under `wiki/` changed, run `python .claude/skills/llm-wiki/scripts/wiki_graph_lint.py wiki/`. If it is not clean, revert the last edit and fix — a malformed `mcp_tools:` block or registry row must never be committed. (The subagent already does this as step 6; for the script fast-paths the main agent does it.)

For every other aspect (and for **defect** reports), delegate the read/write to ONE `general-purpose` subagent (keeps main context clean), then the main agent does the checklist + commit. Spawn it with this task (substitute `{FILE}` and the `{CONTEXT}` chosen above):

```
You are ingesting ONE source file into the project's LLM Wiki at `wiki/` (relative to the repo root; run from there). Work autonomously — no questions. Only create/edit files inside `wiki/`; NEVER edit files under `docs/`.

Source file: {FILE}
Context (aspect routing): {CONTEXT}

Be economical: prefer REUSING/cross-linking existing pages over new ones. Create a dedicated concept page ONLY if the topic is genuinely distinct, wiki-worthy, and not yet represented.

1. Read `wiki/SCHEMA.md` first.
2. Read the source (chunk-read if large).
3. Survey the wiki (`wiki/index.md`, `wiki/sources/`, `wiki/concepts/`, `wiki/gotchas/`, `wiki/backlog.md`); read the candidate pages this source touches. Before creating ANY new concept, run the duplicate guard: `python .claude/skills/wiki-ingest-doc/scripts/find_concept.py --wiki wiki --method <Class.Method>` (or `--name "<title>"`). If it says EXISTS → ENRICH that page, never fork a second one.
4. Apply the aspect routing above:
   - Write a concise source-summary in `wiki/sources/` ONLY when the aspect calls for one (NOT for a CLEAN validation report, NOT for a single MCP tool doc). Surgically update touched pages via str_replace. New concept/entity pages ONLY if warranted (each with >=1 inbound link) AND only after `find_concept.py` reported NONE.
   - For a CLEAN report: fold a one-line verdict into the concept + add the report to `sources:`; no page.
   - For a governance doc (LIMITATIONS/FINDINGS_REGISTER/VALIDATION_LEDGER/VALIDATION_PROTOCOL): create or ENRICH its canonical governance concept per the routing; do not just make a bare source page.
   - For an MCP tool doc: add the tool to the concept's `mcp_tools:` frontmatter + a catalog line; if `wiki/concepts/mcp-tool-catalog.md` does not exist, create it minimally (title + this tool's server section + MCP-hub links) seeded with this tool; no per-tool page.
5. GOTCHA (do NOT skip — this is why `wiki/gotchas/` was nearly empty). **First run the bundled extractor** on an algorithm/spec doc — it resolves the concept + its `mcp_tools:`, surfaces the trap-bearing sections, and tells you whether a gotcha already covers those tools:
   ```
   python .claude/skills/wiki-ingest-doc/scripts/gotcha_candidate.py --wiki wiki --doc {FILE}
   ```
   - **VERDICT `COVERED`** → an existing gotcha already binds these tools; do nothing (or enrich it).
   - **VERDICT `NO_SIGNAL`** → no deterministic sharp edge; still eyeball the trigger checklist, then record **"no gotcha"** in your report — don't silently omit.
   - **VERDICT `LIKELY_GOTCHA`** → read the `*`-marked sections it surfaced and write the gotcha with the **full wiring** below (all four steps — a gotcha without `mcp_tools:` will NOT light up the `(!)` marker in `seqeron-discovery`, which is the whole point).

   **Trigger checklist** (write one gotcha per distinct trap that applies):
   - **Heuristic / profile / *Simplified*, not the trained/calibrated/full method the name implies** — the dominant trap here: a composition/propensity/rule heuristic (disorder, coiled-coil, Chou-Fasman, chromatin-state) or a *Simplified*-status subset (CNV depth-only, VEP-without-reference-window, PROSITE-not-Pfam), sold by a name that sounds like the trained/complete tool. Say what it is NOT and what to use for decision-grade work.
   - **Inverted / counter-intuitive direction** — high score = *low* quality/complexity (e.g. DUST), lower = better.
   - **Silent edge case** — a value/branch that returns a plausible-but-wrong default instead of erroring (empty → 0, NaN window → neutral CN2, non-integer input → observed richness).
   - **Divergence from a reference tool** — result differs from Jellyfish `-C` / BLAST-RBH / DESeq2 / cutadapt / EMBOSS / Biopython / the cited paper under some input.
   - **Clamped / bounded range** — output pinned to `[min,max]` (ENC 20–61); values outside are silently clipped.
   - **"Not the exact published variant"** — a lower bound or simplification (⌈b/2⌉ vs exact reversal distance; linear gap not affine; greedy-SCS not optimal).
   - **Cross-domain name collision** — same term means different things (SV `genotype` vs popgen `genotype`; tumor-`normal` vs `normalize`).
   - **Unit / coordinate trap** — 0- vs 1-based, ×100 vs fraction, amino-acids vs nucleotides, which Tm formula, inclusive vs exclusive ends.

   **Full wiring for each gotcha you write** (match the ~24 existing pages exactly):
   1. Create `wiki/gotchas/<slug>.md` with frontmatter: `type: gotcha`, `title`, `tags: [<area>, gotcha]`, **`mcp_tools:`** (the concept's tool list from the extractor — REQUIRED for the discovery `(!)` marker), `sources:` (the algorithm doc), `source_commit:` (`git rev-parse HEAD`), `created`/`updated`. Body = **The trap.** / **Why it bites.** / **What to rely on instead.**, ending with `Full model: [[<concept>]].`
   2. Backlink from the concept: append `\n**Sharp edge:** [[<gotcha-slug>]] — <one-liner>.\n` to `wiki/concepts/<concept>.md` (this is what stops the gotcha being an orphan — `indexes/` links do NOT count as inbound).
   3. Add one line to `wiki/indexes/gotchas.md`: `- [[<slug>]] — <one-liner>.`
   4. Bump the count in `wiki/index.md` (`- [Gotchas](indexes/gotchas.md) - N ...`).
6. Update `wiki/index.md`; append ONE line to `wiki/log.md`. If you covered a `wiki/backlog.md` slug, move it pending -> covered. Add typed concept-to-concept edges only if the source supports them (edges on concept pages, never on source pages). Run `python .claude/skills/llm-wiki/scripts/wiki_graph_lint.py wiki/` and `python .claude/skills/llm-wiki/scripts/wiki_graph_extract.py wiki/`; append `   graph: +N nodes, +M typed edges` to the log (harmless cp1252 print may appear on Windows AFTER writing — ignore; do not redirect its output).
7. In every wiki page you create/update set `sources:` to {FILE} and `source_commit:` to `git rev-parse HEAD`. Reference docs by path; do not copy long passages.

Do NOT run git add/commit. Do NOT modify WIKI_INGEST_CHECKLIST.md.
Report: pages created, pages updated, gotcha yes/no, contradictions flagged, follow-ups.
```

After the subagent reports (main agent):

1. Mark the processed file done in the checklist, if it appears there:
   ```
   python tools/wiki-ingest/mark_done.py "{FILE}"
   ```
2. `git add wiki WIKI_INGEST_CHECKLIST.md` (exactly this — not `-A`; leave `.claude/` and `docs/` alone).
3. `git commit -m "docs(wiki): ingest <basename without .md>"` and end the message with:
   ```
   Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>
   ```
4. If several files were added for one algorithm, repeat per file (each routed by its own aspect), then summarize once.

## Guardrails

- **Route first, then ingest.** The aspect decides the treatment — never process an Evidence file, a CLEAN report, and an MCP tool the same way.
- **Reuse over create.** New concept page only when genuinely distinct and unrepresented; enrichment is the default — and `find_concept.py` must report NONE before you create one.
- **Verify after write.** Any change under `wiki/` is followed by `wiki_graph_lint.py`; a non-clean lint is reverted, not committed.
- **Gotcha is not optional.** Start step 5 with `gotcha_candidate.py --doc {FILE}`; on `LIKELY_GOTCHA` write the page with all four wiring steps (`mcp_tools:` frontmatter included, or the discovery `(!)` marker won't fire); on `NO_SIGNAL` report "no gotcha" explicitly.
- **CLEAN reports and single MCP tools get NO standalone page** — a one-line verdict / an `mcp_tools:` frontmatter entry respectively.
- **Never edit files under `docs/`** — the wiki is derived; source docs are read-only.
- **Provenance mandatory:** every touched wiki page carries `sources:` (real doc paths) and `source_commit:` (`git rev-parse HEAD`).
- **Not the campaign.** This is one doc (or the few docs of one change), routed. For a large `- [ ]` backlog, use the `/wiki:ingest` campaign instead.

## Bundled scripts

`scripts/mcp_map.py` — travels with this skill. Deterministic MCP-tool → concept mapper: for each `docs/mcp/tools/**/*.md` it finds the wiki concept that documents the wrapped method and idempotently appends the tool name to that concept's `mcp_tools:` frontmatter. Stdlib only, dry-run by default.

- Modes: `--tools <path…>` (one/few — the per-doc path), `--all` (whole MCP layer), `--pending <CHECKLIST>` (only `- [ ]` tools), `--check` (CI gate), `--catalog PATH` (regenerate the catalog JSON). Writes with `--apply` (confirmed+matched); lexical `proposed` only with `--trust-proposed`. A curated `REJECT` set drops known cross-domain false matches (e.g. `normalize_variant`→tumor-normal) the automatic class-veto can't.
- Confidence tiers in the JSON report: **confirmed** (already in `mcp_tools:` — zero risk), **matched** (Method-ID **provenance bridge** — the method's algorithm doc reaches the concept via the backlog map or the concept's `sources:`, or the Method ID appears literally — trustworthy, written on `--apply`), **proposed** (slug/title lexical match + class-family veto — REVIEW, needs `--trust-proposed`), **ambiguous**, **unmatched**. The provenance bridges are ported from the canonical `tools/wiki-ingest/mcp_map.py`, so this matcher is a superset (provenance precision + token recall + class-family veto + the per-file / `--check` interface). On the full layer today: ~211 confirmed / 6 matched / 8 proposed / 13 ambiguous / 189 unmatched.
- This script **owns the MCP layer end-to-end** — it folded in and retired the old one-shot `tools/wiki-ingest/mcp_map.py` (which was dormant: not called by any live workflow, only ran once to wire the 427-tool layer). It is a true superset: the provenance bridges + the curated `REJECT` overrides + `wrapped_source` capture + `--catalog PATH` (regenerates the per-server `mcp_catalog_data.json` behind `wiki/concepts/mcp-tool-catalog.md`) + the new per-file / `--check` / tiered interface. `batch.py` drives it for a bulk sweep; the live campaign helpers (`next_pending.py`, `mark_done.py`, `build_checklist.py`) remain in `tools/wiki-ingest/`.

`scripts/report_verdict.py` — travels with this skill. Classifies a validation report and registers the cheap ones, implementing "full subagent only when a report found a defect":

- It reads ONLY the Stage-B verdict field (not the whole body, which mentions "not a defect" in prose) and splits: **clean** (plain PASS, or PASS-WITH-NOTES that only closed test-coverage / test-quality / code-echo gaps — the concept is unaffected) vs **defect** (a real code defect was found/fixed — `FAIL`, off-by-one, rounding/probability/fidelity defect, code gap, a named `\`Method\`` — the concept may need correcting).
- With `--apply` it appends the **clean** ones as one-line rows to the registry `wiki/sources/validation-verdicts.md` (created if missing, idempotent) — no page, no subagent. The **defect** ones it lists on stderr for the main agent to hand to the full ingest subagent.
- Modes mirror `mcp_map.py`: `--reports <path…>` / `--all` / `--pending <CHECKLIST>`. On the full set today: **239 clean / 16 defect** (so ~94% skip the subagent entirely). Registry rows carry a `[[concept]]` link resolved via `mcp_map`'s matcher.

`scripts/find_concept.py` — the **duplicate-concept guard**. Before creating any new concept, ask "does one already exist?": `--method <Class.Method>` / `--name "<title>"` / `--slug <slug>`. Exit 0 = EXISTS (enrich it), 2 = NONE (create OK), 1 = ambiguous (pick by hand). Read-only. Prevents the single worst incremental-wiki failure — forked duplicate concepts.

`scripts/gotcha_candidate.py` — the **sharp-edge extractor** behind step 5 (makes "gotcha is not optional" mechanical instead of a manual prose scan). For one algorithm doc it resolves the concept + its `mcp_tools:` (reusing `mcp_map`'s provenance bridges), extracts the doc's trap-bearing sections (Deviations / Assumptions / Scope / Limitations / *Simplified* — `*`-marked when they set the verdict; weak Edge-case tables surfaced but only counting when they hit a trap phrase), lists trap phrases, checks whether an existing gotcha already binds those tools, and prints a verdict: **COVERED** / **LIKELY_GOTCHA** / **NO_SIGNAL**. Read-only, stdlib, UTF-8-safe, `--json`, and a `--check` CI gate (exit 1 when a doc shows uncovered gotcha signal). It surfaces candidates high-recall; the subagent still judges keep-vs-skip and writes the prose + the four wiring steps.

`scripts/batch.py` — one-time backlog driver (`reports` | `mcp`). Runs the mapper over `--all`, marks the safely-done rows in the checklist, and **prints** (never runs) the git commit line. `reports --apply` registers the 239 CLEAN verdicts + marks them done; `mcp --apply` marks the 211 CONFIRMED tools done. The risky tails (proposed / ambiguous / unmatched / defect) are left pending for a subagent.

**Shared script conventions:** stdlib only; **dry-run by default** (`--apply` to write); UTF-8-safe; **fail loud** (exit 2) if run outside the repo root (no concepts/tools found); every mapper has a **`--check`** CI gate (exit 1 if the given tools/reports are not yet reflected in the wiki). `mcp_map` uses the Method-ID **class family** to veto cross-domain false matches (a `StructuralVariantAnalyzer` tool can't map to a popgen `genotype` concept).
