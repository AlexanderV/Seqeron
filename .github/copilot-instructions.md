# Seqeron — repo conventions for skill authoring

Shared brief for the Seqeron skill-development pipeline (`.github/agents/`). Every pipeline
agent reads this first. Source of the overall design: [`docs/skills/STRATEGY.md`](../docs/skills/STRATEGY.md).

## What this repo is

`Seqeron.Genomics` (+ `SuffixTree`) — a C#/.NET 10 bioinformatics library: 258 algorithm units,
exposed both as a C# API and as **11 MCP servers = 427 tools**. Heavily documented and validated
(10 test methodologies, an independent validation campaign, `docs/Validation/LIMITATIONS.md`, and a
runtime `LimitationPolicy`). Pre-1.0 / alpha — **not for clinical use**.

## Skills: purpose and catalog

Skills are a thin **routing + discipline** layer that helps solve biological tasks with the library —
they do **not** re-document algorithms. Two runtimes consume them; installs are **byte-identical**
in both trees:

- `.claude/skills/<name>/` — Claude Code
- `.github/skills/<name>/` — GitHub Copilot / VS Code

Current skills: cross-cutting `bio-rigor`, `seqeron-discovery`, `seqeron-dev`; domain
`bio-qc`, `bio-alignment`, `bio-assembly`, `bio-annotation`, `bio-moldesign`, `bio-phylo-popgen`,
`bio-metagenomics`, `bio-chromosome`.

## Non-negotiable authoring conventions

1. **Point, don't duplicate.** Link to `docs/mcp/tools/<server>/<tool>.md` (schemas + `Method ID`),
   `docs/algorithms/<Area>/*.md` (contracts), `docs/Validation/LIMITATIONS.md`. Never copy their content.
   The only generated tool index is `_generated/tools.md` (built by `scripts/skills/gen-catalog.py`) —
   never hand-edit it.
2. **Dual-mode.** Every recipe reads correctly for BOTH an MCP tool-caller and a C# API caller
   (real tool name **and** real `Method ID`). Verify every tool name / Method ID against the docs —
   never guess. Use `python3 scripts/skills/find-tool.py <kw>` to discover/verify.
3. **Rigor by default.** Delegate scientific discipline to `bio-rigor` (tool-only computation,
   provenance, envelope, cross-checking, units/0-based coords) — do not restate its rules.
4. **Envelope-aware.** If a task hits a guarded unit below its `MinimumMode`, carry the STOP rule
   (report + named alternative), never force output. Guarded units + modes live in `LIMITATIONS.md`.
5. **Progressive disclosure.** `SKILL.md` ≤ ~150 lines (trigger + decision guide + a few pipeline
   skeletons); detail goes in `reference/*.md`.
6. **Provenance.** Every recipe ends with which tools / `Method ID`s + parameters were used.
7. **Relative-link depth (a real past bug):** from `SKILL.md` → repo root is `../../../`; from
   `reference/*.md` → repo root is `../../../../`; to a sibling skill from `reference/` is
   `../../<skill>/`. Every link must resolve (`python3 scripts/skills/check-links.py`).

## SKILL.md frontmatter

YAML `name:` (matches dir) and `description:` — the `description` is the **trigger**: concrete
biology verbs/nouns and example prompts so the skill activates automatically. It must read correctly
for both runtimes (no Claude-only or Copilot-only assumptions).

## The pipeline working area

All intermediate artifacts live under `.skill-dev/` (git-ignored working area), never elsewhere:

```
.skill-dev/
  skill-candidates.json                  # Cartographer
  skill-roadmap.md                       # Orchestrator triage
  <candidate_id>/
    escalation.md                        # (optional) Researcher halt
    research-pack/{conventions.md,citations.md,...}   # Domain Researcher
    draft/{SKILL.md, references/*}        # Skill Author
    trigger-attacks.json                 # Trigger Adversary
    compat-report.md                     # Dual-Target Verifier
    evals/evals.json                     # Eval Architect
    eval-runs/auto/results.json          # Eval Runner
  coherence-report.md                    # Coherence Auditor
  PIPELINE-REPORT.md                     # Orchestrator final report
```

Installed skills (`.claude/skills/<id>/`, `.github/skills/<id>/`) are written **only** by the
orchestrator's install step, and must be byte-identical across the two trees.

## Ground-truth references

- Tool docs + `Method ID`: `docs/mcp/tools/**` · per-server list: `src/**/Seqeron.Mcp.<Server>/README.md`
- Algorithm contracts: `docs/algorithms/**`
- Operating envelope / guarded units: `docs/Validation/LIMITATIONS.md`
- Discovery + catalog: `scripts/skills/find-tool.py`, `docs/skills/_generated/catalog.json`
- Strategy + acceptance criteria: `docs/skills/STRATEGY.md` (§7 = a skill is "done" only when it meets them)
- Golden regression tasks: `docs/skills/golden/tasks.md`
