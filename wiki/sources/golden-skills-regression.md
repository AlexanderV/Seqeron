---
type: source
title: "Golden skills regression set (docs/skills/golden/README.md)"
tags: [skills, testing]
doc_path: docs/skills/golden/README.md
sources:
  - docs/skills/golden/README.md
source_commit: b3f950caf701615bb8a0296df6c5368d26dde7ec
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Golden skills regression set

A curated set of hard, realistic biological tasks used as **manual / eval regression
aids** for the skill catalog — the Phase-3 acceptance instrument from [[skills-strategy]]
§7. Each task is paired with expected routing (which skills should fire), the expected
dual-mode pipeline (MCP tool calls + C# `Method ID`s), the rigor checkpoints `bio-rigor`
must enforce, and an expected-shape output. The task list itself is in
`docs/skills/golden/tasks.md` (coverage-excluded, reference-only).

## What it is and is not

It **is** a frozen reference of expected behaviour for a representative slice of the
catalog. It **is not** an automated test: the skills are markdown routers, not code, so
nothing executes — a "pass" is a reviewer's judgement. It does not duplicate schemas or
formulas; those stay in the source of truth, and it cites only *verified* tool names +
Method IDs.

## Reviewer workflow and pass criteria

Paste a task verbatim into a Claude Code session with the skills installed, then check
six things: **routing** (expected skills trigger, no clearly-wrong skill dominates),
**pipeline** (tool chain matches expected tools + order; a different tool for a graded
step is a finding), **rigor** (tool-only computation, provenance block, units +
coordinate base stated, listed cross-checks done), **envelope** (guarded-unit tasks
STOP and report the limitation + its `MinimumMode` instead of forcing a number),
**disclaimer** (decision-relevant results carry the beta / not-for-clinical caveat), and
**shape** (output structure matches). Exact numbers are not graded — the point is
routing + orchestration + rigor.

## Coverage

12 tasks spanning sequence QC, pairwise + MSA alignment, variant calling +
classification, PCR primer design, CRISPR design, phylogeny, popgen (Fst / HWE /
Tajima's D), metagenomics diversity, assembly / N50, and chromosome centromere +
GC-skew/origin — including 4 cross-domain chains and 3 guarded-unit tasks (→
`PARSE-FASTQ-001`, `META-BIN-001`, `PROBE-DESIGN-001`) plus a variant-pathogenicity
clinical-caveat exercise. Every cited tool + Method ID was verified with
`find-tool.py`; guarded units + their `MinimumMode` verified against `LIMITATIONS.md`.

## Where this fits

The acceptance instrument for the [[skills-strategy]] and the [[skill-layer]]; enforces
the [[scientific-rigor]] discipline; task style follows the worked workflows in
[[mcp-readme]].
