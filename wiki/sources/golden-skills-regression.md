---
type: source
title: "Golden skills regression set (docs/skills/golden/README.md)"
tags: [skills, testing]
doc_path: docs/skills/golden/README.md
sources:
  - docs/skills/golden/README.md
source_commit: a54ba17b2ffb2125fad9712de4ad2cea84ae74a8
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-18
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

## The 10 skills under test

The set exercises **10 skills**: two cross-cutting (`bio-rigor`, `seqeron-discovery`)
and eight domain skills (`bio-qc`, `bio-alignment`, `bio-assembly`, `bio-annotation`,
`bio-moldesign`, `bio-phylo-popgen`, `bio-metagenomics`, `bio-chromosome`).

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
GC-skew/origin — including 4 cross-domain chains (G1, G6, G10, G12) and 3 guarded-unit
tasks (G2 → `PARSE-FASTQ-001`, G8 → `META-BIN-001`, G11 → `PROBE-DESIGN-001`) plus a
variant-pathogenicity clinical-caveat exercise (G4). Every cited tool + Method ID was
verified with `find-tool.py`; guarded units + their `MinimumMode` verified against
`LIMITATIONS.md`. One provenance anomaly is flagged: `detect_rearrangements` exists on
**both** the Analysis server (`ComparativeGenomics.DetectRearrangements`) and the
Chromosome server (`ChromosomeAnalyzer.DetectRearrangements`); the chromosome-scale
synteny task (G10) uses the Chromosome one, matching `bio-chromosome`.

## Where this fits

The acceptance instrument for the [[skills-strategy]] and the [[skill-layer]]; enforces
the [[scientific-rigor]] discipline; task style follows the worked workflows in
[[mcp-readme]].
