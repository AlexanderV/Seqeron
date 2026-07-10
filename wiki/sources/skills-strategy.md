---
type: source
title: "Claude Code skills strategy (docs/skills/STRATEGY.md)"
tags: [architecture, skills]
doc_path: docs/skills/STRATEGY.md
sources:
  - docs/skills/STRATEGY.md
source_commit: b3f950caf701615bb8a0296df6c5368d26dde7ec
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Claude Code skills strategy

The plan of record (DRAFT status, last updated 2026-07-01) for turning the
`Seqeron.Genomics` + `SuffixTree` library and its 11 MCP servers (427 tools) into a
system that reliably solves complex biological tasks — both via MCP orchestration and
via the direct C# API. It is the design document behind the [[skill-layer]].

## Fixed decisions

- **Audience: both profiles** — end user via MCP, and developer via C# API. Every
  domain skill gives a recipe for both modes.
- **Many independent skills, not one hub** — each fires on its own trigger.
- **Anti-drift by construction** — skills *reference* the source of truth and share a
  single generator for their tool tables.

The premise: the bottleneck is not documentation or coverage (already exhaustive) but
three things a reference alone doesn't close — **discovery/routing** (pick the right
tools from 427 without bloating context), **orchestration** (correct step order,
parameters, format hand-off, 0-based coordinates, units), and **scientific validity**
(compute with tools, stay inside the validated envelope, cross-check, keep
provenance). Skills are a thin routing + discipline layer over the source of truth,
delivered through the [[mcp-readme|MCP front door]].

## Design principles

Point-don't-duplicate (link to `docs/mcp/tools/`, `docs/algorithms/`, `LIMITATIONS.md`;
never copy schemas/formulas — they drift); self-updating generated tool tables;
dual-mode (MCP tool names + C# `Method ID`); rigor-by-default delegated to `bio-rigor`
(see [[scientific-rigor]]); progressive disclosure (short `SKILL.md`, on-demand
reference files); envelope-aware STOP on out-of-mode guarded units; reproducible
provenance-carrying output.

## Catalog and anti-drift mechanism

Cross-cutting skills (`bio-rigor`, `seqeron-discovery`, `seqeron-dev`) plus a per-domain
family (`bio-qc`, `bio-alignment`, `bio-assembly`, `bio-annotation`, `bio-phylo-popgen`,
`bio-metagenomics`, `bio-moldesign`, `bio-chromosome`) organized as workflow families
rather than server=skill. A single `gen-catalog` generator reads `docs/mcp/tools/**` +
`MCP_STATUS.md` to emit each domain's `_generated/tools.md`; generation markers wall
off hand-written text; a `check-catalog-fresh` CI test fails on drift; a new MCP tool
is "done" only after the catalog is regenerated.

## Status

All phases marked complete (2026-07-01): Ф0 foundation (`bio-rigor`, `seqeron-discovery`),
Ф1 highest-value domains, Ф2 remaining domains, Ф3 developer profile + polish +
skill-development pipeline, with a byte-identical `.github/skills/` mirror kept in parity
by CI. Acceptance is checked against the [[golden-skills-regression|golden regression set]].

## Where this fits

The design doc under the [[skill-layer]] concept; consumes the [[mcp-readme|MCP front door]];
acceptance-tested by [[golden-skills-regression]]; embodies [[scientific-rigor]].
