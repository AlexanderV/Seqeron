---
type: source
title: "MCP servers front door (docs/mcp/README.md)"
tags: [architecture, mcp]
doc_path: docs/mcp/README.md
sources:
  - docs/mcp/README.md
source_commit: b3f950caf701615bb8a0296df6c5368d26dde7ec
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# MCP servers front door

The user-facing guide to Seqeron's Model Context Protocol layer: it turns any
MCP-aware LLM into a bioinformatics analyst by exposing the genomics algorithms and
file parsers as tools with strict JSON schemas, so the model **calls a real,
validated algorithm** instead of guessing a number. Same math as the C# API — the
tool call is just a different front door, one of the [[three-front-doors]].

## Why MCP here

Real computation (every tool wraps a tested algorithm; deterministic, reproducible),
structured I/O (explicit input/output schemas), portability across MCP clients,
local/self-contained execution over stdio (no network, no data leaves the machine),
and auditability (each answer records the exact tools called in order). This is the
delivery mechanism behind the [[scientific-rigor]] discipline — compute with tools,
never by hand.

## Shape: 427 tools across 11 servers

Tools are split into focused per-domain servers so a task attaches only what it
needs (smaller surface loads faster, easier to reason about):

| Server | Tools |
| --- | ---: |
| `Seqeron.Mcp.Sequence` | 35 |
| `Seqeron.Mcp.Parsers` | 41 |
| `Seqeron.Mcp.Alignment` | 22 |
| `Seqeron.Mcp.Analysis` | 91 |
| `Seqeron.Mcp.Annotation` | 97 |
| `Seqeron.Mcp.Phylogenetics` | 13 |
| `Seqeron.Mcp.Population` | 18 |
| `Seqeron.Mcp.Metagenomics` | 19 |
| `Seqeron.Mcp.Chromosome` | 32 |
| `Seqeron.Mcp.MolTools` | 47 |
| `SuffixTree.Mcp.Core` | 12 |

Per-server rollout status lives in the campaign ledger `docs/mcp/MCP_STATUS.md`
(coverage-excluded, not ingested). The tool-count and server naming here are the
current ground truth, superseding the 12-server/241-tool plans in [[mcp-plan]] and
[[mcp-checklist]].

## Connecting and using it

A stdio server registers as `dotnet run --project <ServerProject>`; tool names,
descriptions, and JSON schemas are advertised at runtime via MCP discovery. The
guide walks the generic stdio path and a Codex (`~/.codex/config.toml`) setup, then
shows two worked workflows — cloning-insert QC (GC% + EcoRI/BamHI sites) and PCR
primer QC (validity, GC%, Tm, ΔTm) — each with the exact tool chain, underscoring
that the model parses FASTA and computes every number *with tools*.

Most Claude Code / Copilot users never touch MCP directly: the [[skill-layer]]
attaches the right server on demand and keeps the 427 schemas out of context. Each
tool ships a per-tool doc + schema under `docs/mcp/tools/` (coverage-excluded,
generated), with tool→algorithm→validation traceability in `docs/mcp/traceability.md`.

## Where this fits

The MCP front door alongside the C# API and the [[skill-layer]] ([[three-front-doors]]);
the delivery layer for [[scientific-rigor]]; current successor to the [[mcp-plan]] /
[[mcp-checklist]] design docs. The skills strategy that sits on top is [[skills-strategy]].
