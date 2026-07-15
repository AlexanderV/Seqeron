---
type: concept
title: "Three front doors (skills / C# API / MCP)"
tags: [architecture]
sources:
  - README.md
  - docs/MCP-Methods-Audit.md
source_commit: d48405df40409c9f3fd6d5386f3b19c2bf2188df
created: 2026-07-09
updated: 2026-07-15
graph:
  relationships:
    - predicate: relates_to
      object: concept:skill-layer
      source: readme
      evidence: "One engine, three front doors ... The same validated algorithm answers whichever door you walk through."
      confidence: high
      status: current
---

# Three front doors (skills / C# API / MCP)

Seqeron exposes one algorithm engine through three interchangeable entry points — plain-language [[skill-layer|Agent Skills]], a normal C# API, and [MCP](https://modelcontextprotocol.io) tools — such that the same validated algorithm produces the same result regardless of which door is used. The [[readme|project README]] frames this as "one engine, three front doors."

## The three doors

1. **Plain-language skills** — describe a biology task; the matching skills load, pick tools, and chain a correct pipeline. See [[skill-layer]].
2. **The C# library** — call the algorithms directly (`using Seqeron.Genomics`), with `TryCreate` validation on sequence types. A 2026-01-23 census counts **277 public static methods across 54 classes** on this surface — see [[mcp-methods-audit]] (note this is a smaller, different surface from the README's 427 MCP tools).
3. **MCP integration** — any LLM calls the tools with strict schemas and reproducible outputs; each call is a real algorithm, so results are deterministic and auditable.

## Why it matters

The equivalence is a correctness guarantee, not just convenience: skill recipes are **dual-mode** (the MCP tool and the equivalent C# `Method ID` are identical), so a result never depends on the path taken to reach it. This underpins the project's [[scientific-rigor|rigor-by-construction]] stance.

A separate knowledge layer — the repository's LLM Wiki (this wiki) — *complements* the three doors rather than adding a fourth: the doors compute biology results, while the wiki answers questions about how the repository is designed, validated, connected, and constrained. `README.md` frames the division as skills routing biology tasks, MCP executing algorithms, and the wiki as the repository-knowledge layer over the same engine.

Concrete entry-point details (code snippets, MCP connection steps) live in `README.md`.
