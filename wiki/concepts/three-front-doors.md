---
type: concept
title: "Three front doors (skills / C# API / MCP)"
tags: [architecture]
sources:
  - README.md
source_commit: 6a14170477c9472c0be07e3b7c7f7123e31eddcf
created: 2026-07-09
updated: 2026-07-09
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

Seqeron exposes one algorithm engine through three interchangeable entry points — plain-language [[skill-layer|Agent Skills]], a normal C# API, and [MCP](https://modelcontextprotocol.io) tools — such that the same validated algorithm produces the same result regardless of which door is used. The README frames this as "one engine, three front doors."

## The three doors

1. **Plain-language skills** — describe a biology task; the matching skills load, pick tools, and chain a correct pipeline. See [[skill-layer]].
2. **The C# library** — call the algorithms directly (`using Seqeron.Genomics`), with `TryCreate` validation on sequence types.
3. **MCP integration** — any LLM calls the tools with strict schemas and reproducible outputs; each call is a real algorithm, so results are deterministic and auditable.

## Why it matters

The equivalence is a correctness guarantee, not just convenience: skill recipes are **dual-mode** (the MCP tool and the equivalent C# `Method ID` are identical), so a result never depends on the path taken to reach it. This underpins the project's [[scientific-rigor|rigor-by-construction]] stance.

Concrete entry-point details (code snippets, MCP connection steps) live in `README.md`.
