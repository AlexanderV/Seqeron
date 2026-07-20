---
type: concept
title: "Skill routing + discipline layer"
tags: [architecture, skills]
sources:
  - README.md
source_commit: 32523eb3576cde5f7a1713f1b168e5808435e3ec
created: 2026-07-09
updated: 2026-07-20
graph:
  relationships:
    - predicate: relates_to
      object: concept:scientific-rigor
      source: readme
      evidence: "stay scientifically honest (compute with tools — never guess; respect each algorithm's validated envelope; carry provenance)"
      confidence: high
      status: current
---

# Skill routing + discipline layer

A thin layer of Agent Skills that turns the Seqeron library into an agent which solves whole biological tasks rather than making single tool calls. It is one of the [[three-front-doors]] and the reason a user can work in plain language without wrangling schemas.

## The problem it solves

With **427 tools**, attaching every schema drowns an LLM's context. The skill layer keeps tool descriptions **out of the model's context** and instead teaches the model to:

- **discover** the right tool (e.g. via `seqeron-discovery`, which searches 427 tools without loading schemas),
- **orchestrate** a correct multi-step pipeline, and
- **stay scientifically honest** — see [[scientific-rigor]].

Every recipe is **dual-mode**: it works whether the MCP tool or the equivalent C# `Method ID` is called, so MCP is not required.

## Shape

The README describes **21 skills** — cross-cutting ones (`seqeron-setup`, `seqeron-discovery`, `bio-rigor`, `seqeron-dev`, `seqeron-python-client`) plus per-domain skills (`bio-qc`, `bio-alignment`, `bio-annotation`, and the rest). They live under `.claude/skills/` with a byte-identical mirror under `.github/skills/`, and an auto-generated catalog plus a CI guardrail keep them in sync with the tools. The plan of record is [[skills-strategy]], acceptance-tested by the [[golden-skills-regression|golden regression set]]; the MCP delivery layer beneath the skills is [[mcp-readme]], and the full tool→concept map of all 427 tools is [[mcp-tool-catalog]].
