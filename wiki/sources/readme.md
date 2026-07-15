---
type: source
title: "README — Seqeron project overview"
tags: [overview]
doc_path: README.md
sources:
  - README.md
source_commit: c4df9b3a138915b4803ed7b97fdeab6658963d04
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-15
---

# README — Seqeron project overview

The repository's front-page `README.md`: it frames Seqeron as a from-scratch bioinformatics toolkit for .NET 10 (250+ algorithms) that can be driven in plain English through AI *agent skills*, while the same algorithms remain a normal C# API and a set of MCP tools. Its recurring thesis is that every number is **computed by the library, never guessed**, and carries provenance.

## What it establishes

- **One engine, three front doors** — the same validated algorithm answers whether invoked through a plain-language skill, a C# call, or an MCP tool. See [[three-front-doors]].
- **A skill routing + discipline layer** — 21 Agent Skills that discover the right tool among 427, orchestrate a correct multi-step pipeline, and keep tool schemas out of the model's context. See [[skill-layer]].
- **A strictly-layered library** — modules across Levels 0–4 with a str* dependency rule enforced by architecture tests. See [[layered-architecture]].
- **Rigor by construction** — a runtime `LimitationPolicy` guards each algorithm's validated scope; results are tool-computed with provenance. See [[scientific-rigor]].
- **A ten-methodology validation campaign** — 22,000+ executed test cases plus a per-unit validation ledger. See [[validation-and-testing]].
- **Honest status** — beta, research-grade, not for clinical or diagnostic use, with many algorithms documented as simplified subsets. See [[research-grade-limitations]].
- **An LLM Wiki retrieval layer** — repository documentation remains authoritative, while a provenance-tracked Markdown index, BM25 search, backlinks, and a typed graph give agents bounded discovery context and exact routes back to source files.

Headline facts cited by the README: **427 MCP tools** across **11 MCP servers**, **~250 algorithms** over **258 algorithm units**, **47 projects**, and a Ukkonen suffix tree substrate. Its versioned LLM Wiki section records **529 curated pages**, **4,594 wikilinks**, and a compiled graph of **529 nodes / 4,214 edges**, while explicitly framing the numbers as a context-size comparison rather than a correctness benchmark. Detailed sub-topics (performance/NativeAOT, the worked triage example, repository layout) live in `README.md` directly and are referenced by path rather than copied here.

## Where this fits

Concept and gotcha pages this source touches:
- [[three-front-doors]]
- [[skill-layer]]
- [[layered-architecture]]
- [[scientific-rigor]]
- [[validation-and-testing]]
- [[research-grade-limitations]]
