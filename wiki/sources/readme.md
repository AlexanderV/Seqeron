---
type: source
title: "README — Seqeron project overview"
tags: [overview]
doc_path: README.md
sources:
  - README.md
source_commit: 585e15dfe3a54881f74e86e74b5015b7cd4ec14f
ingested: 2026-07-20
created: 2026-07-09
updated: 2026-07-20
---

# README — Seqeron project overview

The repository's front-page `README.md`: it frames Seqeron as a from-scratch bioinformatics toolkit for .NET 10 (250+ algorithms) that can be driven in plain English through AI *agent skills*, while the same algorithms remain a normal C# API and a set of MCP tools. Its recurring thesis is that every number is **computed by the library, never guessed**, and carries provenance.

## What it establishes

- **One engine, three front doors** — the same validated algorithm answers whether invoked through a plain-language skill, a C# call, or an MCP tool. See [[three-front-doors]].
- **A skill routing + discipline layer** — 21 Agent Skills that discover the right tool among 427, orchestrate a correct multi-step pipeline, and keep tool schemas out of the model's context. See [[skill-layer]].
- **A strictly-layered library** — modules across Levels 0–4 with a strict dependency rule enforced by architecture tests. See [[layered-architecture]].
- **Rigor by construction** — a runtime `LimitationPolicy` guards each algorithm's validated scope; results are tool-computed with provenance. See [[scientific-rigor]].
- **A ten-methodology validation campaign** — 22,000+ executed test cases plus a per-unit validation ledger. See [[validation-and-testing]].
- **Honest status** — beta, research-grade, not for clinical or diagnostic use, with many algorithms documented as simplified subsets. See [[research-grade-limitations]].
- **An LLM Wiki retrieval layer** — repository documentation remains authoritative, while a provenance-tracked Markdown index, BM25 search, backlinks, and a typed graph give agents bounded discovery context and exact routes back to source files. A static conceptual map presents the overlap of biological meaning, computational methods, and evidence/limits without treating pages as rigid categories. Its maintenance commands enforce the parser's business rules and the configured branch-coverage floor.

Headline facts cited by the README: **427 MCP tools** across **11 MCP servers**, **~250 algorithms** over **258 algorithm units**, **47 projects**, and a Ukkonen suffix tree substrate. Its versioned LLM Wiki section records **546 curated pages**, **5,059 wikilinks**, and a compiled graph of **546 nodes / 4,489 edges**, while explicitly framing the volume numbers as a context-size comparison rather than a correctness benchmark. A separate fixed 30-intent bilingual benchmark compares BM25 retrieval `Hit@1/3/10` over the authoritative source corpus (gold = any local `sources:` document) against concept-page retrieval with the wiki, for direct English, direct Ukrainian, and manually reviewed Ukrainian-to-English normalization; it measures retrieval, not answer correctness. Detailed sub-topics (performance/NativeAOT, the worked triage example, repository layout) live in `README.md` directly and are referenced by path rather than copied here.

## Where this fits

Concept and gotcha pages this source touches:
- [[three-front-doors]]
- [[skill-layer]]
- [[layered-architecture]]
- [[scientific-rigor]]
- [[validation-and-testing]]
- [[research-grade-limitations]]
