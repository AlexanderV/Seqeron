---
type: source
title: "MCP Implementation Plan v4 — the (superseded) multi-server design + surviving standards"
tags: [mcp, planning, superseded]
doc_path: docs/mcp-plan.md
source_commit: 374bfab439bbf45712e68c2e4fbb8b255df69e4d
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# MCP Implementation Plan v4 — the (superseded) multi-server design + surviving standards

The architectural **design document** for exposing the genomics engine through the **MCP door** — the third of the [[three-front-doors]]. Where its sibling [[mcp-checklist]] is the per-tool build *tracker*, this is the *plan*: the v3→v4 migration rationale, the complete 241-tool inventory rostered across all 12 planned servers, an implementation roadmap, and the engineering standards. `docs/mcp-plan.md` is referenced by path; the long per-tool tables are not copied here.

## Superseded — read this first

The document opens with a banner (dated **2026-07-01**) declaring itself **SUPERSEDED**: the 12-server / 241-tool design (`SuffixTree.Mcp.*`, with dedicated Variants / Assembly / Epigenetics / Structure servers) was **never built**. The repository actually ships **11 servers / ~427 tools** with a different decomposition, and the live source of truth is `docs/mcp/MCP_STATUS.md`.

The banner carves out one explicit exception: **sections 6–8 (standards, testing policy, documentation policy) "remain valid as standards."** So this page treats the *architecture and counts* as historical while preserving the *standards* as still-current engineering contracts.

## The design rationale (v3 → v4)

The plan documents the pivot from a **single monolithic 254-tool server (v3)** to **12 focused servers (v4)**. The stated driver is context economy, citing MCP best practices: tool definitions alone consume ~5–7% of the context window before any prompt, and a monolithic server cost ~24% context per connection versus ~2–8% per focused server. Users "pick needed servers" instead of all-or-nothing. This is the same context-budget pressure that the shipped system answers differently — via the [[skill-layer]], which keeps tool schemas out of the model's context entirely rather than sharding them into servers.

Separation criteria: domain cohesion, ≤8% context budget per server, independent deployment, minimal cross-server dependencies.

## Complete planned tool inventory (historical)

Unlike the checklist (which fully expanded only Core + Sequence), this plan **rosters all 241 tools across all 12 servers** with per-tool `MethodId`, `HasDocs`, and `DocRef`. The distribution: Core 12, Sequence 35, Parsers 45, Alignment 15, Variants 22, MolBio 38, Assembly 13, Phylogenetics 10, Population 15, Epigenetics 14, Structure 14, Annotation 8. Every listed tool is marked `HasDocs=true` / `stable` / `1.0.0`. Each row's `DocRef` points either at an XML-doc anchor in source (`File.cs:Lnn#xml`) or an algorithm markdown page under `docs/algorithms/**`, making this the fullest map from planned tool → C# method → doc. It traces back to the [[mcp-methods-audit|C# method census]].

Notable: the planned per-server groupings do **not** match the shipped skill/server boundaries — e.g. this plan folds miRNA into an "Epigenetics" server and puts comparative-genomics tools under Phylogenetics/Population, whereas the shipped surface distributes them differently.

## Standards that survived (sections 6–8)

These are the reusable, still-valid part of the document:

- **Naming conventions.** Server `SuffixTree.Mcp.{Domain}` (the old namespace — see rename note below), assembly `suffixtree-mcp-{domain}`, tool `{domain}_{action}` or bare `{action}` (`align_global`, `gc_content`), doc files `{toolName}.mcp.json` + `{toolName}.md`.
- **Schema rules.** JSON Schema **draft 2020-12**; required fields explicit (no implicit); nullability via `["string","null"]`; defaults specified and documented; string enums with descriptions.
- **Error mapping (1000–5999).** Ranged taxonomy: 1000–1999 input validation, 2000–2999 sequence format, 3000–3999 file I/O, 4000–4999 algorithm, 5000–5999 resource limits. (The [[mcp-checklist]] enumerates specific codes within these ranges.)
- **Versioning.** SemVer 2.0 — breaking→major, new optional param→minor, fix→patch, experimental→`0.x.y`.
- **Testing policy.** Minimum **2 tests per tool** — a **Schema** test and a **Binding** test — with extra tests only when justified (per overload, per union input variant, complex validation). Deliberately no business-logic asserts (algorithm correctness is validated elsewhere; see [[validation-and-testing]]). **Note:** the current operational workflow ([[mcp-prompt]]) revises this — its Binding test must assert the algorithm's exact documented values (evidence-based; a wrong wrapper must fail).
- **Documentation contract.** Every tool ships a machine-readable `{toolName}.mcp.json` (name, server, description, methodId, input/output schema, references, examples, version, stability) **and** a human-readable `{toolName}.md` (Description / Parameters / Returns / Errors / Examples / References).
- **Quality gates G1–G5.** Coverage (count = 241, per-server distribution), Docs-first ordering, Traceability-matrix completeness, Documentation-file completeness, Tests. **Definition of Done** is layered Tool → Server → Delivery, matching the checklist.

## Tool-count reconciliation

This plan is one of three distinct numbers in the recurring "how many MCP tools?" question, none of which contradict — they count different things at different times:

- **241 tools / 12 servers** — this superseded *plan* (and its [[mcp-checklist]] tracker).
- **277 public static methods / 54 classes** — the C# surface census in [[mcp-methods-audit]].
- **427 MCP tools / 11 servers** — the shipped figure headlined in the README (see [[three-front-doors]]).

The real, current surface is 11 servers; per-tool status lives in `docs/mcp/MCP_STATUS.md`, not here.

## Naming migration marker

Sections 6 still use the pre-rebrand `SuffixTree.Mcp.*` namespace, a marker of the codebase's `SuffixTree` → `Seqeron` rename (the checklist page captures the mid-rename state where later servers already read `Seqeron.Mcp.*`).

## Where this fits

- [[mcp-prompt]] — the current, live per-tool completion workflow that supersedes this plan's operative status (its standards §6–8 survive; per-tool status now lives in `docs/mcp/MCP_STATUS.md`); it names the shipped 11 servers this plan only alluded to.
- [[mcp-checklist]] — the sibling *build tracker* for this same superseded plan; that page documents the DoD/gate mechanics, this one the design rationale and full inventory.
- [[three-front-doors]] — this plan operationalizes the **MCP door**.
- [[skill-layer]] — solves the same tool-context problem this plan raised (5–7% per tool set), but by keeping schemas out of context rather than sharding servers.
- [[mcp-methods-audit]] — the C# method inventory the plan's per-tool `MethodId`/`DocRef` columns trace back to.
- [[research-grade-limitations]] — the architecture and counts here are a superseded snapshot; treat `docs/mcp/MCP_STATUS.md` as ground truth.
