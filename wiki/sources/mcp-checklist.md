---
type: source
title: "MCP Implementation Checklist v4 — the (superseded) MCP build tracker"
tags: [mcp, planning, superseded]
doc_path: docs/mcp-checklist.md
source_commit: 15ced726eacf8af040cedb6969d34719c4e1009c
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# MCP Implementation Checklist v4 — the (superseded) MCP build tracker

A per-tool build-and-Definition-of-Done tracker for exposing the genomics engine through the **MCP door** — the third of the [[three-front-doors]]. It is a planning artifact: a hierarchy of DoD gates, per-server task lists, a 241-row per-tool checklist, and five quality gates for standing up an MCP tool surface. `docs/mcp-checklist.md` is referenced by path; the long checkbox lists are not copied here.

## Superseded — read this first

The document carries a banner (dated **2026-07-01**) declaring itself **SUPERSEDED**: it tracks a **12-server / 241-tool design that was never built**. Live per-tool status (Built / Tested / Documented) for the **real 11 servers** lives in `docs/mcp/MCP_STATUS.md` — the ledger, not this checklist, is authoritative. The instruction is explicit: do not tick items here.

So this page documents the plan's *shape and conventions*, which partly survived into the shipped servers, while treating its counts and per-tool status as historical.

## What the plan defined

- **Layered Definition of Done.** Three nested DoD bars: **Tool DoD** (MethodId linked; tool/server name frozen; input+output schema; error mapping; 2 tests; `.mcp.json` + `.md`; traceability row), **Server DoD** (all tools pass; `Program.cs`; builds; integration test; README), and **Delivery DoD** (all servers pass; gates G1–G5; cross-server tests; release notes).
- **Two tests per tool, by design.** Every tool gets exactly two minimal tests — a **Schema** test (`*_Schema_ValidatesCorrectly`) and a **Binding** test (`*_Binding_InvokesSuccessfully`) — deliberately *no business-logic asserts* (gate G5.3 audits for that; the algorithm correctness is validated elsewhere — see [[validation-and-testing]]). **Superseded on this point:** the current operational workflow ([[mcp-prompt]]) *reverses* this, requiring the Binding test to assert the algorithm's exact documented values so a wrong wrapper fails.
- **Error-code catalog (1000–5999).** A shared numeric error taxonomy: e.g. `1001` empty text, `1002` null pattern, `1003` empty second sequence, `1004` length mismatch, `1006` invalid nucleotide, `1007` invalid k, `1008` invalid IUPAC bases, `1009` invalid reading frame, `2001` invalid sequence, `2002` invalid protein, `2003` invalid RNA.
- **Traceability to the C# surface.** Each tool row pins a **MethodId** (e.g. `SuffixTree.Contains`, `SequenceStatistics.CalculateMeltingTemperature`), a **HasDocs** flag, and a **DocRef** pointing at either the XML doc in source (`File.cs:Lnn#xml`) or an algorithm markdown page under `docs/algorithms/**`. This is the plan's bridge to the [[mcp-methods-audit|C# method inventory]].
- **Quality gates G1–G5.** Coverage (count = 241, per-server distribution), Docs-first ordering, Traceability-matrix completeness, Documentation-file completeness, and Tests (482 = 241 × 2, all passing, no business asserts).

## Planned server / tool distribution (historical)

12 servers totaling **241 tools**: Core 12, Sequence 35, Parsers 45, Alignment 15, Variants 22, MolBio 38, Assembly 13, Phylogenetics 10, Population 15, Epigenetics 14, Structure 14, Annotation 8. Only Core (12) and Sequence (35) are fully expanded per-tool with frozen schemas; Parsers is partially expanded; servers 4–12 are listed as tool-name rosters only. This uneven detail is itself a sign the plan was abandoned mid-fill.

## Tool-count reconciliation

This is a third, distinct number in the recurring "how many MCP tools?" question:

- **241 tools / 12 servers** — this superseded *plan*.
- **277 public static methods / 54 classes** — the C# surface census in [[mcp-methods-audit]].
- **427 MCP tools / 11 servers** — the shipped figure headlined in the README (see [[three-front-doors]]).

None of these contradict each other; they count different things at different times. The real, current surface is 11 servers, and per-tool status is tracked in `docs/mcp/MCP_STATUS.md`, not here.

## Naming migration visible in the doc

The plan is mid-rename from the old `SuffixTree.Mcp.*` project namespace to `Seqeron.Mcp.*`: Core is still `SuffixTree.Mcp.Core` and the shared/test projects are `SuffixTree.Mcp.Shared` / `SuffixTree.Mcp.Tests`, while Sequence and Parsers are already `Seqeron.Mcp.Sequence` / `Seqeron.Mcp.Parsers`. A useful marker of the codebase's `SuffixTree` → `Seqeron` rebranding.

## Where this fits

- [[mcp-prompt]] — the current, live per-tool completion workflow that declares this checklist superseded and points to `docs/mcp/MCP_STATUS.md` for real status; it also names the shipped 11 servers and evolves the testing policy.
- [[mcp-plan]] — the sibling *design document* (MCP Implementation Plan v4) for this same superseded plan: it carries the v3→v4 rationale, the full 241-tool inventory across all 12 servers, and the standards (sections 6–8) declared still-valid.
- [[three-front-doors]] — this tracker operationalizes the **MCP door** (schemas, error mapping, deterministic tool calls).
- [[mcp-methods-audit]] — the C# method inventory the checklist's MethodId/DocRef columns trace back to.
- [[validation-and-testing]] — the two-tests-per-tool discipline (Schema + Binding, no business asserts) complements the deeper per-unit algorithm validation.
- [[research-grade-limitations]] — counts and per-tool status here are a superseded snapshot; treat `docs/mcp/MCP_STATUS.md` as ground truth.
