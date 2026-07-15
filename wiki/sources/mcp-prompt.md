---
type: source
title: "MCP Tool Completion — the current one-tool-per-session subagent prompt"
tags: [mcp, workflow, testing]
doc_path: docs/mcp-prompt.md
source_commit: 23e16e74d4fa795924e0351a3f577fb9c2ddc3ee
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
graph:
  relationships:
    - predicate: contradicts
      object: source:mcp-checklist
      source: mcp-prompt
      evidence: "mcp-prompt DoD §2 requires the Binding test to assert EXACT expected values ('a deliberately wrong wrapper must FAIL them'), whereas mcp-checklist mandates two minimal tests with deliberately NO business-logic asserts (gate G5.3 audits for that)."
      confidence: high
      status: current
    - predicate: contradicts
      object: source:mcp-plan
      source: mcp-prompt
      evidence: "mcp-prompt DoD §2 requires evidence-based Binding tests asserting the algorithm's documented values, reversing the mcp-plan §7 testing policy of two minimal Schema+Binding tests with no business-logic asserts."
      confidence: high
      status: current
---

# MCP Tool Completion — the current one-tool-per-session subagent prompt

The **live, operational** playbook for finishing MCP tools — a self-contained subagent prompt to paste into a fresh session that brings **exactly one** MCP tool to Definition of Done (binding audit + tests + docs), runs the full green-test gate, commits, and stops. Each tool runs in its own cold, isolated context. Unlike its superseded planning cousins [[mcp-plan]] and [[mcp-checklist]], this is the process actually in use. `docs/mcp-prompt.md` is referenced by path; the prompt text is not copied here.

The tool is a **thin wrapper**: the prompt is emphatic that bioinformatics algorithms already exist in `Seqeron.Genomics` and the MCP binding must NOT re-implement them — it validates inputs and calls the real method. This is the [[three-front-doors|MCP door]] discipline made procedural.

## Source of truth: MCP_STATUS.md, not the superseded plan

The prompt names `docs/mcp/MCP_STATUS.md` as the authoritative ledger — the real per-tool **B/T/D** (Built / Tested / Documented) status across **11 real servers** — and explicitly instructs the engineer to **ignore `docs/mcp-plan.md` and `docs/mcp-checklist.md`** (the v4 12-server / 241-tool design that was never built). This settles the recurring tool-count / server-decomposition question from the operating side: the shipped surface is the 11 servers below, tracked in the ledger.

## The real 11 servers (server → project → tools file)

This prompt is the first ingested source to **name the shipped 11-server decomposition** concretely (the [[mcp-plan]] and [[mcp-checklist]] pages only note "11 servers, different decomposition" without enumerating them):

| Server | Project | Tools file |
|---|---|---|
| Core | `SuffixTree.Mcp.Core` | `src/SuffixTree/Mcp/SuffixTree.Mcp.Core/Tools/*.cs` |
| Sequence | `Seqeron.Mcp.Sequence` | `.../Tools/SequenceTools.cs` |
| Parsers | `Seqeron.Mcp.Parsers` | `.../Tools/ParsersTools.cs` |
| Alignment | `Seqeron.Mcp.Alignment` | `.../Tools/AlignmentTools.cs` |
| Analysis | `Seqeron.Mcp.Analysis` | `.../Tools/AnalysisTools.cs` |
| Annotation | `Seqeron.Mcp.Annotation` | `.../Tools/AnnotationTools.cs` |
| Chromosome | `Seqeron.Mcp.Chromosome` | `.../Tools/ChromosomeTools.cs` |
| Metagenomics | `Seqeron.Mcp.Metagenomics` | `.../Tools/MetagenomicsTools.cs` |
| MolTools | `Seqeron.Mcp.MolTools` | `.../Tools/MolToolsTools.cs` |
| Phylogenetics | `Seqeron.Mcp.Phylogenetics` | `.../Tools/PhylogeneticsTools.cs` |
| Population | `Seqeron.Mcp.Population` | `.../Tools/PopulationTools.cs` |

Note how this differs from the superseded 12-server plan: there is **no** dedicated Variants / MolBio / Assembly / Epigenetics / Structure server; instead the annotation-domain and analysis-domain tools consolidate under **Annotation** and **Analysis**, and the single **MolTools** server carries the molecular-biology tools. **Core** still lives under the pre-rebrand `SuffixTree.Mcp.Core` project (a lingering marker of the `SuffixTree` → `Seqeron` rename). Tests live in `tests/Seqeron/Seqeron.Mcp.<Server>.Tests/` (Core under `tests/SuffixTree/...`); per-tool docs in `docs/mcp/tools/<server>/{tool}.md` + `{tool}.mcp.json`.

**Gold-standard reference servers** (already Done, to mirror exactly): **Sequence, Parsers, Core** — with `GcContent` / `DnaValidate` as the model binding, `GcContentTests.cs` as the model test, and `gc_content.md` / `gc_content.mcp.json` as the model docs.

## Definition of Done for an MCP tool (three parts, in order)

This is a **tool-wrapper** DoD, distinct from the algorithm-level [[definition-of-done|test-unit DoD]]:

1. **Binding audit + fix.** Attribute `[McpServerTool(Name = "{tool}", Title = "{Domain} — {Human Title}", ReadOnly = true)]` plus a `[Description(...)]` telling the LLM *when* to call it. `ReadOnly = true` for pure/query tools; omit for file-writing tools. Parameters carry `[Description]`; the return is a structured `record` (in the server's `Models/`), never a tuple or raw string. The body validates inputs (`throw new ArgumentException(..., nameof(x))`) and calls the real `Seqeron.Genomics` method. Errors map to the 1000–5999 ranged catalog (input / format / file-IO / algorithm / limits). MolTools tools currently use the bare `[McpServerTool, Description(...)]` form and must be **normalized** to carry an explicit `Name` (keeping the SDK-derived snake_case name so no client breaks).
2. **Tests (≥2, NUnit, one file).** `{Tool}_Schema_ValidatesCorrectly` asserts the input guards; `{Tool}_Binding_InvokesSuccessfully` invokes a known input and asserts the **exact** values from the algorithm's *documented* behavior (`docs/algorithms/**`, the Genomics XML doc, or `Seqeron.Genomics.Tests`) — **not** whatever the code happens to return. "Tests must be evidence-based: a deliberately wrong wrapper (swapped arg, off-by-one) must FAIL them." (See the contradiction below.)
3. **Docs.** `{tool}.mcp.json` (toolName, serverName, version, stability, methodId, real `docRef`, JSON-Schema 2020-12 input+output, `errors[]`, ≥2 real `examples[]`) **and** `{tool}.md` (Description / Parameters / Returns / Errors / ≥1 worked Example / References).

## Execution flow and discipline

Select the first `☐` tool top-to-bottom in `MCP_STATUS.md` server order (or the orchestrator's explicit `<Server>/<tool>`) → fix binding → write docs from real behavior → write tests → **full green gate** (whole unfiltered `dotnet test` + `dotnet build`, `Failed: 0`, never commit on red, never weaken/skip a test) → flip B/T/D to `☑` and adjust roll-ups → commit only when green with explicit `git add -- <paths>` (never `-A`), message `feat(MCP/<Server>): finish <tool> ...`, no push. If correct expected values can't be established from the algorithm's own docs/tests, mark the tool **⛔ Blocked**, revert, do not commit. A hard operational caveat: **no in-place source-mutation experiments** trusted against an incremental `dotnet test` — a stale `bin/obj` DLL yields false results; use a scratch copy or `--no-incremental`.

## Contradiction flagged: evidence-based tests reverse the "no business asserts" policy

The [[mcp-checklist]] and [[mcp-plan]] both specify two *minimal* tests per tool with **deliberately no business-logic asserts** (algorithm correctness validated elsewhere; checklist gate G5.3 audits for their absence). This prompt **reverses that**: its Binding test must assert the **exact** documented algorithm values, and is explicitly required to fail on a wrong wrapper. This is a genuine policy evolution on the operational side, not a copy error — captured as a typed `contradicts` edge to both superseded pages. It is consistent with the project's broader [[scientific-rigor|rigor]] / [[validation-and-testing|evidence-based testing]] stance.

## Where this fits

- [[mcp-plan]] — the superseded design doc this prompt tells you to ignore (except its still-valid standards); it names MCP_STATUS.md as the live successor authority.
- [[mcp-checklist]] — the superseded build tracker, likewise ignored in favor of the ledger.
- [[three-front-doors]] — this prompt is the build procedure for the **MCP door**.
- [[definition-of-done]] — the algorithm-level test-unit DoD; the tool-wrapper DoD here is a separate, thinner bar for the binding layer.
- [[validation-and-testing]] — the evidence-based Binding-test requirement aligns the MCP wrapper layer with the project's assert-real-values testing culture.
- [[research-grade-limitations]] — per-tool status is governed by `docs/mcp/MCP_STATUS.md`; treat that ledger, not the older plans, as ground truth.
