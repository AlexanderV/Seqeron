---
type: source
title: "MCP Methods Audit — the public-static-method inventory"
tags: [api-surface, inventory]
doc_path: docs/MCP-Methods-Audit.md
source_commit: 6abaea7ca1fdd97060f8b0130b193f07c8ffaa0b
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# MCP Methods Audit — the public-static-method inventory

A flat, method-level census of `Seqeron.Genomics` (`docs/MCP-Methods-Audit.md`, audit date 2026-01-23): it enumerates **277 public static methods** across **54 classes** (one file each), each row giving class, method name, and source line number, plus a per-class count summary. It is the raw reconciliation artifact behind the C# API — the second of the [[three-front-doors]] — and the working list for deciding what to expose as MCP tools.

## What it establishes

- **277 public static methods / 54 classes.** The full table (class → method → line) and the per-class rollup live in the doc; they are referenced by path, not copied here.
- **Counting rules that shape the number:**
  - **Overloads count separately** — e.g. `SequenceAligner.GlobalAlign` contributes 4, `Translator.Translate` 3.
  - **Extension methods are included** — `SequenceExtensions` alone adds 12.
  - **SAM-flag one-liners are included** — `SequenceIO` has 16 methods, of which 12 are single-line SAM-flag predicates (`IsPaired`, `IsProperPair`, `IsUnmapped`, …).
  - **`StatisticsHelper` (2 methods) is an internal helper** the audit flags as a candidate for exclusion from the MCP surface.
- **Largest surfaces:** `QualityScoreAnalyzer` (13), `SequenceExtensions` (12), then `PopulationGeneticsAnalyzer` / `SequenceStatistics` / `VcfParser` (11 each), and `FastqParser` / `PrimerDesigner` / `SequenceAligner` / `SequenceComplexity` (10 each).

## Where this fits

- [[three-front-doors]] — this audit quantifies the **C# API door**: 277 static methods over 54 classes.
- [[algorithms-checklist-v2]] — an independent view of the same surface. Its Appendix A is a method index and Appendix D reports class coverage of **44/57 classes (77%)**. The denominators differ from this audit (54 vs 57 classes) — see the caveat below.

## Caveats and discrepancies

- **"Public static methods" ≠ "MCP tools."** This is a C#-surface census (277 static methods); the README's headline **427 MCP tools** counts a different, larger surface (instance methods, per-server tool wrappers, and parameter variants also become tools). The two numbers are not comparable and neither contradicts the other. A third, historical figure — **241 tools across 12 servers** — appears in the superseded [[mcp-checklist]] build plan; it too counts a different thing at a different time (a design that was never built).
- **Class-count denominators disagree.** This audit finds 54 classes exposing public static methods; the checklist's Appendix D coverage table uses 57. The gap is almost certainly an inclusion-criteria difference (static-method-bearing files here vs. all analyzer classes there) rather than a true conflict — worth confirming if the two inventories are ever reconciled.
- **Snapshot, not a contract.** The line numbers and counts are a 2026-01-23 snapshot of an unstable public API (see [[research-grade-limitations]]); they drift as the code changes.
