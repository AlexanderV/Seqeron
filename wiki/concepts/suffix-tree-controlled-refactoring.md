---
type: concept
title: "Suffix-tree controlled refactoring — contract freeze + per-cycle verification gates"
tags: [refactoring, methodology, testing, architecture]
sources:
  - docs/refactoring/suffix-tree-contract-freeze-cycle3.md
  - docs/refactoring/suffix-tree-contract-freeze.md
source_commit: c8af090404e38f4567fe9a24aa912242cc1f0500
created: 2026-07-17
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:suffix-tree
      source: suffix-tree-controlled-refactoring
      evidence: "suffix-tree-contract-freeze-cycle3.md — the campaign freezes and then refactors the Ukkonen suffix-tree engine ([[suffix-tree]]): the frozen invariants are its persistent v6 binary format, the ISuffixTree public contract (Contains/FindAllOccurrences/CountOccurrences/LongestRepeatedSubstring/LongestCommonSubstring/Traverse/FindExactMatchAnchors), empty-pattern behaviour, traversal branch-balance, and the six suffix_tree_* MCP tool names."
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:characterization-testing
      source: suffix-tree-controlled-refactoring
      evidence: "suffix-tree-contract-freeze-cycle3.md 'Frozen invariants' + 'Guardrails' — freezing behaviour/format/contract and running the SuffixTree test projects (incl. Category=Parity / Category=Performance) as mandatory per-phase gates is the applied form of [[characterization-testing]]: prove each refactor phase behaviour-preserving against a pinned contract."
      confidence: medium
      status: current
---

# Suffix-tree controlled refactoring

A disciplined, multi-cycle refactoring campaign on the Ukkonen [[suffix-tree]] engine
(`SuffixTree`, `SuffixTree.Persistent`, `SuffixTree.Mcp.Core`). Each cycle opens by
**freezing a contract** — an explicit, written list of the behaviours, formats, and public
surfaces that the refactor must NOT change — and then reshapes internal architecture behind
that frozen contract, proving each phase behaviour-preserving with **mandatory test gates**.
This is [[characterization-testing]] applied at campaign scale: the freeze is the pinned
"as-is", the gates are the invariance check. This page is the **shared campaign record**;
the individual cycle plans, freeze docs, and verification reports reconcile onto it.

## The contract-freeze method

Before touching code in a cycle, write down what stays fixed. The **Cycle 3** freeze
(`suffix-tree-contract-freeze-cycle3.md`, 2026-02-25) locks five classes of invariant:

1. **Persistent binary format** stays **v6-compatible** — `MAGIC`, header layout, the hybrid
   compact/large zones, and jump-table semantics (see the v6 header in [[suffix-tree]]).
2. **`ISuffixTree` public contract** stays stable for `Contains`, `FindAllOccurrences`,
   `CountOccurrences`, `LongestRepeatedSubstring`, `LongestCommonSubstring`, `Traverse`, and
   `FindExactMatchAnchors`.
3. **Empty-pattern behaviour** stays stable: `Contains("") == true`,
   `CountOccurrences("") == Text.Length`, and `FindAllOccurrences("")` returns every valid
   start index.
4. **Traversal branch-balance** stays stable: every `EnterBranch` has a matching `ExitBranch`
   (the invariant the serializer's structural-hash visitor depends on).
5. **MCP suffix-tree tool names** stay stable: `suffix_tree_contains`, `suffix_tree_count`,
   `suffix_tree_find_all`, `suffix_tree_lrs`, `suffix_tree_lcs`, `suffix_tree_stats`.

The frozen surface is deliberately the *externally observable* one — file format, public
interface, empty-input edge cases, and the tool names other systems bind to. Everything
below it is fair game for the cycle.

## What Cycle 3 was allowed to change

With the contract frozen, Cycle 3's architecture focus was internal-only:

1. Reduce `PersistentSuffixTreeBuilder` coupling by **extracting explicit state** and passing
   components in.
2. **De-duplicate** streaming traversal logic into shared algorithms.
3. Improve lifecycle ergonomics for the persistent factory APIs.
4. Keep the MCP external contracts stable while **clarifying bounded contexts** internally.

None of these touch a frozen invariant — they refactor the persistent builder and shared
traversal code while leaving format, interface, empty-pattern semantics, branch-balance, and
tool names byte-for-byte identical.

## The base contract-freeze doc and its `MaxDepth` divergence

The original freeze (`suffix-tree-contract-freeze.md`, dated 2026-02-25, self-titled
"Contract Freeze (Cycle 2)") locks the **same five invariant classes** as the Cycle 3 freeze
above — v6 binary format, the `ISuffixTree` public contract, empty-pattern behaviour,
traversal branch-balance, and the six `suffix_tree_*` MCP tool names — and enforces the same
per-phase gates. What it adds that Cycle 3 does not is an explicit **"known contract
divergence to resolve in this cycle"**: a documented spot where the frozen contract is *not
yet* internally consistent, made a first-class work item rather than an accident to trip over
mid-refactor.

- **The divergence:** `MaxDepth` semantics **differ between the in-memory and persistent
  implementations** — the two backends behind the one `ISuffixTree` do not agree on what
  `MaxDepth` means.
- **The cycle target:** **unify `MaxDepth` under a single `ISuffixTree` contract** while
  **keeping the parity tests green** (the persistent `Category=Parity` in-memory-vs-persistent
  equivalence gate below is exactly what proves the unification did not break either backend).

This is the freeze method used in its harder mode: freezing a contract does not require the
contract to already be self-consistent — a known divergence is written down alongside the
frozen surface and scheduled for reconciliation within the cycle, under the same gates.

## Guardrails (per-phase gates)

The freeze is only enforceable if every phase is checked against it. Cycle 3's guardrails:

- **Mandatory gates after each phase** — the three SuffixTree test projects must pass:
  `tests/SuffixTree/SuffixTree.Tests`, `tests/SuffixTree/SuffixTree.Persistent.Tests`,
  `tests/SuffixTree/SuffixTree.Mcp.Core.Tests`.
- **Format/perf-sensitive phases** additionally run the targeted categories:
  `Category=Parity` in the persistent tests (in-memory vs persistent equivalence) and
  `Category=Performance` in the in-memory tests.
- **One commit per phase; no unrelated changes** in a phase commit — so any regression the
  gates catch bisects to a single, reviewable step.

This mirrors the [[build-quality-gate]] philosophy (a gate that must be green to proceed) but
scoped to behavioural equivalence rather than static analysis.

## Why this matters

The suffix-tree engine has two backends behind one interface and a persistent on-disk format
that other artifacts already depend on ([[suffix-tree]]). A naive refactor risks silently
changing the binary format, an empty-input edge case, or an MCP tool name — breakage that
unit correctness tests would not necessarily surface. Freezing the contract turns those
risks into explicit, testable assertions and makes each cycle's diff safe to bisect.

## Campaign cluster

This concept is the reconciliation target for the six suffix-tree controlled-refactoring
documents under `docs/refactoring/`:

- `suffix-tree-controlled-refactoring-plan.md` — the overall plan.
- `suffix-tree-contract-freeze.md` — the initial contract freeze (ingested here; adds the
  `MaxDepth` divergence work item).
- `suffix-tree-contract-freeze-cycle3.md` — the Cycle 3 freeze (ingested here).
- `suffix-tree-controlled-refactoring-verification.md`,
  `...-cycle2-verification.md`, `...-cycle3-verification.md` — the per-cycle verification
  reports.

As the remaining four are ingested, add each to this page's `sources:`, bump `source_commit`
to HEAD, and enrich only genuinely distinct content (plan rationale, earlier frozen
contracts, and the verification outcomes per cycle).
