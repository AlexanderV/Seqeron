# Suffix Tree Test Hierarchy

## Goal
Keep tests discoverable, non-duplicated, and contract-driven across all suffix-tree modules.

## Resulting structure

### `tests/SuffixTree/SuffixTree.Tests` (in-memory tree)
- `Core/`
  - Construction, invariants, diagnostics (`Traverse`, `PrintTree`, `ToString`).
- `Search/`
  - `Contains`, `FindAllOccurrences`, `CountOccurrences` API behavior.
- `Algorithms/`
  - LRS/LCS/anchors, brute-force and metamorphic checks, algorithmic guards.
- `Performance/`
  - Complexity/perf guardrails only.
- `Properties/`, `Regression/`, `Robustness/`, `Compatibility/`
  - Fuzz/property invariants, regressions, stress/threading, encoding/binary contracts.

### `tests/SuffixTree/SuffixTree.Persistent.Tests` (persistent tree)
- `Core/`
  - Factory/build lifecycle, heap/MMF variants, guard checks.
- `Validation/`
  - Public contract parity and semantic invariants (`Span` overloads, load contracts, diagnostics contracts).
- `Format/`
  - Binary/hybrid/jump-table format behavior.
- `Parity/`
  - Cross-implementation differential tests (in-memory vs persistent, safe vs unsafe paths).
- `Safety/`
  - Disposal, concurrency, provider safety, bounds/overflow behavior.
- `Serialization/`
  - Import/export/hash and truncation contracts.
- `Algorithms/`
  - Algorithm-level fuzzing on persistent backend.

### `tests/SuffixTree/SuffixTree.Mcp.Core.Tests` (MCP tools)
- `Core/`
  - Suffix-tree MCP tools behavior + cross-tool consistency/oracle checks.
- `Genomics/`
  - Genomics-oriented MCP tools (repeat/common-region/similarity/edit/hamming/approximate count).
- `Contracts/`
  - Tool registration-name stability, wrapper/core parity, compatibility constraints.

## Anti-duplication rules
1. Cross-tool consistency checks live in dedicated contract tests (not duplicated per endpoint test).
2. Input guard patterns should be asserted once per contract surface, with representative endpoint cases.
3. Formatter tests assert structural invariants (depth/indent/order/cardinality), not brittle full snapshots.
4. New tests must be placed by concern first (`Validation` vs `Safety` vs `Format`) rather than by implementation detail.
