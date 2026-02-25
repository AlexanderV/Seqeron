# Suffix Tree Contract Freeze (Cycle 3)

## Date
2026-02-25

## Purpose
Lock behavior and architecture boundaries before the next controlled refactoring cycle.

## Frozen invariants
1. Persistent binary format remains v6-compatible (`MAGIC`, header layout, hybrid compact/large zones, jump-table semantics).
2. `ISuffixTree` public contract remains stable for:
   - `Contains`
   - `FindAllOccurrences`
   - `CountOccurrences`
   - `LongestRepeatedSubstring`
   - `LongestCommonSubstring`
   - `Traverse`
   - `FindExactMatchAnchors`
3. Empty-pattern behavior remains stable:
   - `Contains("") == true`
   - `CountOccurrences("") == Text.Length`
   - `FindAllOccurrences("")` returns every valid start index.
4. Traversal branch-balance remains stable: each `EnterBranch` must have matching `ExitBranch`.
5. MCP suffix-tree tool names remain stable:
   - `suffix_tree_contains`
   - `suffix_tree_count`
   - `suffix_tree_find_all`
   - `suffix_tree_lrs`
   - `suffix_tree_lcs`
   - `suffix_tree_stats`

## Cycle 3 architecture refactoring focus
1. Reduce `PersistentSuffixTreeBuilder` coupling by extracting explicit state and pass components.
2. De-duplicate streaming traversal logic in shared algorithms.
3. Improve lifecycle ergonomics for persistent factory APIs.
4. Keep MCP external contracts stable while clarifying bounded contexts internally.

## Guardrails
1. Mandatory gates after each phase:
   - `tests/SuffixTree/SuffixTree.Tests`
   - `tests/SuffixTree/SuffixTree.Persistent.Tests`
   - `tests/SuffixTree/SuffixTree.Mcp.Core.Tests`
2. For perf/format-sensitive phases, additionally run:
   - `Category=Parity` in persistent tests
   - `Category=Performance` in in-memory tests
3. One commit per phase; no unrelated changes in phase commits.
