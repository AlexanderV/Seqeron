# Suffix Tree Contract Freeze (Cycle 2)

## Date
2026-02-25

## Purpose
Lock externally observable behavior and storage invariants before the next controlled refactoring cycle.

## Frozen invariants
1. Persistent binary format compatibility remains v6 (`MAGIC`, header layout, compact/large hybrid semantics, jump-table semantics).
2. `ISuffixTree` API surface and method names stay stable for:
   - `Contains`
   - `FindAllOccurrences`
   - `CountOccurrences`
   - `LongestRepeatedSubstring`
   - `LongestCommonSubstring`
   - `Traverse`
   - `FindExactMatchAnchors`
3. Empty-pattern contracts remain stable:
   - `Contains("") == true`
   - `CountOccurrences("") == Text.Length`
   - `FindAllOccurrences("")` returns every valid start index.
4. Traversal branch-balance contract remains stable: every `EnterBranch` must have a matching `ExitBranch`.
5. MCP tool names remain stable:
   - `suffix_tree_contains`
   - `suffix_tree_count`
   - `suffix_tree_find_all`
   - `suffix_tree_lrs`
   - `suffix_tree_lcs`
   - `suffix_tree_stats`

## Known contract divergence to resolve in this cycle
1. `MaxDepth` semantics differ today between in-memory and persistent implementations.
2. Target for this cycle: unify `MaxDepth` semantics under a single `ISuffixTree` contract and keep parity tests green.

## Guardrails for upcoming phases
1. Mandatory test gates must pass after every phase:
   - `tests/SuffixTree/SuffixTree.Tests`
   - `tests/SuffixTree/SuffixTree.Persistent.Tests`
   - `tests/SuffixTree/SuffixTree.Mcp.Core.Tests`
2. For format/perf-sensitive changes, additionally run:
   - persistent `Category=Parity`
   - in-memory `Category=Performance`
3. Each phase is committed separately with no unrelated changes.
