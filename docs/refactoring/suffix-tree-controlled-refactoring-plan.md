# Suffix Tree Controlled Refactoring Plan

## Goal
Refactor suffix-tree modules without regressions in correctness, performance profile, or storage format behavior.

## Mandatory gates after each phase
1. `dotnet test tests/SuffixTree/SuffixTree.Tests/SuffixTree.Tests.csproj -c Release --nologo`
2. `dotnet test tests/SuffixTree/SuffixTree.Persistent.Tests/SuffixTree.Persistent.Tests.csproj -c Release --nologo`
3. `dotnet test tests/SuffixTree/SuffixTree.Mcp.Core.Tests/SuffixTree.Mcp.Core.Tests.csproj -c Release --nologo`

## Additional gates for format/perf-sensitive phases
1. `dotnet test tests/SuffixTree/SuffixTree.Persistent.Tests/SuffixTree.Persistent.Tests.csproj -c Release --nologo --filter Category=Parity`
2. `dotnet test tests/SuffixTree/SuffixTree.Tests/SuffixTree.Tests.csproj -c Release --nologo --filter Category=Performance`

## Phases
1. Baseline and guardrail artifacts.
2. Decompose persistent load: extract header read/validation.
3. Normalize text matching contract to avoid concrete-type checks.
4. Split mapped storage responsibilities into dedicated partials/helpers.
5. De-recursify in-memory diagnostics traversal.
6. Split MCP tools by bounded context while preserving external API behavior.
7. Align target-framework and nullable policy for suffix-tree modules/tests.
8. Final cleanup and invariant verification.

## Rules
1. Keep each phase small and reviewable.
2. Commit only changes related to the current phase.
3. Do not merge phase N+1 until phase N gates are green.
4. Preserve binary format v6 compatibility and parity contracts.
