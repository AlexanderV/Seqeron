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

### `tests/SuffixTree/SuffixTree.Esa.Tests` (ESA backend: in-memory `EsaSuffixTree` and the file-backed tree of `SuffixIndexFactory.Open`)
- `Facade/`
  - The `ISuffixTree`-level assertions of the other two projects re-run through a backend factory (`Ported/`), `FacadeChecks` (every member against the definitions, the Ukkonen tree and the v6 tree), the shared-index tests
    (`FacadeConcurrencyTests`, `FacadeSharedIndexTests`: four threads, 20 repetitions, checksums equal to the single-thread run, once-only build of the navigation pack) and the Dispose race (`FacadeDisposeRaceTests`, a child process).
- `Contract/`, `External/`
  - The index contract (ADR-001) and its conformance suite over the in-memory and the file-backed factories; the container (SESA) builder, the reader, the validation levels and the damaged-file tests.
    The two codecs of the CHILD section (`SesaChildCodecTests`, work package 1b.9a) and the readers of the container are tested on both; the damaged-file tests alternate them.
- `Search/`
  - The pattern descent of the in-memory index against the definitions and the Ukkonen tree (`EsaSearchTests`, `EsaSearchExhaustiveTests`) and the descent over a container (`MappedSearchTests`, 1b.9a: the specialised loop against the in-memory index and the readable reference `EsaIndex.DescendMapped`
    for every packed width, LCP representation, child codec and pattern length, and on damaged containers between guard pages).
- `Fuzz/` (category `Fuzz`, `[Explicit]`, time budget from the environment; the S8 CI job runs the same code)
  - The 1a.7 harness (seeded generator families x alphabet sizes x size bands, minimizer, regression corpus `regressions.txt`) and the full differential fuzz of 1b.8 (`FullFuzzRunner`: every engine, every member of the
    interfaces, the container corruption fuzz, the shared-index test; `FullDifferentialFuzzTests` is the budgeted run). The operation coverage list is `docs/refactoring/suffix-tree-upgrade/analysis/s1b/harness_member_coverage.tsv`
    (one row per interface member; a reflection test in the default run fails when a member has no row).

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
