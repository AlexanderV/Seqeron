# Seqeron Bioinformatics

C#/.NET 8 library for bioinformatics algorithms, sequence utilities, and file format parsing. The repository also includes a high-performance suffix tree implementation and MCP servers that expose the APIs as tools.

[![.NET](https://img.shields.io/badge/.NET-8.0-512BD4)](https://dotnet.microsoft.com/)
[![License](https://img.shields.io/badge/license-MIT-blue)](LICENSE)

## Scope

- Sequence models and utilities (DNA/RNA/Protein types, validation, composition, complexity).
- Algorithms for pattern matching, k-mer analysis, repeat analysis, annotation, and molecular tools.
- Parsers and writers for common genomics formats.
- MCP servers for Core, Sequence, and Parsers toolsets.
- Benchmarks, stress harnesses, and extensive unit tests.

## Repository Layout

```
Seqeron.sln
SuffixTree/                 # Suffix tree implementation
Seqeron.Genomics/           # Core algorithms and sequence model
Seqeron.Mcp.Sequence/       # MCP server: sequence analysis tools
Seqeron.Mcp.Parsers/        # MCP server: parser and format tools
SuffixTree.Mcp.Core/        # MCP server: core algorithms (suffix tree, distances)
SuffixTree.Tests/           # Suffix tree tests
Seqeron.Genomics.Tests/     # Genomics tests
Seqeron.Mcp.Sequence.Tests/ # MCP sequence tool tests
Seqeron.Mcp.Parsers.Tests/  # MCP parser tool tests
SuffixTree.Benchmarks/      # Benchmarks
SuffixTree.Console/         # Stress and verification harness
TestSpecs/                  # Algorithm test specifications
docs/                       # Documentation
```

## Documentation

- Algorithms index: [docs/algorithms/README.md](docs/algorithms/README.md).
- Algorithm areas: [Annotation](docs/algorithms/Annotation), [K-mer Analysis](docs/algorithms/K-mer_Analysis), [Pattern Matching](docs/algorithms/Pattern_Matching), [Repeat Analysis](docs/algorithms/Repeat_Analysis), [Sequence Composition](docs/algorithms/Sequence_Composition), [Molecular Tools](docs/algorithms/Molecular_Tools), [MolTools](docs/algorithms/MolTools).
- Suffix tree algorithm: [Suffix Tree (Ukkonen)](docs/algorithms/Pattern_Matching/Suffix_Tree.md).
- MCP tool docs: [Core](docs/mcp/tools/core), [Sequence](docs/mcp/tools/sequence), [Parsers](docs/mcp/tools/parsers), [Traceability](docs/mcp/traceability.md).
- Plans and coverage: [Algorithm checklist](ALGORITHMS_CHECKLIST_V2.md), [Test specs](TestSpecs), [MCP plan](docs/mcp-plan.md), [MCP server plan](docs/MCP-Server-Plan.md), [MCP methods audit](docs/MCP-Methods-Audit.md), [Genomics implementation plan](Seqeron.Genomics/IMPLEMENTATION_PLAN.md), [Integration plan](Seqeron.Genomics.Integration/INTEGRATION_PLAN.md).

## Build and Test

```bash
dotnet build

dotnet test
```

## License

See [LICENSE](LICENSE).
