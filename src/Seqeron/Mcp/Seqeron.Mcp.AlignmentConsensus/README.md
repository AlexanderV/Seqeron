# Seqeron.Mcp.AlignmentConsensus

MCP server — **Alignment consensus (DECIPHER ConsensusSequence port) - separately licensed (GPL-3.0-only).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **1 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/alignmentconsensus/`](../../../../docs/mcp/tools/alignmentconsensus).

## Run

```bash
dotnet run --project Seqeron.Mcp.AlignmentConsensus
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.AlignmentConsensus"]`).

## Tools (1)

| Tool | Description |
|------|-------------|
| `generate_decipher_consensus` | Bioconductor DECIPHER `ConsensusSequence` consensus (DNA, RNA or protein). |
