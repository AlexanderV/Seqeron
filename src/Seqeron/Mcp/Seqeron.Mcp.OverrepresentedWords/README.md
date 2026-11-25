# Seqeron.Mcp.OverrepresentedWords

MCP server — **Over-represented words and spaced dyads (RSAT oligo-analysis / dyad-analysis port) - separately licensed (AGPL-3.0-only).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **3 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/overrepresentedwords/`](../../../../docs/mcp/tools/overrepresentedwords).

## Run

```bash
dotnet run --project Seqeron.Mcp.OverrepresentedWords
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.OverrepresentedWords"]`).

## Tools (3)

| Tool | Description |
|------|-------------|
| `dyad_analysis` | RSAT dyad-analysis: over-represented spaced dyads (monad pairs with a fixed spacing). |
| `oligo_analysis` | RSAT oligo-analysis k-mer over-representation with binomial significance. |
| `shared_motifs_significance` | RSAT matching-sequence significance of k-mers shared by several sequences. |
