# Seqeron.Mcp.EmbossRepeats

MCP server — **EMBOSS word-match and tandem-repeat programs (wordmatch, seqmatchall, etandem, equicktandem) - separately licensed (GPL-2.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **3 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/embossrepeats/`](../../../../docs/mcp/tools/embossrepeats).

## Run

```bash
dotnet run --project Seqeron.Mcp.EmbossRepeats
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.EmbossRepeats"]`).

## Tools (3)

| Tool | Description |
|------|-------------|
| `find_tandem_repeats_emboss` | EMBOSS etandem / equicktandem tandem-repeat programs (diagonal-score models), ported exactly. |
| `seqmatchall` | EMBOSS seqmatchall: all-against-all comparison of a set of sequences by exact word matches. |
| `wordmatch` | EMBOSS 6.6.0 wordmatch: every exact match of at least wordSize symbols between a set of pattern sequences (a FASTA text) and large sequen... |
