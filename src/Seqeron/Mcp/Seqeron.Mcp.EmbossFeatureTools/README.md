# Seqeron.Mcp.EmbossFeatureTools

MCP server — **EMBOSS feature programs (showfeat, extractfeat, twofeat, featcopy, ...) - separately licensed (GPL-2.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **9 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/embossfeaturetools/`](../../../../docs/mcp/tools/embossfeaturetools).

## Run

```bash
dotnet run --project Seqeron.Mcp.EmbossFeatureTools
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.EmbossFeatureTools"]`).

## Tools (9)

| Tool | Description |
|------|-------------|
| `coderet` | EMBOSS 6.6.0 'coderet': the coding sequences, messenger RNAs, translations and non-coding regions of the sequences of a file with feature... |
| `extractfeat` | EMBOSS 6.6.0 'extractfeat': the features of the sequences of a file (any input format of seqret with features) that pass the source, type... |
| `featcopy` | placeholder |
| `featmerge` | placeholder |
| `featreport` | EMBOSS 6.6.0 'featreport': the first feature table of a feature file reported for a sequence in a report format. |
| `feattext` | placeholder |
| `maskfeat` | EMBOSS 6.6.0 maskfeat: the regions of the features of the wanted types are masked in the sequences of a file (any input format of seqret ... |
| `showfeat` | EMBOSS 6.6.0 'showfeat': the features of the sequences of a file as lines of ASCII graphics under a scale. |
| `twofeat` | EMBOSS 6.6.0 'twofeat': pairs of features of a sequence in a required relation (overlap, distance, sense) reported. |
