# Seqeron.Mcp.EmbossOpenReadingFrames

MCP server — **EMBOSS open reading frame programs (getorf, plotorf, tcode, checktrans) - separately licensed (GPL-2.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **4 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/embossopenreadingframes/`](../../../../docs/mcp/tools/embossopenreadingframes).

## Run

```bash
dotnet run --project Seqeron.Mcp.EmbossOpenReadingFrames
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.EmbossOpenReadingFrames"]`).

## Tools (4)

| Tool | Description |
|------|-------------|
| `checktrans` | EMBOSS 6.6.0 checktrans: Reports the STOP codons and the open reading frames of a protein translation, as checktrans of EMBOSS 6.6.0. |
| `getorf` | EMBOSS 6.6.0 getorf: Finds and extracts the open reading frames (ORFs) of nucleotide sequences, as getorf of EMBOSS 6.6.0. |
| `plotorf` | EMBOSS 6.6.0 plotorf: Plots the potential open reading frames of a nucleotide sequence, as plotorf of EMBOSS 6.6.0. |
| `tcode` | EMBOSS 6.6.0 tcode: Identifies protein-coding regions of a nucleotide sequence with the Fickett TESTCODE statistic, as tcode of EMBOSS 6.... |
