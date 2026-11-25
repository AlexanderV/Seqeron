# Seqeron.Mcp.EmbossPairwiseAlignment

MCP server — **EMBOSS pairwise alignment programs (needle, water, stretcher, matcher, supermatcher, merger, est2genome, diffseq, ...) - separately licensed (GPL-2.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **16 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/embosspairwisealignment/`](../../../../docs/mcp/tools/embosspairwisealignment).

## Run

```bash
dotnet run --project Seqeron.Mcp.EmbossPairwiseAlignment
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.EmbossPairwiseAlignment"]`).

## Tools (16)

| Tool | Description |
|------|-------------|
| `alignment_distance_matrix` | Pairwise distance matrix of a multiple alignment as EMBOSS 6.6.0 distmat (every correction and option, verified on 1 100 + 600 runs of th... |
| `diffseq` | EMBOSS 6.6.0 diffseq: the regions of difference between two nearly identical sequences (336 nucleotide / protein pairs, every -rformat: r... |
| `dotmatcher` | EMBOSS dotmatcher: dot plot segments of two sequences. |
| `est2genome` | EMBOSS 6.6.0 est2genome: the spliced alignment of ESTs to a genomic sequence (exons, introns with donor / acceptor sites, both strands). |
| `matcher_align` | Waterman–Eggert (1987, J Mol Biol 197:723) declumping: the best local alignment, then the best local alignment that shares no ALIGNED RES... |
| `megamerger` | EMBOSS megamerger: merge two large overlapping DNA sequences by word-match anchors. |
| `merger` | EMBOSS merger: merge two overlapping nucleotide sequences. |
| `multiple_alignment_statistics` | Statistics of a gapped multiple alignment: per-sequence comparison with a reference row or the EMBOSS cons consensus (EMBOSS infoalign: s... |
| `needle_align` | EMBOSS 6.6 needle: global alignment with the program's single-precision arithmetic and end-gap options. |
| `needleall_align` | EMBOSS 6.6 needleall: global alignment of every sequence of one set against every sequence of another. |
| `plotcon` | EMBOSS 6.6.0 plotcon: the similarity of a multiple alignment along its length (the data of the plot, -graph data). |
| `showalign` | EMBOSS showalign: pretty display of a multiple alignment. |
| `stretcher_align` | EMBOSS 6.6 stretcher: Myers-Miller linear-space global alignment with integer scores. |
| `supermatcher_align` | EMBOSS 6.6 supermatcher: Rabin-Karp word seeds and a banded Smith-Waterman of large sequences. |
| `water_align` | EMBOSS 6.6 water: local (Smith-Waterman-Gotoh) alignment as the program reports it. |
| `wordfinder` | EMBOSS 6.6.0 wordfinder: Matches large sequences against one or more other sequences by their words, as wordfinder of EMBOSS 6.6.0. |
