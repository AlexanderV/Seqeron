# Seqeron.Mcp.EmbossMotifs

MCP server — **EMBOSS pattern and motif programs (fuzzpro, fuzztran, patmatmotifs, patmatdb, prosextract, marscan, consambig) - separately licensed (GPL-2.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **12 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/embossmotifs/`](../../../../docs/mcp/tools/embossmotifs).

## Run

```bash
dotnet run --project Seqeron.Mcp.EmbossMotifs
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.EmbossMotifs"]`).

## Tools (12)

| Tool | Description |
|------|-------------|
| `consambig` | EMBOSS consambig: ambiguity-code consensus of an alignment. |
| `emboss_report_file` | The whole report FILE of the EMBOSS 6.6.0 programs fuzznuc, fuzzpro, fuzztran, etandem and equicktandem on a FASTA text of one or several... |
| `find_known_motifs` | Search a DNA sequence for a set of known motifs simultaneously (exact, or IUPAC / both strands / up to k mismatches). |
| `find_motif` | Exact motif occurrences in a DNA sequence (ascending, forward strand). |
| `find_pattern_hits` | EMBOSS fuzznuc pattern search (classes, exclusions, repeats, ranges, anchors, mismatches, both strands). |
| `fuzzpro` | EMBOSS fuzzpro: protein pattern search (PROSITE-like patterns, exact). |
| `fuzztran` | EMBOSS fuzztran: protein pattern search after translation of the chosen frames. |
| `generate_emboss_consensus` | EMBOSS cons scoring-matrix plurality consensus of an alignment. |
| `marscan` | EMBOSS 6.6.0 marscan: Finds matrix/scaffold recognition (MRS) signatures in DNA sequences, as marscan of EMBOSS 6.6.0. |
| `patmatdb` | EMBOSS 6.6.0 patmatdb: Searches protein sequences with a PROSITE sequence motif, as patmatdb of EMBOSS 6.6.0. |
| `patmatmotifs` | EMBOSS 6.6.0 patmatmotifs: Scans a protein sequence with the patterns of the PROSITE database, as patmatmotifs of EMBOSS 6.6.0. |
| `prosextract` | EMBOSS 6.6.0 prosextract: Builds the PROSITE motif database for patmatmotifs from prosite.dat and prosite.doc, as prosextract of EMBOSS 6... |
