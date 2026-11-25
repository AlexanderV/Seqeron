# Seqeron.Mcp.EmbossSequenceFormats

MCP server — **EMBOSS sequence, format-conversion and sequence-database programs (seqret, union, revseq, transeq, infoseq, ...) - separately licensed (GPL-2.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **42 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/embosssequenceformats/`](../../../../docs/mcp/tools/embosssequenceformats).

## Run

```bash
dotnet run --project Seqeron.Mcp.EmbossSequenceFormats
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.EmbossSequenceFormats"]`).

## Tools (42)

| Tool | Description |
|------|-------------|
| `cutseq` | EMBOSS 6.6.0 cutseq: Removes a section from a sequence, as cutseq of EMBOSS 6.6.0. |
| `dbtell` | EMBOSS 6.6.0 dbtell: the definition of a database of the standard EMBOSS configuration. |
| `degapseq` | EMBOSS 6.6.0 degapseq: Removes the gap characters from sequences, as degapseq of EMBOSS 6.6.0. |
| `descseq` | EMBOSS 6.6.0 descseq: the name and / or the description of the FIRST sequence of a file (any input format of seqret) is replaced by name ... |
| `drfinddata` | EMBOSS 6.6.0 drfinddata: finds the data resources of the DRCAT catalogue by the EDAM data term of their output. |
| `entret` | EMBOSS 6.6.0 'entret': the full text of the entries of a sequence file (or of the entries a query selects) exactly as the sequence reader... |
| `extractseq` | EMBOSS 6.6.0 extractseq: Extracts regions of a sequence, as extractseq of EMBOSS 6.6.0. |
| `infoseq` | EMBOSS 6.6.0 infoseq: Displays basic information about sequences, as infoseq of EMBOSS 6.6.0. |
| `listor` | EMBOSS 6.6.0 listor: Writes a list file of the logical OR / AND / XOR / NOT of two sets of sequences, as listor of EMBOSS 6.6.0. |
| `makenucseq` | EMBOSS 6.6.0 makenucseq: Creates random nucleotide sequences, as makenucseq of EMBOSS 6.6.0. |
| `makeprotseq` | EMBOSS 6.6.0 makeprotseq: Creates random protein sequences, as makeprotseq of EMBOSS 6.6.0. |
| `maskambignuc` | EMBOSS 6.6.0 maskambignuc: Masks all ambiguity nucleotide codes with N, as maskambignuc of EMBOSS 6.6.0. |
| `maskambigprot` | EMBOSS 6.6.0 maskambigprot: Masks all ambiguity amino acid codes with X, as maskambigprot of EMBOSS 6.6.0. |
| `maskseq` | EMBOSS 6.6.0 maskseq: Masks regions of a sequence with a character or lower case, as maskseq of EMBOSS 6.6.0. |
| `msbar` | EMBOSS 6.6.0 msbar: Mutates sequences with point, codon and block mutations, as msbar of EMBOSS 6.6.0. |
| `newseq` | EMBOSS 6.6.0 newseq: Makes a sequence from a name, a description, a type and the residues, as newseq of EMBOSS 6.6.0. |
| `noreturn` | EMBOSS 6.6.0 noreturn: Changes the line ends of a text file, as noreturn of EMBOSS 6.6.0. |
| `nospace` | EMBOSS 6.6.0 nospace: Removes whitespace from a text file, as nospace of EMBOSS 6.6.0. |
| `notab` | EMBOSS 6.6.0 notab: Replaces the tabs of a text file with blanks, as notab of EMBOSS 6.6.0. |
| `notseq` | EMBOSS 6.6.0 notseq: Writes the sequences that are not on a list, as notseq of EMBOSS 6.6.0. |
| `nthseq` | EMBOSS 6.6.0 nthseq: Writes the n-th sequence of a file, as nthseq of EMBOSS 6.6.0. |
| `nthseqset` | EMBOSS 6.6.0 nthseqset: Writes the n-th set of sequences of an input of several sets, as nthseqset of EMBOSS 6.6.0. |
| `pasteseq` | EMBOSS 6.6.0 pasteseq: Inserts one sequence into another, as pasteseq of EMBOSS 6.6.0. |
| `prettyseq` | EMBOSS 6.6.0 prettyseq: Writes a nucleotide sequence and its translation, as prettyseq of EMBOSS 6.6.0. |
| `revseq` | EMBOSS 6.6.0 revseq: Reverses and complements nucleotide sequences, as revseq of EMBOSS 6.6.0. |
| `seqcount` | EMBOSS 6.6.0 seqcount: Counts the sequences of a file, as seqcount of EMBOSS 6.6.0. |
| `seqret` | EMBOSS 6.6.0 seqret: the sequences of a file are read by the sequence readers of the EMBOSS library (ajseqread.c) and written in FASTA or... |
| `seqretsplit` | EMBOSS 6.6.0 seqretsplit: the sequences of a file (any input format of seqret) are written each to a file of its own, name.extension in l... |
| `seqxref` | EMBOSS 6.6.0 seqxref: Reports the database cross-references of sequence entries, as seqxref of EMBOSS 6.6.0. |
| `showdb` | EMBOSS 6.6.0 showdb: the table of the databases of the standard EMBOSS configuration. |
| `showorf` | EMBOSS 6.6.0 showorf: Displays a nucleotide sequence and its translations in pretty format, as showorf of EMBOSS 6.6.0. |
| `shuffleseq` | EMBOSS 6.6.0 shuffleseq: Shuffles sequences keeping their composition, as shuffleseq of EMBOSS 6.6.0. |
| `skipseq` | EMBOSS 6.6.0 skipseq: Reads and writes sequences, skipping the first n, as skipseq of EMBOSS 6.6.0. |
| `splitter` | EMBOSS 6.6.0 splitter: Splits sequences into smaller sequences, as splitter of EMBOSS 6.6.0. |
| `transeq` | EMBOSS 6.6.0 transeq: Translates nucleotide sequences to protein in the chosen frames, as transeq of EMBOSS 6.6.0. |
| `trimest` | EMBOSS 6.6.0 trimest: Removes the poly-A tails of EST sequences, as trimest of EMBOSS 6.6.0. |
| `trimseq` | EMBOSS 6.6.0 trimseq: Trims the ambiguity characters and gaps from the ends of sequences, as trimseq of EMBOSS 6.6.0. |
| `trimspace` | EMBOSS 6.6.0 trimspace: Cleans up whitespace in a text file, as trimspace of EMBOSS 6.6.0. |
| `union` | EMBOSS 6.6.0 union: the sequences of a file (any input format of seqret) are joined into ONE sequence that carries the name, description,... |
| `wordcount` | EMBOSS 6.6.0 wordcount: Counts the words of a given size in sequences, as wordcount of EMBOSS 6.6.0. |
| `write_alignment` | Write an alignment in an EMBOSS -aformat output format (pair, markx0-10, srspair, msf, clustal, mega, nexus, phylip, selex, sam, ...). |
| `yank` | EMBOSS 6.6.0 yank: Writes the USA (the universal sequence address) of a sequence, as yank of EMBOSS 6.6.0. |
