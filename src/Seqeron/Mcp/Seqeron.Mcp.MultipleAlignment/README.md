# Seqeron.Mcp.MultipleAlignment

MCP server — **Multiple sequence alignment programs (Kalign, POA, ClustalW, T-Coffee, ProbCons) - separately licensed (GPL / LGPL).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **25 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/multiplealignment/`](../../../../docs/mcp/tools/multiplealignment).

## Run

```bash
dotnet run --project Seqeron.Mcp.MultipleAlignment
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.MultipleAlignment"]`).

## Tools (25)

| Tool | Description |
|------|-------------|
| `clustalw_command` | ClustalW 2.1 as a function of its command line and the text of the files it names. |
| `clustalw_guide_tree` | ClustalW 2.1 neighbour-joining guide-tree file (.dnd) from pairwise distances. |
| `clustalw_pairwise_align` | ClustalW 2.1 full pairwise stage (FullPairwiseAlign) for DNA or protein pairs. |
| `clustalw_profile_align` | ClustalW 2.1 profile-profile alignment of two existing alignments (clustalw -PROFILE). |
| `clustalw_quick_pairwise_align` | ClustalW 2.1 fast k-tuple pairwise stage (-QUICKTREE). |
| `clustalw_sequences_to_profile` | ClustalW 2.1 alignment of sequences to an existing profile (clustalw -SEQUENCES). |
| `kalign_command` | Kalign 3.3.5 as a function of its command line and the text of the files it names. |
| `multiple_align_clustalw` | ClustalW 2.1's multiple alignment of DNA (clustalw -TYPE=DNA -ALIGN) or PROTEIN (type='protein', clustalw -TYPE=PROTEIN -ALIGN), ported s... |
| `multiple_align_consistency` | Consistency-based MSA (T-Coffee; Notredame, Higgins & Heringa 2000): a primary library of global (Needleman–Wunsch) and local (Smith–Wate... |
| `multiple_align_iterative` | Progressive MSA followed by MUSCLE-style tree-dependent restricted-partitioning refinement (Edgar 2004, stage 3): each guide-tree edge sp... |
| `multiple_align_kalign` | Kalign 3.3.5's multiple sequence alignment, ported from the program. |
| `multiple_align_poa` | POA 2.0 partial-order alignment (port of the program) |
| `multiple_align_probcons` | ProbCons 1.12 (port of the program) |
| `multiple_align_progressive` | Guide-tree progressive MSA (Feng & Doolittle 1987 / Clustal-style): the ClustalW pairwise distance of every pairwise alignment, d = 1 − (... |
| `multiple_align_tcoffee` | The multiple alignment of the T-Coffee 13.41 program (t_coffee; Notredame, Higgins & Heringa 2000) reproduced end to end: its own primary... |
| `poa_command` | command line (FASTA / matrix / pair-score files) |
| `probcons_command` | command line (options, MFA / MSF / Clustal input files) |
| `tcoffee_blast_msa_library` | T-Coffee blastp_msa library built from the XML of an all-against-all BLAST search. |
| `tcoffee_evaluate` | T-Coffee 13.41 evaluation of an alignment (score_ascii) |
| `tcoffee_guide_tree` | The guide tree T-Coffee builds for a set of sequences, as Newick. |
| `tcoffee_ktup_similarity` | The k-tuple similarity matrix of the T-Coffee guide-tree stage. |
| `tcoffee_lalign_top_alignments` | T-Coffee lalign_id_pair: the top-N non-intersecting local alignments (Huang-Miller SIM). |
| `tcoffee_primary_library` | The primary library T-Coffee builds for a set of sequences. |
| `tcoffee_slow_pair_align` | T-Coffee slow_pair: Myers-Miller global alignment of two sequences. |
| `tcoffee_tblastx_library` | T-Coffee tblastx_msa library built from the XML of an all-against-all tblastx search. |
