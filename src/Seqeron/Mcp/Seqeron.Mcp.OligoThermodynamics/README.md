# Seqeron.Mcp.OligoThermodynamics

MCP server — **Primer3-style primer and probe design and ntthal thermodynamic structure analysis (port of Primer3 thal.c, dpal.c, masker.c, libprimer3) - separately licensed (GPL-2.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The tool names and behaviour are the ones the tools had on the MIT `Seqeron.Mcp.Alignment` / `Seqeron.Mcp.Analysis` servers before the ports of GPL / LGPL programs were separated into their own assemblies
(inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **11 tools** — the same validated algorithms as the C# API of the assemblies this server references, callable over MCP. Every tool carries an explicit JSON input/output schema and a Schema+Binding test, with a
per-tool doc under [`docs/mcp/tools/oligothermodynamics/`](../../../../docs/mcp/tools/oligothermodynamics).

## Run

```bash
dotnet run --project Seqeron.Mcp.OligoThermodynamics
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.OligoThermodynamics"]`).

## Tools (11)

| Tool | Description |
|------|-------------|
| `analyze_oligo` | Analyze the basic physical properties of a short oligonucleotide. |
| `design_antisense_probes` | Design antisense probes against an mRNA-sense sequence. |
| `design_molecular_beacon` | Design a hairpin molecular-beacon probe for real-time detection. |
| `design_primers` | Design a forward/reverse PCR primer pair flanking a target region. |
| `design_probes` | Design ranked hybridization probes for a target sequence. |
| `design_probes_primer3` | Pick hybridization probes exactly as Primer3 does for `PRIMER_TASK=pick_hyb_probe_only`. |
| `design_tiling_probes` | Tile a target with fixed-length, overlapping probes for full coverage. |
| `evaluate_primer` | Evaluate a single primer against quality criteria. |
| `generate_primer_candidates` | Enumerate and evaluate every primer candidate in a region. |
| `primer_dimer` | Primer3 alignment-mode 3′-end primer-dimer check between two primers (PRIMER_PAIR_COMPL_END). |
| `validate_probe` | Validate a probe's specificity against reference sequences. |
