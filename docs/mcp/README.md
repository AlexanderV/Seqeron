# Seqeron MCP servers

**Turn any MCP-aware LLM into a bioinformatics analyst.** These servers expose Seqeron's genomics
algorithms and file parsers as [Model Context Protocol](https://modelcontextprotocol.io/specification/)
tools — each with a strict JSON schema, so the model *calls a real, validated algorithm* instead of
guessing a number. Same math as the C# API; the tool call is just a different front door.

## Contents

- [Why MCP here](#why-mcp-here)
- [The servers](#the-servers)
- [Connect a server](#connect-a-server)
  - [Any MCP client (stdio)](#any-mcp-client-stdio)
  - [Codex (CLI or IDE)](#codex-cli-or-ide)
- [Use it: two worked workflows](#use-it-two-worked-workflows)
- [Tool catalog & schemas](#tool-catalog--schemas)

## Why MCP here

- **Real computation, not hallucination** — every tool wraps a tested algorithm; outputs are
  deterministic and reproducible.
- **Structured in and out** — explicit input/output schemas mean the model can't fumble the format.
- **Portable** — works across any MCP-aware client (Claude Code, Codex, Copilot, custom agents).
- **Local & self-contained** — servers run over stdio; no network service, no data leaves your machine.
- **Auditable** — each answer comes with the exact tools called, in order (see the worked workflows).

> **New here?** If you use Claude Code or Copilot, you usually don't touch MCP directly — the
> [Agent Skills](../../.claude/skills) attach the right server *on demand* and keep the 427 tool
> schemas out of the model's context. This guide is for wiring the servers into any MCP client
> yourself, or understanding what runs underneath.

## The servers

Tools are split into focused, per-domain servers so you attach only what a task needs — a smaller
tool surface loads faster and is easier for the model to reason about. **427 tools across 11 servers:**

| Server | What it's for | Tools |
|--------|---------------|------:|
| [`Seqeron.Mcp.Sequence`](../../src/Seqeron/Mcp/Seqeron.Mcp.Sequence) | DNA/RNA/Protein models, composition, complexity, k-mers, Tm | 35 |
| [`Seqeron.Mcp.Parsers`](../../src/Seqeron/Mcp/Seqeron.Mcp.Parsers) | FASTA/FASTQ/GenBank/GFF/VCF/BED/EMBL parsing & utilities | 41 |
| [`Seqeron.Mcp.Alignment`](../../src/Seqeron/Mcp/Seqeron.Mcp.Alignment) | Pairwise & multiple sequence alignment | 22 |
| [`Seqeron.Mcp.Analysis`](../../src/Seqeron/Mcp/Seqeron.Mcp.Analysis) | K-mer, motif, repeat, complexity, comparative & structural genomics | 91 |
| [`Seqeron.Mcp.Annotation`](../../src/Seqeron/Mcp/Seqeron.Mcp.Annotation) | Genes/ORFs/promoters, variants, epigenetics, miRNA, splicing, SVs, transcriptomics | 97 |
| [`Seqeron.Mcp.Phylogenetics`](../../src/Seqeron/Mcp/Seqeron.Mcp.Phylogenetics) | Distances, tree building, phylogenetic statistics | 13 |
| [`Seqeron.Mcp.Population`](../../src/Seqeron/Mcp/Seqeron.Mcp.Population) | Population genetics (Fst, diversity, LD, selection) | 18 |
| [`Seqeron.Mcp.Metagenomics`](../../src/Seqeron/Mcp/Seqeron.Mcp.Metagenomics) | Taxonomic classification & community profiling | 19 |
| [`Seqeron.Mcp.Chromosome`](../../src/Seqeron/Mcp/Seqeron.Mcp.Chromosome) | Chromosome-scale analysis (karyotype, centromere, synteny) | 32 |
| [`Seqeron.Mcp.MolTools`](../../src/Seqeron/Mcp/Seqeron.Mcp.MolTools) | Primer/probe/CRISPR design, codon optimization, restriction | 47 |
| [`SuffixTree.Mcp.Core`](../../src/SuffixTree/Mcp/SuffixTree.Mcp.Core) | Suffix-tree search, edit/Hamming distance, k-mer similarity | 12 |

Per-server rollout status lives in [`MCP_STATUS.md`](MCP_STATUS.md).

## Connect a server

### Any MCP client (stdio)

Most clients register a stdio server as a **command** plus **args**. Point the command at `dotnet`
and the args at the server project:

```jsonc
{
  "command": "dotnet",
  "args": ["run", "--project", "Seqeron.Mcp.Sequence"]
}
```

Then restart the client, and the server's tools appear — names, descriptions, and JSON schemas are
advertised at runtime via MCP discovery, so any MCP-aware client can introspect them directly. Swap
the project name for any server in the table above. (HTTP-based clients take a URL instead of a
command.)

### Codex (CLI or IDE)

Codex keeps MCP config in `~/.codex/config.toml` (shared by the CLI and IDE extension; scope it to a
project with `.codex/config.toml`). Add servers from the repo root with the CLI:

```bash
codex mcp add seqeron-sequence -- dotnet run --project Seqeron.Mcp.Sequence
codex mcp add seqeron-parsers  -- dotnet run --project Seqeron.Mcp.Parsers
codex mcp add seqeron-core     -- dotnet run --project SuffixTree.Mcp.Core
codex mcp list          # verify
```

…or write them straight into `config.toml`:

```toml
[mcp_servers.seqeron_sequence]
command = "dotnet"
args = ["run", "--project", "Seqeron.Mcp.Sequence"]

[mcp_servers.seqeron_parsers]
command = "dotnet"
args = ["run", "--project", "Seqeron.Mcp.Parsers"]
```

Codex MCP guide: <https://developers.openai.com/codex/mcp/>

## Use it: two worked workflows

Real prompts, real outputs, and the **exact tools each one called** — note that the model parses the
FASTA and computes every number *with tools*, never by hand.

### 1. Cloning-insert QC — GC% + restriction sites

A standard pre-cloning check: does the insert carry an EcoRI or BamHI site, and what's its GC%?

**Prompt:**

```
Use tools only; no manual parsing or calculations. I have a cloning insert in FASTA below.
Read the sequence with tools, then report GC% (2 decimals) and any EcoRI (GAATTC) / BamHI
(GGATCC) sites as 0-based positions, as a Markdown table (id, length, gc_percent, EcoRI_sites,
BamHI_sites). Sites as JSON arrays. Output only the table.

>seq1
GCGCGAATTCATGGATCCATAT
```

**Result:**

```
| id   | length | gc_percent | EcoRI_sites | BamHI_sites |
|------|-------:|-----------:|-------------|-------------|
| seq1 |     22 |      45.45 | [4]         | [12]        |
```

**Tools, in order:** `fasta_parse` → `gc_content` (45.45; 10/22) →
`suffix_tree_find_all "GAATTC"` → `[4]` → `suffix_tree_find_all "GGATCC"` → `[12]`.

### 2. PCR primer QC — validity, GC%, Tm, ΔTm

A routine pre-screen of a primer pair before ordering.

**Prompt:**

```
Use tools only. These are PCR primers in FASTA. Read them with tools, confirm each is valid DNA,
report GC% (2 decimals) and Tm in °C (1 decimal) as a Markdown table (id, length, gc_percent,
tm_c), then a line: tm_diff_c = |Tm_FWD - Tm_REV|.

>FWD
ATGCGATCGATCGATCGTAG
>REV
GCGCGATCGATCGATCGCAA
```

**Result:**

```
| id  | length | gc_percent | tm_c |
|-----|-------:|-----------:|-----:|
| FWD |     20 |      50.00 | 51.8 |
| REV |     20 |      60.00 | 55.9 |

tm_diff_c = 4.1
```

**Tools, in order:** `fasta_parse` → `dna_validate` ×2 → `gc_content` ×2 → `melting_temperature` ×2.

## Tool catalog & schemas

Every one of the 427 tools ships a per-tool doc (`{tool}.md`) and machine-readable schema
(`{tool}.mcp.json`) under [`tools/<server>/`](tools), and each server's own `README.md` lists its
full tool table.

[core](tools/core) · [sequence](tools/sequence) · [parsers](tools/parsers) · [alignment](tools/alignment) ·
[analysis](tools/analysis) · [annotation](tools/annotation) · [chromosome](tools/chromosome) ·
[metagenomics](tools/metagenomics) · [moltools](tools/moltools) · [phylogenetics](tools/phylogenetics) ·
[population](tools/population)

Traceability from each tool back to its algorithm and validation unit:
[`traceability.md`](traceability.md).
