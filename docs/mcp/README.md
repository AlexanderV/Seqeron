# MCP Servers (Seqeron)

## What is MCP?

MCP (Model Context Protocol) is a standard way for LLM clients to call tools with structured inputs/outputs. In this repo, MCP servers expose Seqeron bioinformatics algorithms and parsers as tools with explicit schemas.

See the MCP spec: https://modelcontextprotocol.io/specification/

## What you get in practice

- **Tool access for LLMs** (no copy‑paste code).
- **Structured inputs/outputs** for reproducible results.
- **Portable integration** across MCP‑aware clients.
- **Local execution** via stdio (no extra services required).

## How to connect a server to your LLM tool

Most MCP clients let you register a stdio server with a **command** and **args**. Example:

```text
command: dotnet
args: ["run", "--project", "Seqeron.Mcp.Sequence"]
```

For HTTP-based MCP, the client usually needs a URL instead of a command.

## How to use in practice

1) Start the server you need (Sequence / Parsers / Core).  
2) Register it in your MCP client.  
3) Use the tool docs to build prompts and tool calls.

## Example workflows (prompts and outputs)

### Example (real workflow): cloning insert QC

**Task:** given an insert in FASTA, report GC% and whether it contains EcoRI (GAATTC) or BamHI (GGATCC) sites (0‑based positions).  
This is a standard cloning QC step to avoid cutting sites and to estimate PCR behavior.

**User prompt (in your LLM chat):**

```
Use tools only; no manual parsing or calculations, and no code or command execution.
I have a cloning insert in FASTA below. Please read the sequence from the FASTA using tools (do not interpret FASTA manually),
then check GC% (report with 2 decimals)
and whether it contains EcoRI (GAATTC) or BamHI (GGATCC) sites (report positions as 0-based).
Return a Markdown table with columns: id, length, gc_percent, EcoRI_sites, BamHI_sites,
using this exact header and separator:
| id | length | gc_percent | EcoRI_sites | BamHI_sites |
|---|---:|---:|---|---|
Format EcoRI_sites and BamHI_sites as JSON arrays (e.g., [4] or []).
Output only the table, no extra commentary.

>seq1
GCGCGAATTCATGGATCCATAT
```

**LLM result:**

```
| id | length | gc_percent | EcoRI_sites | BamHI_sites |
|---|---:|---:|---|---|
| seq1 | 22 | 45.45 | [4] | [12] |
```

**Tools used (in order):**
1) `fasta_parse` → sequence = "GCGCGAATTCATGGATCCATAT"  
2) `gc_content` → gcContent = 45.45, gcCount = 10, totalCount = 22  
3) `suffix_tree_find_all` (pattern = "GAATTC") → positions = [4]  
4) `suffix_tree_find_all` (pattern = "GGATCC") → positions = [12]

### Example (real workflow): PCR primer QC

**Task:** validate primer sequences, compute GC% and Tm, and report the Tm difference.  
This is a routine pre‑screen for primer pairs before PCR.

**User prompt (in your LLM chat):**

```
Use tools only; no manual parsing or calculations, and no code or command execution.
These are PCR primers in FASTA below. Please read the sequences from the FASTA using tools (do not interpret FASTA manually),
confirm each primer is valid DNA (A/C/G/T only), report GC% (2 decimals) and Tm in °C (1 decimal), and the Tm difference for the pair.

>FWD
ATGCGATCGATCGATCGTAG
>REV
GCGCGATCGATCGATCGCAA

Return a Markdown table with columns: id, length, gc_percent, tm_c,
using this exact header and separator:
| id | length | gc_percent | tm_c |
|---|---:|---:|---:|
Then add a blank line and one line: tm_diff_c = |Tm_FWD - Tm_REV| (1 decimal).
Output only the table and tm_diff line, no extra commentary.
```

**LLM result:**

```
| id | length | gc_percent | tm_c |
|---|---:|---:|---:|
| FWD | 20 | 50.00 | 51.8 |
| REV | 20 | 60.00 | 55.9 |

tm_diff_c = 4.1
```

**Tools used (in order):**
1) `fasta_parse` → sequences for FWD/REV  
2) `dna_validate` (FWD) → valid, length = 20  
3) `dna_validate` (REV) → valid, length = 20  
4) `gc_content` (FWD) → 50.00%  
5) `gc_content` (REV) → 60.00%  
6) `melting_temperature` (FWD) → 51.8°C  
7) `melting_temperature` (REV) → 55.9°C

## Servers in this repo (and why split)

Servers are split so you can attach only what you need:

- **Seqeron.Mcp.Sequence** — DNA/RNA/Protein models, statistics, complexity, k‑mers.
- **Seqeron.Mcp.Parsers** — FASTA/FASTQ/GenBank/GFF/VCF/BED/EMBL parsing and utilities.
- **SuffixTree.Mcp.Core** — suffix tree, distances, similarity.

This keeps the tool surface focused, faster to load, and easier to reason about.

## Tool catalog and schemas

Tool schemas and examples live here:

- `docs/mcp/tools/core`
- `docs/mcp/tools/sequence`
- `docs/mcp/tools/parsers`

## Connect to Codex (CLI or IDE)

Codex stores MCP configuration in `~/.codex/config.toml`. The CLI and the IDE extension share this file, so you only configure once. You can also scope MCP servers to a trusted project using `.codex/config.toml`.

### Add servers with the CLI (stdio)

Run from the repo root:

```bash
codex mcp add seqeron-sequence -- dotnet run --project Seqeron.Mcp.Sequence
codex mcp add seqeron-parsers  -- dotnet run --project Seqeron.Mcp.Parsers
codex mcp add seqeron-core     -- dotnet run --project SuffixTree.Mcp.Core
```

Verify:

```bash
codex mcp list
```

### Or add to config.toml

```toml
[mcp_servers.seqeron_sequence]
command = "dotnet"
args = ["run", "--project", "Seqeron.Mcp.Sequence"]

[mcp_servers.seqeron_parsers]
command = "dotnet"
args = ["run", "--project", "Seqeron.Mcp.Parsers"]

[mcp_servers.seqeron_core]
command = "dotnet"
args = ["run", "--project", "SuffixTree.Mcp.Core"]
```

Codex MCP guide: https://developers.openai.com/codex/mcp/
