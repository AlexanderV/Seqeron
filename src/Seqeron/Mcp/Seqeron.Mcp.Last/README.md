# Seqeron.Mcp.Last

MCP server - **the script mafft with the LAST programs (LAST 1256 port) - separately licensed (GPL-3.0-or-later).**

**Licence: see [LICENSE](LICENSE) (copyleft; NOT the MIT licence of the other Seqeron MCP servers).** The MIT `Seqeron.Mcp.Alignment` tool `mafft_command` is unchanged: it has no LAST installed, so `--lastpair` / `--lastmultipair` stop there with the "Install LAST" text of the script;
this server runs the same script with the port of LAST behind it (inventory: [docs/Validation/review-2026-09/LICENSES.md](../../../../docs/Validation/review-2026-09/LICENSES.md)).

Exposes **1 tool** - callable over MCP, with an explicit JSON input schema and a per-tool doc under [`docs/mcp/tools/last/`](../../../../docs/mcp/tools/last).

## Run

```bash
dotnet run --project Seqeron.Mcp.Last
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Last"]`).

## Tools (1)

| Tool | Description |
|---|---|
| `mafft_command_last` | `mafft [options] in.fa` of MAFFT 7.505 with LAST 1256 installed (`--lastpair`, `--lastmultipair`): returns exit status, standard output, standard error and the files of the run |
