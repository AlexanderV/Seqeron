---
name: seqeron-python-client
description: >-
  Write a small Python tool/script that CALLS Seqeron algorithms over MCP. Use
  when the task is "wrap a Seqeron MCP tool in Python", "make a python function
  that calls gc_content / this Seqeron tool", "call the Seqeron algorithm from
  python", "generate a python wrapper for <tool>", or "build a python script /
  pipeline using Seqeron over MCP". The wrapper WRAPS the validated tool (never
  reimplement the algorithm in Python) using the repo helper
  scripts/skills/seqeron_mcp_client.py (call_tool / SeqeronToolClient) and the
  codegen scripts/skills/gen-python-stub.py. For the in-process C# API path use
  seqeron-dev; for tool discipline use bio-rigor.
allowed-tools: Bash, Read, Grep, Glob
---

# seqeron-python-client — Python tools that call Seqeron over MCP

For writing a **small Python program that calls a Seqeron algorithm** by driving the MCP server as a
subprocess. The Python side is a thin **wrapper**: it forwards arguments to the validated tool and
returns the tool's result. It never reimplements the algorithm — that would discard the library's
validation and envelope guards (the wrap-don't-reimplement discipline of
[`bio-rigor`](../bio-rigor/SKILL.md); the surface-envelope rules live there, not restated here).

Two Python paths exist. This skill covers **MCP over stdio** (any language client, no CLR). For the
**in-process** path (pythonnet over the C# API, same process) see the C# mechanics in
[`seqeron-dev`](../seqeron-dev/SKILL.md).

## Prerequisites

- **No pip dependencies** — the helper speaks MCP (stdio JSON-RPC) with the Python stdlib only.
- `dotnet` (the .NET SDK) on PATH. The shipped servers set `<PublishAot>true</PublishAot>`, which
  disables reflection-based JSON — so a plain `dotnet run` / Debug DLL **crashes at startup** (the MCP
  SDK builds tool schemas via reflection). The helper therefore builds each server once in a
  reflection-enabled (JIT) config into a temp cache and runs that DLL. `SERVER_PROJECT` maps the 11
  server names → their dotnet projects; override the launch per server with `SEQERON_MCP_CMD_<SERVER>`
  (e.g. point it at an AOT-published native binary).

## The recipe

1. **Find the tool.** Use [`seqeron-discovery`](../seqeron-discovery/SKILL.md) /
   `python3 scripts/skills/find-tool.py <keywords>` or `docs/skills/_generated/catalog.json`
   (tool → server → methodId → doc). Get the exact tool name and its **server**.
2. **Read its contract.** Open `docs/mcp/tools/<server>/<tool>.mcp.json` for the exact input
   properties (names, types, required, defaults) and the output schema. Do **not** guess arg names.
3. **Get the wrapper.** Either generate it —
   `python3 scripts/skills/gen-python-stub.py <tool>` emits a typed function built from that
   `.mcp.json` — **or** hand-write it against `call_tool` (below). Codegen is preferred; it keeps the
   signature honest to the schema.
4. **One-shot vs pooled.** `call_tool(server, tool, arguments, project_root=".")` spins a server,
   makes one call, tears it down — fine for a single call. For **many** calls use
   `SeqeronToolClient`, which holds one live stdio session per server and reuses it — this avoids
   paying `dotnet` startup on every call.
5. **Rigor.** The wrapper wraps the validated tool; it must **surface `SeqeronToolError`**
   (the guarded-unit / limitation / validation envelope), never swallow it. Follow the STOP rule in
   [`bio-rigor`](../bio-rigor/SKILL.md): report the limitation, don't force a number.

The helper is `scripts/skills/seqeron_mcp_client.py` (forthcoming — a sibling task builds it):
`call_tool(...)`, `class SeqeronToolClient`, `SERVER_PROJECT`, and `SeqeronToolError`.

## Worked example — a gc_content wrapper (one-shot)

`gc_content` (server **Sequence**) takes `{"sequence": str}` and returns
`{gcContent, gcCount, totalCount}` — verified in
[`docs/mcp/tools/sequence/gc_content.mcp.json`](../../../docs/mcp/tools/sequence/gc_content.mcp.json).

```python
from scripts.skills.seqeron_mcp_client import call_tool, SeqeronToolError

def gc_content(sequence: str, project_root: str = ".") -> float:
    """GC% of a DNA/RNA sequence, computed by the Seqeron `gc_content` tool."""
    result = call_tool("Sequence", "gc_content", {"sequence": sequence}, project_root)
    return result["gcContent"]   # 0-100, per the tool's outputSchema

try:
    print(gc_content("ATGC"))    # -> 50.0
except SeqeronToolError as e:
    print(f"envelope/validation error: {e}")   # surface it, don't swallow
```

## Worked example — two calls, pooled session

Reuse one live Sequence session across calls instead of restarting `dotnet` each time. Both
`gc_content` and `melting_temperature` live on the **Sequence** server (verified `.mcp.json`).

```python
from scripts.skills.seqeron_mcp_client import SeqeronToolClient

with SeqeronToolClient(project_root=".") as client:
    seq = "ATGCGATCGATCG"
    gc = client.call("Sequence", "gc_content", {"sequence": seq})["gcContent"]
    tm = client.call("Sequence", "melting_temperature", {"sequence": seq})["tm"]  # °C
    print(f"GC%={gc}  Tm={tm}°C")
```

The `SeqeronToolError` envelope still propagates from `client.call`; let it — see the STOP rule.

## Provenance

Report which tool(s) + server(s) + arguments each wrapper called, in order (a `bio-rigor` block). A
Python wrapper does not change provenance: the numbers come from the tool, not from Python.

## Reference

- [`reference/recipes.md`](reference/recipes.md) — server→project map note, one-shot vs pooled,
  error/envelope handling, batching many tools, a codegen example command.
- [`reference/examples.md`](reference/examples.md) — hand-written example wrappers (gc_content,
  melting_temperature, a two-tool pipeline).
- Helper (forthcoming): [`scripts/skills/seqeron_mcp_client.py`](../../../scripts/skills/seqeron_mcp_client.py) ·
  codegen (forthcoming): [`scripts/skills/gen-python-stub.py`](../../../scripts/skills/gen-python-stub.py).
- Discovery: [`scripts/skills/find-tool.py`](../../../scripts/skills/find-tool.py) ·
  catalog [`docs/skills/_generated/catalog.json`](../../../docs/skills/_generated/catalog.json) ·
  tool schemas [`docs/mcp/tools/`](../../../docs/mcp/tools/).
- Rigor: [`bio-rigor`](../bio-rigor/SKILL.md) · discovery: [`seqeron-discovery`](../seqeron-discovery/SKILL.md) ·
  in-process C# path: [`seqeron-dev`](../seqeron-dev/SKILL.md).
