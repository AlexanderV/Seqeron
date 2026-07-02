# recipes — patterns for Python tools that call Seqeron over MCP

Fuller patterns for [`SKILL.md`](../SKILL.md). The Python side is a thin wrapper over a **validated**
tool; scientific discipline (provenance, cross-check, the STOP rule, the alpha/clinical caveat) is
delegated to [`bio-rigor`](../../bio-rigor/SKILL.md) — not restated. All paths below are repo-root
relative from here as `../../../../`.

## Prerequisites (recap)

- `pip install mcp` (the SDK is imported lazily inside the helper — installing it late is fine).
- `dotnet` on PATH. Each server is launched as `dotnet run --project <Project>`; the first call to a
  server pays a JIT/restore startup cost, which is exactly why the pooled client exists.

## The helper contract

`scripts/skills/seqeron_mcp_client.py` (built by a sibling task — reference it by this agreed path):

- `call_tool(server, tool, arguments, project_root=".") -> dict` — one-shot: start the server for
  `server`, open one stdio session, call `tool` with `arguments`, return the parsed result dict, tear
  the server down. Convenient for a single call; wasteful in a loop.
- `class SeqeronToolClient` — pooled. Holds **one live stdio session per server** and reuses it
  across calls (`client.call(server, tool, arguments)`). Use as a context manager so sessions close
  cleanly. Reuse across many calls to avoid per-call `dotnet` startup.
- `SERVER_PROJECT` — dict mapping the 11 server names to their dotnet project paths, so callers pass
  a **server name** (e.g. `"Sequence"`), not a project path.
- `SeqeronToolError` — raised on a guarded-unit / limitation / validation error returned by a tool.

## Server → project map

You address a server by **name**; the helper's `SERVER_PROJECT` resolves it to
`dotnet run --project <Project>`. The 11 server names match the tool-doc directories under
[`docs/mcp/tools/`](../../../../docs/mcp/tools/): `core, sequence, parsers, alignment, analysis,
annotation, chromosome, metagenomics, moltools, phylogenetics, population` (the helper keys are the
capitalized server names, e.g. `Sequence`, `MolTools`, as recorded in each tool's `serverName`). Get
the right server for a tool from `find-tool.py` / `catalog.json` or the tool's `.mcp.json`
`serverName` field — never assume it.

## One-shot vs pooled

```python
from scripts.skills.seqeron_mcp_client import call_tool, SeqeronToolClient, SeqeronToolError

# One-shot — a single call, server started and stopped for you.
res = call_tool("Sequence", "gc_content", {"sequence": "ATGC"})
print(res["gcContent"])   # 50.0

# Pooled — many calls, one live session per server reused across them.
with SeqeronToolClient(project_root=".") as client:
    for s in ("ATGC", "GGGCCC", "ATATAT"):
        print(client.call("Sequence", "gc_content", {"sequence": s})["gcContent"])
```

Rule of thumb: **1 call → `call_tool`; a loop or a multi-step pipeline → `SeqeronToolClient`.** A
pipeline that touches two servers keeps one session per server open inside the same `with` block.

## Error / envelope handling

A tool can return a guarded-unit / limitation / validation error; the helper raises it as
`SeqeronToolError`. **Surface it — do not swallow it.** Widening a limitation mode just to force a
number is the STOP rule violation from [`bio-rigor`](../../bio-rigor/SKILL.md).

```python
try:
    res = call_tool("Sequence", "gc_content", {"sequence": user_input})
except SeqeronToolError as e:
    # Report the limitation + its documented alternative; do NOT retry with a looser mode.
    raise SystemExit(f"Seqeron refused this input (envelope/validation): {e}")
```

## Batching many tools

For a fan-out over inputs, or a chain of different tools, open the pooled client once:

```python
with SeqeronToolClient(project_root=".") as client:
    def gc(seq):  return client.call("Sequence", "gc_content", {"sequence": seq})["gcContent"]
    def tm(seq):  return client.call("Sequence", "melting_temperature", {"sequence": seq})["tm"]
    rows = [(s, gc(s), tm(s)) for s in sequences]   # one Sequence session for all of them
```

Different servers in the same block are fine — the client keeps one session per distinct server.

## Codegen example

Prefer generating the wrapper from the tool's schema so the signature stays honest:

```bash
python3 scripts/skills/gen-python-stub.py gc_content            # emits a typed gc_content(...) wrapper
python3 scripts/skills/gen-python-stub.py melting_temperature   # includes the optional useWallaceRule arg
```

`gen-python-stub.py` (built by a sibling task) reads `docs/mcp/tools/<server>/<tool>.mcp.json` and
emits a Python function that calls `call_tool` with the schema's argument names/types and returns the
documented output. If you hand-write instead, copy the argument names verbatim from the `.mcp.json`
`inputSchema.properties` — a typo there is a silent bug.

## Provenance

Each wrapper's result carries a `bio-rigor` provenance block: the tool name(s) + server(s) + the
arguments passed, in call order. The Python layer forwards; the validated tool computes.
