# examples — hand-written Python wrappers over Seqeron MCP tools

Three worked wrappers using the documented `call_tool` / `SeqeronToolClient` API from
`scripts/skills/seqeron_mcp_client.py`. Every tool, its server, arguments and output field are
verified against the tool's `.mcp.json` (linked). These wrap validated tools — they do not
reimplement any algorithm ([`bio-rigor`](../../bio-rigor/SKILL.md)). Repo root from here is
`../../../../`.

## 1. gc_content

Tool `gc_content`, server **Sequence**; input `{"sequence": str}`; output
`{gcContent, gcCount, totalCount}` — verified in
[`docs/mcp/tools/sequence/gc_content.mcp.json`](../../../../docs/mcp/tools/sequence/gc_content.mcp.json).

```python
from scripts.skills.seqeron_mcp_client import call_tool, SeqeronToolError

def gc_content(sequence: str, project_root: str = ".") -> float:
    """GC% (0-100) of a DNA/RNA sequence via Seqeron's `gc_content` tool."""
    res = call_tool("Sequence", "gc_content", {"sequence": sequence}, project_root)
    return res["gcContent"]

if __name__ == "__main__":
    try:
        print(gc_content("ATGC"))    # -> 50.0
    except SeqeronToolError as e:
        raise SystemExit(f"Seqeron envelope/validation error: {e}")
```

## 2. melting_temperature (with the optional arg)

Tool `melting_temperature`, server **Sequence**; input `{"sequence": str, "useWallaceRule": bool =
true}`; output `{tm, unit}` — verified in
[`docs/mcp/tools/sequence/melting_temperature.mcp.json`](../../../../docs/mcp/tools/sequence/melting_temperature.mcp.json).
The optional flag matches the schema default; only send it when overriding.

```python
from scripts.skills.seqeron_mcp_client import call_tool

def melting_temperature(sequence: str, use_wallace_rule: bool = True,
                        project_root: str = ".") -> float:
    """Tm in °C via Seqeron's `melting_temperature` tool (Wallace rule or GC formula)."""
    args = {"sequence": sequence}
    if not use_wallace_rule:                 # default is true; only send when overriding
        args["useWallaceRule"] = False
    res = call_tool("Sequence", "melting_temperature", args, project_root)
    return res["tm"]                          # res["unit"] is "°C"

if __name__ == "__main__":
    print(melting_temperature("ATGCGATCGATCG"))                       # Wallace rule
    print(melting_temperature("ATGCGATCGATCG", use_wallace_rule=False))  # GC formula
```

## 3. A two-tool pipeline (pooled session)

Chain `dna_reverse_complement` → `gc_content`, both on the **Sequence** server, reusing one live
session. `dna_reverse_complement` input `{"sequence": str}` → output `{reverseComplement}` — verified
in
[`docs/mcp/tools/sequence/dna_reverse_complement.mcp.json`](../../../../docs/mcp/tools/sequence/dna_reverse_complement.mcp.json).
(GC% of a sequence and its reverse complement should match — a natural `bio-rigor` cross-check.)

```python
from scripts.skills.seqeron_mcp_client import SeqeronToolClient, SeqeronToolError

def revcomp_then_gc(sequence: str, project_root: str = ".") -> dict:
    """Reverse-complement a DNA sequence, then GC% both — one reused Sequence session."""
    with SeqeronToolClient(project_root=project_root) as client:
        rc = client.call("Sequence", "dna_reverse_complement",
                         {"sequence": sequence})["reverseComplement"]
        gc_fwd = client.call("Sequence", "gc_content", {"sequence": sequence})["gcContent"]
        gc_rc  = client.call("Sequence", "gc_content", {"sequence": rc})["gcContent"]
    return {"reverseComplement": rc, "gcForward": gc_fwd, "gcReverseComplement": gc_rc}

if __name__ == "__main__":
    try:
        out = revcomp_then_gc("ATGCATGC")
        print(out)   # gcForward == gcReverseComplement (cross-check passes)
    except SeqeronToolError as e:
        raise SystemExit(f"Seqeron envelope/validation error: {e}")
```

## Notes carried from every example

- Argument names (`sequence`, `useWallaceRule`) and output fields (`gcContent`, `tm`,
  `reverseComplement`) are copied verbatim from each tool's `.mcp.json` — never invented.
- `SeqeronToolError` is always allowed to propagate; the caller reports it (STOP rule).
- Prefer generating these with `gen-python-stub.py <tool>` (see [`recipes.md`](recipes.md)); the
  hand-written forms above show what the generated wrapper looks like.
