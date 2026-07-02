"""Reusable, dependency-free Python MCP client for the Seqeron servers.

What this is
------------
A small helper that lets Python tools call the shipped Seqeron MCP servers
(``Seqeron.Mcp.*`` / ``SuffixTree.Mcp.Core``) **through the real MCP protocol**
— JSON-RPC over stdio — with **no third-party packages** (Python stdlib only)
and **without registering MCP anywhere** (the server is spawned on demand and
torn down after). Computation happens inside the validated MCP tool; this is a
thin transport wrapper, never a re-implementation.

Why not the ``mcp`` SDK?
------------------------
The official ``mcp`` SDK needs Python >= 3.10 and a network install. To keep
this usable everywhere (and to keep MCP genuinely optional), we speak the stdio
JSON-RPC protocol directly with the standard library.

Launching a server (important)
------------------------------
The shipped servers set ``<PublishAot>true</PublishAot>``, which disables
reflection-based JSON by default — so ``dotnet run`` / a plain Debug DLL
*crash at startup* (the MCP SDK builds tool schemas via reflection). This helper
therefore builds each server once in a reflection-enabled (JIT) configuration
into a cache dir and runs that DLL. Requires the .NET SDK (``dotnet``) on PATH.
Override the whole launch with the env var ``SEQERON_MCP_CMD_<SERVER>`` (e.g.
point it at an AOT-published native binary): its value is run via the shell.

Usage
-----
One-shot (spawn + call + tear down)::

    from seqeron_mcp_client import call_tool
    r = call_tool("Sequence", "gc_content", {"sequence": "ACGT"},
                  project_root="/path/to/Seqeron")

Pooled (one live session per server, reused across calls)::

    from seqeron_mcp_client import SeqeronToolClient
    with SeqeronToolClient(project_root="/path/to/Seqeron") as c:
        a = c.call("Sequence",  "gc_content",   {"sequence": "ACGT"})
        b = c.call("Alignment", "global_align", {"sequence1": "...", "sequence2": "..."})

Both return the tool's ``structuredContent`` dict, falling back to the first
text-content block parsed as JSON (Seqeron servers return their JSON there).
An error envelope (guarded-unit / LimitationPolicy limitation, validation error)
is surfaced as :class:`SeqeronToolError`, never swallowed (bio-rigor).
"""

from __future__ import annotations

import hashlib
import json
import os
import queue
import subprocess
import tempfile
import threading
import time
from typing import Any, Dict, Optional

# Server name -> dotnet project (relative to project_root). Names match the
# `server` field in docs/skills/_generated/catalog.json.
SERVER_PROJECT: Dict[str, str] = {
    "Sequence": "src/Seqeron/Mcp/Seqeron.Mcp.Sequence",
    "Parsers": "src/Seqeron/Mcp/Seqeron.Mcp.Parsers",
    "Alignment": "src/Seqeron/Mcp/Seqeron.Mcp.Alignment",
    "Analysis": "src/Seqeron/Mcp/Seqeron.Mcp.Analysis",
    "Annotation": "src/Seqeron/Mcp/Seqeron.Mcp.Annotation",
    "Chromosome": "src/Seqeron/Mcp/Seqeron.Mcp.Chromosome",
    "Metagenomics": "src/Seqeron/Mcp/Seqeron.Mcp.Metagenomics",
    "MolTools": "src/Seqeron/Mcp/Seqeron.Mcp.MolTools",
    "Phylogenetics": "src/Seqeron/Mcp/Seqeron.Mcp.Phylogenetics",
    "Population": "src/Seqeron/Mcp/Seqeron.Mcp.Population",
    "Core": "src/SuffixTree/Mcp/SuffixTree.Mcp.Core",
}

PROTOCOL_VERSION = "2024-11-05"


class SeqeronToolError(Exception):
    """Raised when an MCP tool call returns an error envelope, or launch fails."""

    def __init__(self, message: str, code: Optional[str] = None):
        self.message = message
        self.code = code
        super().__init__(message if code is None else f"[{code}] {message}")


def _project_for(server: str) -> str:
    try:
        return SERVER_PROJECT[server]
    except KeyError:
        known = ", ".join(sorted(SERVER_PROJECT))
        raise SeqeronToolError(
            f"Unknown MCP server {server!r}. Known servers: {known}",
            code="UNKNOWN_SERVER",
        )


def _server_dll(server: str, project_root: str) -> str:
    """Build the server once (reflection-enabled JIT) and return its DLL path.

    Cached under the system temp dir keyed by (project_root, server) so the repo
    is never polluted. Rebuilds only when the DLL is missing.
    """
    project = _project_for(server)
    proj_dir = os.path.join(project_root, project)
    name = os.path.basename(project)
    key = hashlib.sha1(os.path.abspath(project_root).encode()).hexdigest()[:12]
    out = os.path.join(tempfile.gettempdir(), "seqeron-mcp", key, server)
    dll = os.path.join(out, f"{name}.dll")
    if os.path.exists(dll):
        return dll
    csproj = os.path.join(proj_dir, f"{name}.csproj")
    proc = subprocess.run(
        ["dotnet", "build", csproj, "-c", "Release",
         "-p:PublishAot=false", "-p:JsonSerializerIsReflectionEnabledByDefault=true",
         "-o", out],
        cwd=project_root, capture_output=True, text=True,
    )
    if proc.returncode != 0 or not os.path.exists(dll):
        raise SeqeronToolError(
            f"Failed to build MCP server {server!r}:\n{proc.stdout[-1500:]}{proc.stderr[-500:]}",
            code="SERVER_BUILD_FAILED",
        )
    return dll


class _Session:
    """One live stdio JSON-RPC session to a single Seqeron MCP server."""

    def __init__(self, server: str, project_root: str):
        cmd_env = os.environ.get(f"SEQERON_MCP_CMD_{server.upper()}")
        if cmd_env:
            args, shell = cmd_env, True
        else:
            args, shell = ["dotnet", _server_dll(server, project_root)], False
        self._p = subprocess.Popen(
            args, cwd=project_root, shell=shell, text=True, bufsize=1,
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
        )
        self._q: "queue.Queue[str]" = queue.Queue()
        threading.Thread(
            target=lambda: [self._q.put(l) for l in self._p.stdout], daemon=True
        ).start()
        self._id = 0
        self._rpc("initialize", {
            "protocolVersion": PROTOCOL_VERSION, "capabilities": {},
            "clientInfo": {"name": "seqeron-mcp-client", "version": "0.2"}})
        self._notify("notifications/initialized")

    def _write(self, obj: dict) -> None:
        self._p.stdin.write(json.dumps(obj) + "\n")
        self._p.stdin.flush()

    def _notify(self, method: str, params: Optional[dict] = None) -> None:
        self._write({"jsonrpc": "2.0", "method": method, "params": params or {}})

    def _rpc(self, method: str, params: dict, timeout: float = 90.0) -> dict:
        self._id += 1
        rid = self._id
        self._write({"jsonrpc": "2.0", "id": rid, "method": method, "params": params})
        end = time.time() + timeout
        while time.time() < end:
            try:
                line = self._q.get(timeout=timeout).strip()
            except queue.Empty:
                break
            if not line or line[0] not in "{[":   # skip host log lines on stdout
                continue
            try:
                msg = json.loads(line)
            except ValueError:
                continue
            if msg.get("id") != rid:
                continue
            if "error" in msg:
                err = msg["error"]
                raise SeqeronToolError(err.get("message", str(err)),
                                       code=str(err.get("code")))
            return msg.get("result", {})
        raise SeqeronToolError(f"Timed out waiting for {method!r}", code="TIMEOUT")

    def call(self, tool: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        result = self._rpc("tools/call", {"name": tool, "arguments": arguments or {}})
        if result.get("isError"):
            raise SeqeronToolError(*_envelope(result))
        return _payload(result)

    def close(self) -> None:
        try:
            self._p.terminate()
        except Exception:
            pass


def _payload(result: dict) -> Dict[str, Any]:
    structured = result.get("structuredContent")
    if isinstance(structured, dict):
        return structured
    for block in result.get("content") or []:
        text = block.get("text") if isinstance(block, dict) else None
        if text is not None:
            try:
                parsed = json.loads(text)
            except ValueError:
                return {"text": text}
            return parsed if isinstance(parsed, dict) else {"result": parsed}
    return {}


def _envelope(result: dict):
    for block in result.get("content") or []:
        text = block.get("text") if isinstance(block, dict) else None
        if text is None:
            continue
        try:
            parsed = json.loads(text)
        except ValueError:
            parsed = None
        if isinstance(parsed, dict):
            code = parsed.get("code") or parsed.get("errorCode")
            msg = parsed.get("message") or parsed.get("error") or text
            return msg, (str(code) if code is not None else None)
        return text, None
    return "MCP tool call returned an error with no content.", None


def call_tool(server: str, tool: str, arguments: Dict[str, Any],
              project_root: str = ".") -> Dict[str, Any]:
    """One-shot: launch ``server``, call ``tool``, return its result dict."""
    s = _Session(server, project_root)
    try:
        return s.call(tool, arguments)
    finally:
        s.close()


class SeqeronToolClient:
    """Pooled client: one live session per server, reused across calls."""

    def __init__(self, project_root: str = "."):
        self.project_root = project_root
        self._sessions: Dict[str, _Session] = {}

    def call(self, server: str, tool: str, arguments: Dict[str, Any]) -> Dict[str, Any]:
        _project_for(server)
        if server not in self._sessions:
            self._sessions[server] = _Session(server, self.project_root)
        return self._sessions[server].call(tool, arguments)

    def close(self) -> None:
        for s in self._sessions.values():
            s.close()
        self._sessions = {}

    def __enter__(self) -> "SeqeronToolClient":
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


# Offline self-test: no `dotnet` needed (import + wiring only).
if __name__ == "__main__":
    import sys

    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    assert len(SERVER_PROJECT) == 11, f"expected 11 servers, found {len(SERVER_PROJECT)}"
    missing = [f"{n} -> {r}" for n, r in sorted(SERVER_PROJECT.items())
               if not os.path.isdir(os.path.join(repo_root, r))]
    assert not missing, "missing project dirs:\n  " + "\n  ".join(missing)
    try:
        _project_for("NoSuchServer")
    except SeqeronToolError as e:
        assert e.code == "UNKNOWN_SERVER", e.code
    else:  # pragma: no cover
        print("FAIL: unknown server did not raise", file=sys.stderr)
        sys.exit(1)
    print(f"OK: {len(SERVER_PROJECT)} servers, all project dirs exist (repo root: {repo_root})")
