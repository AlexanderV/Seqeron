"""Reusable Python MCP-client helper for the Seqeron repo.

What this is
------------
A small, dependency-light helper that lets Python tools call the Seqeron
MCP servers. Each server is a .NET (stdio) MCP server launched with

    dotnet run --project <Project>

from the repo root. This module maps the 11 Seqeron MCP server names (as
they appear in ``docs/skills/_generated/catalog.json`` under the ``server``
field) to their dotnet project, spawns the server as a stdio subprocess,
initializes an MCP session, calls a tool, and returns its result.

Requirements
------------
- ``pip install mcp``   (the official Python MCP SDK; imported lazily -- see below)
- the .NET SDK on PATH  (``dotnet``), so the servers can be built/run.

Neither ``mcp`` nor ``dotnet`` is required merely to *import* this module or
to run its offline self-test; they are only needed to actually call a tool.
The ``mcp`` SDK is imported lazily (inside functions/methods), so the module
top level has no hard third-party imports.

Usage
-----
One-shot (spawns a fresh subprocess per call -- simplest, but pays the
dotnet startup cost every time)::

    from seqeron_mcp_client import call_tool
    result = call_tool("Sequence", "reverse_complement",
                       {"sequence": "ACGT"}, project_root="/path/to/Seqeron")

Pooled (keeps ONE live session per server and reuses it across calls --
use this when you make many calls)::

    from seqeron_mcp_client import SeqeronToolClient
    with SeqeronToolClient(project_root="/path/to/Seqeron") as client:
        a = client.call("Sequence", "gc_content", {"sequence": "ACGT"})
        b = client.call("Alignment", "needleman_wunsch", {...})

Both return the tool's ``structuredContent`` dict (falling back to parsing
the first text-content block as JSON when ``structuredContent`` is absent).

Bio-rigor note
--------------
This helper is a THIN transport wrapper around the already-validated Seqeron
tools. It does NOT reimplement any bioinformatics algorithm -- computation
happens inside the validated C#/MCP tool. If a tool returns an MCP error
(e.g. a guarded-unit / LimitationPolicy limitation, or an input-validation
error), that envelope is SURFACED as a ``SeqeronToolError`` rather than
swallowed, so callers never mistake a refused/limited result for a real one.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

# ---------------------------------------------------------------------------
# Server name -> dotnet project (relative to the repo root / project_root).
# Names are exactly the `server` values in docs/skills/_generated/catalog.json.
# ---------------------------------------------------------------------------
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


class SeqeronToolError(Exception):
    """Raised when an MCP tool call returns an error envelope.

    Carries the tool/transport error ``message`` and an optional ``code``
    (e.g. a LimitationPolicy / guarded-unit code, or a validation code).
    """

    def __init__(self, message: str, code: Optional[str] = None):
        self.message = message
        self.code = code
        super().__init__(message if code is None else f"[{code}] {message}")


# ---------------------------------------------------------------------------
# Internal helpers (require the `mcp` SDK; imported lazily).
# ---------------------------------------------------------------------------
def _project_for(server: str) -> str:
    """Return the dotnet project path for a server name, or raise."""
    try:
        return SERVER_PROJECT[server]
    except KeyError:
        known = ", ".join(sorted(SERVER_PROJECT))
        raise SeqeronToolError(
            f"Unknown MCP server {server!r}. Known servers: {known}",
            code="UNKNOWN_SERVER",
        )


def _server_params(server: str, project_root: str):
    """Build StdioServerParameters that launch the server via `dotnet run`."""
    from mcp import StdioServerParameters  # lazy import

    project = _project_for(server)
    return StdioServerParameters(
        command="dotnet",
        args=["run", "--project", project],
        cwd=project_root,
    )


def _extract_result(call_result: Any) -> Dict[str, Any]:
    """Turn an MCP CallToolResult into a plain dict, surfacing errors.

    - If ``isError`` is set, raise ``SeqeronToolError`` with the envelope text.
    - Prefer ``structuredContent``; otherwise parse the first text block as JSON.
    """
    import json

    # Error envelope -> surface, do not swallow (bio-rigor).
    if getattr(call_result, "isError", False):
        message, code = _envelope_error(call_result)
        raise SeqeronToolError(message, code=code)

    structured = getattr(call_result, "structuredContent", None)
    if structured is not None:
        return dict(structured)

    # Fall back to the first text content block parsed as JSON.
    content = getattr(call_result, "content", None) or []
    for block in content:
        text = getattr(block, "text", None)
        if text is not None:
            try:
                parsed = json.loads(text)
            except (ValueError, TypeError):
                return {"text": text}
            return parsed if isinstance(parsed, dict) else {"result": parsed}

    return {}


def _envelope_error(call_result: Any):
    """Extract (message, code) from an error CallToolResult."""
    import json

    content = getattr(call_result, "content", None) or []
    for block in content:
        text = getattr(block, "text", None)
        if text is None:
            continue
        code = None
        message = text
        try:
            parsed = json.loads(text)
        except (ValueError, TypeError):
            parsed = None
        if isinstance(parsed, dict):
            code = parsed.get("code") or parsed.get("errorCode")
            message = (
                parsed.get("message")
                or parsed.get("error")
                or parsed.get("detail")
                or text
            )
        return message, (str(code) if code is not None else None)
    return "MCP tool call returned an error with no content.", None


async def _acall_tool(server: str, tool: str, arguments: Dict[str, Any],
                      project_root: str) -> Dict[str, Any]:
    """Async: spawn server, init session, call one tool, tear down."""
    from mcp import ClientSession  # lazy import
    from mcp.client.stdio import stdio_client  # lazy import

    params = _server_params(server, project_root)
    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            await session.initialize()
            result = await session.call_tool(tool, arguments or {})
            return _extract_result(result)


def call_tool(server: str, tool: str, arguments: Dict[str, Any],
              project_root: str = ".") -> Dict[str, Any]:
    """One-shot: launch ``server``, call ``tool``, return its result dict.

    Spawns ``dotnet run --project <project>`` (cwd=``project_root``) as a
    stdio MCP subprocess, initializes a session, calls ``tool`` with
    ``arguments``, returns the tool's ``structuredContent`` (or the first
    text block parsed as JSON), then cleanly shuts the subprocess down.

    Raises
    ------
    SeqeronToolError
        If ``server`` is unknown, or the tool returns an error envelope
        (guarded-unit / LimitationPolicy limitation, validation error, ...).
    """
    import asyncio

    _project_for(server)  # fail fast on unknown server before spawning asyncio
    return asyncio.run(_acall_tool(server, tool, arguments, project_root))


class SeqeronToolClient:
    """Pooled MCP client: ONE live session per server, reused across calls.

    Use as a context manager. Each distinct server is launched on first use
    and its session kept open until :meth:`close` (or context exit). This
    avoids paying the dotnet startup cost on every call.

    Example
    -------
        with SeqeronToolClient(project_root="/path/to/Seqeron") as client:
            client.call("Sequence", "gc_content", {"sequence": "ACGT"})
    """

    def __init__(self, project_root: str = "."):
        self.project_root = project_root
        # Lazily created dedicated event loop + background sessions.
        self._loop = None
        # server -> (session, async-exit-stack) live handles.
        self._sessions: Dict[str, Any] = {}

    # -- lifecycle ----------------------------------------------------------
    def _ensure_loop(self):
        import asyncio

        if self._loop is None:
            self._loop = asyncio.new_event_loop()
        return self._loop

    async def _aget_session(self, server: str):
        from contextlib import AsyncExitStack

        from mcp import ClientSession
        from mcp.client.stdio import stdio_client

        if server in self._sessions:
            return self._sessions[server][0]

        params = _server_params(server, self.project_root)
        stack = AsyncExitStack()
        read, write = await stack.enter_async_context(stdio_client(params))
        session = await stack.enter_async_context(ClientSession(read, write))
        await session.initialize()
        self._sessions[server] = (session, stack)
        return session

    async def _acall(self, server: str, tool: str,
                     arguments: Dict[str, Any]) -> Dict[str, Any]:
        session = await self._aget_session(server)
        result = await session.call_tool(tool, arguments or {})
        return _extract_result(result)

    async def _aclose(self):
        # Tear down every session's exit stack (subprocess shutdown included).
        for server in list(self._sessions):
            _session, stack = self._sessions.pop(server)
            try:
                await stack.aclose()
            except Exception:
                # Best-effort shutdown; never mask the primary flow.
                pass

    # -- public API ---------------------------------------------------------
    def call(self, server: str, tool: str,
             arguments: Dict[str, Any]) -> Dict[str, Any]:
        """Call ``tool`` on ``server`` over the pooled (reused) session.

        Same result contract and error surfacing as :func:`call_tool`.
        """
        _project_for(server)  # fail fast on unknown server
        loop = self._ensure_loop()
        return loop.run_until_complete(self._acall(server, tool, arguments))

    def close(self):
        """Shut down all live sessions and their server subprocesses."""
        if self._loop is None:
            return
        try:
            self._loop.run_until_complete(self._aclose())
        finally:
            self._loop.close()
            self._loop = None
            self._sessions = {}

    def __enter__(self) -> "SeqeronToolClient":
        return self

    def __exit__(self, exc_type, exc, tb):
        self.close()
        return False


# ---------------------------------------------------------------------------
# Offline self-test: no `mcp` / `dotnet` needed.
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    import os
    import sys

    # Repo root = two levels up from this file (scripts/skills/<file>).
    repo_root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))))

    assert len(SERVER_PROJECT) == 11, (
        f"expected 11 servers, found {len(SERVER_PROJECT)}")

    missing = []
    for name, rel in sorted(SERVER_PROJECT.items()):
        path = os.path.join(repo_root, rel)
        if not os.path.isdir(path):
            missing.append(f"{name} -> {rel}")
    assert not missing, "missing project dirs:\n  " + "\n  ".join(missing)

    # Unknown-server handling must raise SeqeronToolError.
    try:
        _project_for("NoSuchServer")
    except SeqeronToolError as e:
        assert e.code == "UNKNOWN_SERVER", e.code
    else:  # pragma: no cover
        print("FAIL: unknown server did not raise", file=sys.stderr)
        sys.exit(1)

    print(f"OK: {len(SERVER_PROJECT)} servers, all project dirs exist "
          f"(repo root: {repo_root})")
