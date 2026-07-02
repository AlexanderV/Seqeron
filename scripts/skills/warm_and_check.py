#!/usr/bin/env python3
"""Build all Seqeron MCP servers and verify the on-demand path works.

Onboarding helper. The skills call the shipped ``Seqeron.Mcp.*`` servers **on
demand** (spawn, call, tear down) via :mod:`seqeron_mcp_client` — they are never
registered in an MCP config, so their 427 tool schemas never enter the model's
context. But the *build* has to happen once so the first real call is instant
and any compile error surfaces at setup time, not mid-task.

This driver reuses the client's own build/cache logic (single source of truth,
identical cache key), so warming here means the client's first call finds the
DLL already built and skips straight to running it.

What it does
------------
1. Build every server in :data:`seqeron_mcp_client.SERVER_PROJECT` into the same
   temp cache the client uses (``--full`` also pings ``tools/list`` on each so we
   know it actually *starts*, not just compiles).
2. Run one real tool call (``gc_content``) end-to-end and check the number.

Exit code is non-zero on the first failure, with a human-readable reason, so
``setup.sh`` can gate on it.

Usage::

    python3 scripts/skills/warm_and_check.py [--full] [--repo-root DIR]
"""

from __future__ import annotations

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import seqeron_mcp_client as smc  # noqa: E402  (after sys.path tweak)


def _default_repo_root() -> str:
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def build_all(repo_root: str, full: bool) -> int:
    """Build (and optionally start-check) every server. Returns a failure count."""
    servers = sorted(smc.SERVER_PROJECT)
    width = max(len(s) for s in servers)
    failures = 0
    for i, server in enumerate(servers, 1):
        label = f"  [{i:>2}/{len(servers)}] {server:<{width}}"
        t0 = time.time()
        try:
            dll = smc._server_dll(server, repo_root)  # builds if missing
            if full:
                sess = smc._Session(server, repo_root)
                try:
                    tools = sess._rpc("tools/list", {}).get("tools", [])
                finally:
                    sess.close()
                note = f"ok  ({len(tools)} tools, {time.time() - t0:4.1f}s)"
            else:
                built = "built" if time.time() - t0 > 1.0 else "cached"
                note = f"ok  ({built}, {os.path.basename(dll)})"
            print(f"{label}  {note}", flush=True)
        except smc.SeqeronToolError as e:
            failures += 1
            print(f"{label}  FAILED: {e}", file=sys.stderr, flush=True)
    return failures


def smoke_test(repo_root: str) -> bool:
    """One real end-to-end call through a live server. Returns True on success."""
    seq = "GCGCGAATTCATGGATCCATAT"  # the README cloning-QC example insert
    try:
        result = smc.call_tool("Sequence", "gc_content", {"sequence": seq}, repo_root)
    except smc.SeqeronToolError as e:
        print(f"  smoke test FAILED: {e}", file=sys.stderr)
        return False
    gc = result.get("gcContent", result.get("gc_percent", result.get("value")))
    try:
        ok = abs(float(gc) - 45.45) < 0.01
    except (TypeError, ValueError):
        ok = False
    if ok:
        print(f"  smoke test ok      gc_content({seq}) = {gc}")
        return True
    print(f"  smoke test FAILED: expected gc_content 45.45, got {gc!r} "
          f"(full result: {result})", file=sys.stderr)
    return False


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--full", action="store_true",
                    help="also start each server and list its tools (slower, thorough)")
    ap.add_argument("--repo-root", default=_default_repo_root(),
                    help="repository root (default: inferred from this script)")
    args = ap.parse_args()
    repo_root = os.path.abspath(args.repo_root)

    print(f"Building {len(smc.SERVER_PROJECT)} MCP servers "
          f"({'full start-check' if args.full else 'compile + cache'})...")
    failures = build_all(repo_root, args.full)
    if failures:
        print(f"\n{failures} server(s) failed to build. See errors above.",
              file=sys.stderr)
        return 1

    print("\nVerifying the on-demand path end-to-end...")
    if not smoke_test(repo_root):
        return 1

    print("\nAll servers built and the on-demand MCP path is live.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
