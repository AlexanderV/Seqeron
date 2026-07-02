---
name: seqeron-setup
description: >-
  One-time install & configuration for a freshly cloned Seqeron repo, so the
  bioinformatics skills are ready to solve real tasks. Use when the user says
  "install and configure", "set up Seqeron", "get me started", "onboard me",
  "prepare the repo", "why aren't the bio tools working", or asks how to make the
  Seqeron MCP tools / skills available. Verifies the .NET 10 SDK and Python are
  present, builds all 11 MCP servers into the on-demand cache (so the first real
  tool call is instant), and runs a live smoke test — without registering MCP
  anywhere (tool schemas stay out of the model's context; skills pull tools in on
  demand). Run this once per clone/machine; after it, just describe a biology task.
allowed-tools: Bash, Read
---

# Seqeron Setup

The one command a new user runs after cloning: **build everything, verify the
path, then get out of the way.** After it, the user describes a biology task in
plain language and the matching domain skill (`bio-qc`, `bio-alignment`,
`bio-annotation`, …) loads itself and runs the real, validated algorithms.

## The model this sets up

- The **20 skills** live in [`.claude/skills/`](../) (Claude Code) and a mirror
  in `.github/skills/` (Copilot / VS Code). Claude Code **auto-discovers** them —
  no install step. They are a thin routing + rigor layer, not the algorithms.
- The algorithms run in the **11 shipped MCP servers** (`Seqeron.Mcp.*` +
  `SuffixTree.Mcp.Core`, 427 tools). Skills call them **on demand** — spawn the
  one server a task needs, call it, tear it down — via the stdlib client
  [`scripts/skills/seqeron_mcp_client.py`](../../../scripts/skills/seqeron_mcp_client.py).
- **No MCP registration.** We deliberately do **not** write a `.mcp.json` or run
  `claude mcp add`. Registering all 427 tools would load every schema into the
  model's context and drown it. On-demand keeps context lean and lets skills pull
  in only what the current task needs. The **only** thing that must happen up
  front is the **build**, so the first on-demand call is instant.

So "setup" = check the toolchain, build all servers, prove one real call works.

## When to use

- A **freshly cloned** repo, or a new machine, before the first biology task.
- The user says "install and configure", "set up Seqeron", "onboard me".
- A bio task just failed with a build/`dotnet`/server-launch error → re-run this;
  it's idempotent and will rebuild what's stale.

Do **not** run it before every task — once per clone/machine is enough.

## How to run

Run the idempotent orchestrator from the repo root:

```bash
scripts/setup.sh
```

It (1) checks the .NET 10 SDK and Python 3 are on PATH, (2) builds all 11 MCP
servers into the same temp cache the on-demand client uses, and (3) runs a live
`gc_content` smoke test end-to-end. Add `--full` to also start each server and
list its tools (slower, but confirms every server actually launches):

```bash
scripts/setup.sh --full
```

First run compiles everything (a few minutes). Later runs are cached and fast.

## Reading the result

- **Ready.** → done. Tell the user they can now just describe a biology task; the
  right skill will load and compute with real numbers. Offer a first example
  (e.g. "Read this FASTA and report GC% and any EcoRI sites").
- **`dotnet` not found / wrong version** → point them to the .NET 10 SDK
  (https://dotnet.microsoft.com/download); the repo targets .NET 10.
- **`python3` not found** → the on-demand client needs Python 3 (stdlib only,
  3.9+). Any system Python 3 works.
- **A server failed to build** → surface the compiler output the script printed;
  that's a real build error in that server's project, not a setup bug.

## After setup

Nothing else to configure. Point the user at a real task and let the domain
skills route it. If they want to know *which* tool does something without loading
schemas, that's [`seqeron-discovery`](../seqeron-discovery/SKILL.md); for the
scientific-honesty rules every result must follow, that's
[`bio-rigor`](../bio-rigor/SKILL.md).
