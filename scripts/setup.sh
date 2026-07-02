#!/usr/bin/env bash
#
# Seqeron onboarding — build everything the Claude Code / Copilot skills need and
# verify the on-demand MCP path works. Idempotent: safe to run any number of times.
#
# What it does NOT do: it does not register any MCP server in a config. The 20
# skills call the shipped Seqeron.Mcp.* servers on demand (spawn → call → tear
# down), so their 427 tool schemas never enter the model's context. This script
# just makes sure the servers are built (first call is then instant) and the path
# is live.
#
# Usage:
#   scripts/setup.sh            # build all servers + smoke test
#   scripts/setup.sh --full     # also start each server and list its tools (slower)
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FULL=""
[ "${1:-}" = "--full" ] && FULL="--full"

bold() { printf '\033[1m%s\033[0m\n' "$1"; }
ok()   { printf '  \033[32m✓\033[0m %s\n' "$1"; }
warn() { printf '  \033[33m!\033[0m %s\n' "$1"; }
die()  { printf '  \033[31m✗ %s\033[0m\n' "$1" >&2; exit 1; }

bold "Seqeron setup"
echo  "Repo: $REPO_ROOT"
echo

# ── 1. Prerequisites ──────────────────────────────────────────────────────────
bold "1. Checking prerequisites"

if ! command -v dotnet >/dev/null 2>&1; then
  die "The .NET SDK ('dotnet') is not on PATH. Install .NET 10 SDK: https://dotnet.microsoft.com/download"
fi
DOTNET_VER="$(dotnet --version 2>/dev/null || echo '0.0.0')"
DOTNET_MAJOR="${DOTNET_VER%%.*}"
if [ "${DOTNET_MAJOR:-0}" -ge 10 ] 2>/dev/null; then
  ok ".NET SDK $DOTNET_VER"
else
  warn ".NET SDK $DOTNET_VER found, but this repo targets .NET 10. Build may fail — install the .NET 10 SDK."
fi

if ! command -v python3 >/dev/null 2>&1; then
  die "python3 is not on PATH. The on-demand MCP client needs Python 3 (stdlib only; 3.9+ is fine)."
fi
ok "$(python3 --version 2>&1)"

# ── 2. Build every MCP server + verify the on-demand path ─────────────────────
echo
bold "2. Building the MCP servers (first run compiles all 11; later runs are cached)"
if ! python3 "$REPO_ROOT/scripts/skills/warm_and_check.py" $FULL --repo-root "$REPO_ROOT"; then
  die "Build or smoke test failed — see the output above."
fi

# ── 3. Confirm the skills are present ─────────────────────────────────────────
echo
bold "3. Checking skills"
SKILL_COUNT="$(find "$REPO_ROOT/.claude/skills" -maxdepth 2 -name SKILL.md 2>/dev/null | wc -l | tr -d ' ')"
if [ "${SKILL_COUNT:-0}" -gt 0 ]; then
  ok "$SKILL_COUNT skills available in .claude/skills (auto-loaded by Claude Code)"
else
  warn "No SKILL.md files found under .claude/skills — skills may not load."
fi

# ── Done ──────────────────────────────────────────────────────────────────────
echo
bold "Ready."
cat <<'EOF'
  Open this repo in Claude Code (or GitHub Copilot / VS Code) and just describe a
  biology task in plain language — the matching skill loads itself and runs the
  real, validated Seqeron algorithms over the on-demand MCP path. For example:

    "Read this FASTA, report GC% and any EcoRI/BamHI sites."
    "Align these two sequences and tell me what changed."
    "Design a PCR primer pair around this variant."

  Tools are pulled in only as a task needs them, so context stays lean.
EOF
