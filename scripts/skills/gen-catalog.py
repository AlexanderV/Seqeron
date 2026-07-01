#!/usr/bin/env python3
"""Anti-drift catalog generator for the Seqeron Claude Code skills system.

WHY THIS EXISTS
---------------
The skills architecture (docs/skills/STRATEGY.md, esp. sections 5 and 6) is
"many independent skills". The risk of that architecture is drift: N skills
whose tool tables slowly diverge from the real set of 427 MCP tools. The
counter-measure is a SINGLE generator (this script) that reads the source of
truth and emits every tool table, wrapped in machine-owned markers.

SOURCE OF TRUTH (never edited by this script):
  * docs/mcp/tools/<server>/*.md   -- 427 per-tool docs (11 servers)
  * docs/mcp/MCP_STATUS.md         -- (referenced by strategy; not parsed here)

GENERATED (owned by this script; do not hand-edit inside the markers):
  * docs/skills/_generated/catalog.json          -- machine-readable catalog
  * docs/skills/_generated/tool-catalog.md       -- human-readable, per-server
  * .claude/skills/<domain>/_generated/tools.md  -- per-domain slice, ONLY if
                                                    the skill dir already exists

ANTI-DRIFT CONTRACT
-------------------
Generated content lives between:
    <!-- BEGIN generated: do not edit by hand -->
    ...
    <!-- END generated -->
Any hand text outside those markers is left alone (this script fully owns the
target files today, but the markers keep the contract explicit and future-proof).

MODES
-----
  (default)   Write/refresh all generated files in place.
  --check     Generate everything to a temp dir, diff against the committed
              files, print every stale/missing file, exit 1 on ANY drift,
              exit 0 when all are fresh. Never touches committed files.
              This is the CI entrypoint (see check-catalog-fresh.sh).

A new MCP tool is only "done" once this catalog is regenerated and --check is
green (STRATEGY section 5, rule 4).

Dependency-free: python3 stdlib only (target 3.9.6, macOS).
"""

import argparse
import json
import os
import re
import sys
import tempfile

# ---------------------------------------------------------------------------
# Paths & config
# ---------------------------------------------------------------------------

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(SCRIPT_DIR, os.pardir, os.pardir))

TOOLS_ROOT = os.path.join(REPO_ROOT, "docs", "mcp", "tools")
DOMAIN_MAP_PATH = os.path.join(REPO_ROOT, "docs", "skills", "domain-map.json")

GENERATED_DIR = os.path.join(REPO_ROOT, "docs", "skills", "_generated")
CATALOG_JSON = os.path.join(GENERATED_DIR, "catalog.json")
TOOL_CATALOG_MD = os.path.join(GENERATED_DIR, "tool-catalog.md")

CLAUDE_SKILLS_ROOT = os.path.join(REPO_ROOT, ".claude", "skills")

# Lowercase server dirs and their expected tool counts (source of truth).
EXPECTED_COUNTS = {
    "core": 12,
    "sequence": 35,
    "parsers": 41,
    "alignment": 22,
    "analysis": 91,
    "annotation": 97,
    "chromosome": 32,
    "metagenomics": 19,
    "moltools": 47,
    "phylogenetics": 13,
    "population": 18,
}
EXPECTED_TOTAL = 427

BEGIN_MARK = "<!-- BEGIN generated: do not edit by hand -->"
END_MARK = "<!-- END generated -->"

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_H1_RE = re.compile(r"^#\s+(.+?)\s*$")
_SERVER_RE = re.compile(r"^\|\s*\*\*Server\*\*\s*\|\s*(.+?)\s*\|")
_METHOD_RE = re.compile(r"^\|\s*\*\*Method ID\*\*\s*\|\s*`(.+?)`\s*\|")
_SOURCE_RE = re.compile(r"^-\s+Source:\s*\[[^\]]*\]\(([^)]+)\)")


def parse_tool_doc(path):
    """Parse one tool .md into a dict, or return None if it isn't a tool doc.

    Fields: tool, server, methodId, docPath (repo-relative), sourceLink.
    sourceLink is "" when the doc has no '- Source:' line.
    """
    tool = None
    server = None
    method_id = None
    source_link = ""

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if tool is None:
                m = _H1_RE.match(line)
                if m:
                    tool = m.group(1).strip()
                    continue
            if server is None:
                m = _SERVER_RE.match(line)
                if m:
                    server = m.group(1).strip()
                    continue
            if method_id is None:
                m = _METHOD_RE.match(line)
                if m:
                    method_id = m.group(1).strip()
                    continue
            if not source_link:
                m = _SOURCE_RE.match(line)
                if m:
                    source_link = m.group(1).strip()

    if tool is None or server is None or method_id is None:
        return None

    doc_path = os.path.relpath(path, REPO_ROOT)
    return {
        "tool": tool,
        "server": server,
        "methodId": method_id,
        "docPath": doc_path,
        "sourceLink": source_link,
    }


def collect_tools():
    """Walk docs/mcp/tools/<server>/*.md and return (records, per_server_counts).

    records are sorted by (server-dir, tool) for deterministic output.
    per_server_counts is keyed by lowercase server dir name.
    """
    records = []
    per_server = {}

    if not os.path.isdir(TOOLS_ROOT):
        die("tools root not found: %s" % TOOLS_ROOT)

    for server_dir in sorted(os.listdir(TOOLS_ROOT)):
        server_path = os.path.join(TOOLS_ROOT, server_dir)
        if not os.path.isdir(server_path):
            continue
        count = 0
        for name in sorted(os.listdir(server_path)):
            if not name.endswith(".md"):
                continue
            full = os.path.join(server_path, name)
            rec = parse_tool_doc(full)
            if rec is None:
                warn("could not parse tool doc (skipped): %s"
                     % os.path.relpath(full, REPO_ROOT))
                continue
            # Attach the lowercase server DIR for grouping/domain mapping;
            # the 'server' field keeps the doc's display value.
            rec["_serverDir"] = server_dir
            records.append(rec)
            count += 1
        per_server[server_dir] = count

    records.sort(key=lambda r: (r["_serverDir"], r["tool"]))
    return records, per_server


def crosscheck_counts(per_server):
    """Warn loudly (non-fatal) if counts diverge from the expected 427."""
    total = sum(per_server.values())
    problems = []

    for server, expected in EXPECTED_COUNTS.items():
        got = per_server.get(server)
        if got is None:
            problems.append("MISSING server dir '%s' (expected %d tools)"
                            % (server, expected))
        elif got != expected:
            problems.append("server '%s': expected %d, found %d"
                            % (server, expected, got))

    for server in per_server:
        if server not in EXPECTED_COUNTS:
            problems.append("UNEXPECTED server dir '%s' (%d tools)"
                            % (server, per_server[server]))

    if total != EXPECTED_TOTAL:
        problems.append("TOTAL: expected %d, found %d"
                        % (EXPECTED_TOTAL, total))

    if problems:
        warn("=" * 60)
        warn("CATALOG COUNT MISMATCH -- source may have drifted:")
        for p in problems:
            warn("  - " + p)
        warn("=" * 60)
    return total

# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def _relink(from_file, repo_rel_target):
    """Markdown link path from `from_file` (a generated file, absolute path)
    to a repo-relative target, using POSIX separators."""
    target_abs = os.path.join(REPO_ROOT, repo_rel_target)
    rel = os.path.relpath(target_abs, os.path.dirname(from_file))
    return rel.replace(os.sep, "/")


def render_catalog_json(records):
    clean = [
        {
            "tool": r["tool"],
            "server": r["server"],
            "methodId": r["methodId"],
            "docPath": r["docPath"].replace(os.sep, "/"),
            "sourceLink": r["sourceLink"],
        }
        for r in records
    ]
    return json.dumps(clean, indent=2, ensure_ascii=False) + "\n"


def render_tool_catalog_md(records):
    lines = []
    lines.append("# MCP Tool Catalog (generated)")
    lines.append("")
    lines.append("> Generated by `scripts/skills/gen-catalog.py` from "
                 "`docs/mcp/tools/**`. Do not edit inside the markers by hand; "
                 "run the generator instead. Freshness is enforced by "
                 "`scripts/skills/check-catalog-fresh.sh`.")
    lines.append("")
    lines.append(BEGIN_MARK)
    lines.append("")

    # Group by lowercase server dir, deterministic order.
    servers = sorted({r["_serverDir"] for r in records})
    for server_dir in servers:
        rows = [r for r in records if r["_serverDir"] == server_dir]
        display = rows[0]["server"] if rows else server_dir
        lines.append("## %s (%d tools)" % (display, len(rows)))
        lines.append("")
        lines.append("| tool | Method ID | doc |")
        lines.append("|------|-----------|-----|")
        for r in rows:
            link = _relink(TOOL_CATALOG_MD, r["docPath"])
            lines.append("| `%s` | `%s` | [doc](%s) |"
                         % (r["tool"], r["methodId"], link))
        lines.append("")

    lines.append(END_MARK)
    lines.append("")
    return "\n".join(lines)


def render_domain_slice_md(domain, server_dirs, records, target_file):
    wanted = set(server_dirs)
    rows = [r for r in records if r["_serverDir"] in wanted]
    # Already globally sorted by (_serverDir, tool).
    lines = []
    lines.append("# `%s` — tool slice (generated)" % domain)
    lines.append("")
    lines.append("> Generated by `scripts/skills/gen-catalog.py` from "
                 "`docs/mcp/tools/**` via `docs/skills/domain-map.json`. "
                 "Covers servers: %s. Do not edit inside the markers by hand."
                 % ", ".join(server_dirs))
    lines.append("")
    lines.append(BEGIN_MARK)
    lines.append("")
    lines.append("| tool | server | Method ID | doc |")
    lines.append("|------|--------|-----------|-----|")
    for r in rows:
        link = _relink(target_file, r["docPath"])
        lines.append("| `%s` | %s | `%s` | [doc](%s) |"
                     % (r["tool"], r["server"], r["methodId"], link))
    lines.append("")
    lines.append(END_MARK)
    lines.append("")
    return "\n".join(lines)

# ---------------------------------------------------------------------------
# Planning: (target_path -> content) for every file we own
# ---------------------------------------------------------------------------

def build_plan(records):
    """Return (plan, skipped_domains).

    plan: dict of absolute target path -> desired content string.
    skipped_domains: list of (domain, reason) for domains whose skill dir
                     does not exist yet.
    """
    plan = {
        CATALOG_JSON: render_catalog_json(records),
        TOOL_CATALOG_MD: render_tool_catalog_md(records),
    }

    skipped = []
    domain_map = load_domain_map()
    for domain in sorted(domain_map):
        if domain == "_comment":
            continue
        server_dirs = domain_map[domain]
        skill_dir = os.path.join(CLAUDE_SKILLS_ROOT, domain)
        if not os.path.isdir(skill_dir):
            skipped.append((domain, "skill dir absent: %s"
                            % os.path.relpath(skill_dir, REPO_ROOT)))
            continue
        target = os.path.join(skill_dir, "_generated", "tools.md")
        plan[target] = render_domain_slice_md(domain, server_dirs,
                                              records, target)
    return plan, skipped


def load_domain_map():
    if not os.path.isfile(DOMAIN_MAP_PATH):
        die("domain map not found: %s" % DOMAIN_MAP_PATH)
    with open(DOMAIN_MAP_PATH, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        die("domain-map.json must be a JSON object")
    return data

# ---------------------------------------------------------------------------
# Write / check
# ---------------------------------------------------------------------------

def write_plan(plan, skipped, total):
    for path, content in sorted(plan.items()):
        os.makedirs(os.path.dirname(path), exist_ok=True)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(content)
        info("wrote %s" % os.path.relpath(path, REPO_ROOT))
    info("catalog: %d tools across %d servers"
         % (total, len(EXPECTED_COUNTS)))
    for domain, reason in skipped:
        info("skipped domain '%s' (%s)" % (domain, reason))


def check_plan(plan, skipped, total):
    """Compare desired content against committed files. Return exit code."""
    stale = []
    missing = []

    # Write the plan to a temp location so --check never touches committed files.
    tmp = tempfile.mkdtemp(prefix="gen-catalog-check-")
    try:
        for path, content in plan.items():
            if not os.path.isfile(path):
                missing.append(path)
                continue
            with open(path, "r", encoding="utf-8") as fh:
                current = fh.read()
            if current != content:
                stale.append(path)
    finally:
        # Nothing persistent was written into tmp; just remove the empty dir.
        try:
            os.rmdir(tmp)
        except OSError:
            pass

    if not stale and not missing:
        info("catalog fresh: %d tools; %d generated file(s) up to date"
             % (total, len(plan)))
        for domain, reason in skipped:
            info("skipped domain '%s' (%s)" % (domain, reason))
        return 0

    warn("CATALOG DRIFT DETECTED -- run: python3 scripts/skills/gen-catalog.py")
    for path in sorted(missing):
        warn("  MISSING:  %s" % os.path.relpath(path, REPO_ROOT))
    for path in sorted(stale):
        warn("  STALE:    %s" % os.path.relpath(path, REPO_ROOT))
    return 1

# ---------------------------------------------------------------------------
# Small IO helpers
# ---------------------------------------------------------------------------

def info(msg):
    print("[gen-catalog] " + msg)


def warn(msg):
    print("[gen-catalog] WARN: " + msg, file=sys.stderr)


def die(msg):
    print("[gen-catalog] ERROR: " + msg, file=sys.stderr)
    sys.exit(2)

# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv):
    parser = argparse.ArgumentParser(
        prog="gen-catalog.py",
        description="Anti-drift catalog generator for the Seqeron skills "
                    "system. Reads docs/mcp/tools/** and emits the generated "
                    "tool catalog (+ per-domain slices).",
        epilog="Default mode writes files; --check verifies freshness "
               "(exit 1 on drift) without modifying anything.",
    )
    parser.add_argument(
        "--check", action="store_true",
        help="Verify committed generated files are fresh; exit 1 on drift. "
             "Does not modify any committed file.",
    )
    args = parser.parse_args(argv)

    records, per_server = collect_tools()
    total = crosscheck_counts(per_server)
    plan, skipped = build_plan(records)

    if args.check:
        return check_plan(plan, skipped, total)
    write_plan(plan, skipped, total)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
