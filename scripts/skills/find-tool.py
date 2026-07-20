#!/usr/bin/env python3
"""find-tool.py — discovery search over Seqeron's MCP tool + algorithm docs.

Answers "is there a Seqeron tool/algorithm for X?" without loading all 427 MCP
tool schemas into context. Greps the per-tool docs under docs/mcp/tools/ (and,
with --algorithms, the algorithm docs under docs/algorithms/) for keywords and
prints a compact, greppable table pointing at the full doc + Method ID.

Python 3.9+ stdlib only. No third-party deps. Works purely by reading the docs;
if docs/skills/_generated/catalog.json exists it is used only as optional
enrichment (never required).

Usage:
    python3 scripts/skills/find-tool.py <keywords...> [--server <name>]
                                        [--algorithms] [--limit N]

Examples:
    python3 scripts/skills/find-tool.py melting temperature
    python3 scripts/skills/find-tool.py crispr --server moltools
    python3 scripts/skills/find-tool.py alignment --algorithms
"""
import argparse
import json
import os
import re
import sys

# Repo root = two levels up from this file (scripts/skills/find-tool.py).
HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))

TOOLS_DIR = os.path.join(REPO_ROOT, "docs", "mcp", "tools")
ALGOS_DIR = os.path.join(REPO_ROOT, "docs", "algorithms")
CATALOG = os.path.join(REPO_ROOT, "docs", "skills", "_generated", "catalog.json")
# LLM Wiki: concept/gotcha pages that curate the science + sharp edges behind a
# tool. The mapping is the wiki's own `mcp_tools:` frontmatter — read live so the
# WIKI column never drifts and lights up automatically as ingests add tools.
WIKI_DIR = os.path.join(REPO_ROOT, "wiki")

SERVERS = [
    "core", "sequence", "parsers", "alignment", "analysis", "annotation",
    "chromosome", "metagenomics", "moltools", "phylogenetics", "population",
]

# Match rows like: | **Server** | MolTools |  and | **Method ID** | `Foo.Bar` |
_SERVER_RE = re.compile(r"\|\s*\*\*Server\*\*\s*\|\s*([^|]+?)\s*\|")
_METHOD_RE = re.compile(r"\|\s*\*\*Method ID\*\*\s*\|\s*`([^`]+)`")
_ALGO_TEST_UNIT_RE = re.compile(r"\|\s*Test Unit ID\s*\|\s*([^|]+?)\s*\|")


def rel(path):
    return os.path.relpath(path, REPO_ROOT)


def parse_tool_doc(path):
    """Return dict with name, desc, server, method_id, path — or None."""
    try:
        with open(path, "r", encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    except (OSError, UnicodeDecodeError):
        return None

    name = ""
    desc = ""
    # H1 title = tool name; first non-blank line after it = one-line description.
    for i, line in enumerate(lines):
        if line.startswith("# "):
            name = line[2:].strip()
            for nxt in lines[i + 1:]:
                s = nxt.strip()
                if s:
                    desc = s
                    break
            break

    server = ""
    method_id = ""
    for line in lines:
        if not server:
            m = _SERVER_RE.search(line)
            if m:
                server = m.group(1).strip()
        if not method_id:
            m = _METHOD_RE.search(line)
            if m:
                method_id = m.group(1).strip()
        if server and method_id:
            break

    if not name:
        return None
    return {
        "name": name,
        "desc": desc,
        "server": server,
        "method_id": method_id,
        "path": rel(path),
    }


def parse_algo_doc(path):
    try:
        with open(path, "r", encoding="utf-8") as fh:
            lines = fh.read().splitlines()
    except (OSError, UnicodeDecodeError):
        return None

    title = ""
    for line in lines:
        if line.startswith("# "):
            title = line[2:].strip()
            break

    test_unit = ""
    for line in lines:
        m = _ALGO_TEST_UNIT_RE.search(line)
        if m:
            test_unit = m.group(1).strip()
            break

    if not title:
        return None
    return {"title": title, "test_unit": test_unit, "path": rel(path)}


def iter_tool_docs(server_filter):
    if not os.path.isdir(TOOLS_DIR):
        return
    dirs = sorted(os.listdir(TOOLS_DIR))
    for server_dir in dirs:
        if server_filter and server_dir != server_filter:
            continue
        full_dir = os.path.join(TOOLS_DIR, server_dir)
        if not os.path.isdir(full_dir):
            continue
        for fname in sorted(os.listdir(full_dir)):
            if not fname.endswith(".md") or fname == "README.md":
                continue
            yield os.path.join(full_dir, fname)


def iter_algo_docs():
    if not os.path.isdir(ALGOS_DIR):
        return
    for root, _dirs, files in os.walk(ALGOS_DIR):
        for fname in sorted(files):
            if not fname.endswith(".md") or fname in ("README.md", "CANONICAL_MAP.md"):
                continue
            yield os.path.join(root, fname)


def score_tool(doc, keywords):
    """All keywords must appear (AND). Rank by *where* they match.

    name/methodId hits weigh more than description hits. Returns (score, ok).
    """
    name_l = doc["name"].lower()
    method_l = doc["method_id"].lower()
    desc_l = doc["desc"].lower()
    haystack = " ".join((name_l, method_l, desc_l))

    total = 0
    for kw in keywords:
        if kw not in haystack:
            return (0, False)
        if kw in name_l:
            total += 100
        elif kw in method_l:
            total += 40
        elif kw in desc_l:
            total += 10
    return (total, True)


def score_algo(doc, keywords):
    title_l = doc["title"].lower()
    unit_l = doc["test_unit"].lower()
    path_l = doc["path"].lower()
    haystack = " ".join((title_l, unit_l, path_l))
    total = 0
    for kw in keywords:
        if kw not in haystack:
            return (0, False)
        if kw in title_l:
            total += 100
        elif kw in unit_l:
            total += 40
        else:
            total += 10
    return (total, True)


def load_catalog():
    """Optional enrichment. Returns dict keyed by method_id -> extra, or {}."""
    if not os.path.isfile(CATALOG):
        return {}
    try:
        with open(CATALOG, "r", encoding="utf-8") as fh:
            data = json.load(fh)
    except (OSError, ValueError):
        return {}
    index = {}
    entries = data.get("tools", data) if isinstance(data, dict) else data
    if isinstance(entries, list):
        for e in entries:
            if isinstance(e, dict):
                key = e.get("method_id") or e.get("name")
                if key:
                    index[key] = e
    return index


def build_tool_wiki_map():
    """tool name -> (wiki concept/gotcha slug, is_gotcha).

    Inverts the `mcp_tools:` frontmatter list of every wiki/concepts and
    wiki/gotchas page. This is the wiki's own curated tool->page binding, so the
    map stays in sync with the wiki with no generated artifact to refresh. A
    gotcha binding is flagged so the caller can surface it as a known-trap hint.
    """
    mapping = {}
    for sub, is_gotcha in (("gotchas", True), ("concepts", False)):
        d = os.path.join(WIKI_DIR, sub)
        if not os.path.isdir(d):
            continue
        for fname in sorted(os.listdir(d)):
            if not fname.endswith(".md"):
                continue
            slug = fname[:-3]
            try:
                with open(os.path.join(d, fname), "r", encoding="utf-8") as fh:
                    text = fh.read()
            except (OSError, UnicodeDecodeError):
                continue
            if not text.startswith("---"):
                continue
            end = text.find("\n---", 3)
            fm = text[: end if end != -1 else len(text)]
            tools = []
            inline = re.search(r"(?m)^mcp_tools:\s*\[([^\]]*)\]", fm)
            if inline:
                tools = re.findall(r"[A-Za-z0-9_]+", inline.group(1))
            else:
                block = re.search(r"(?m)^mcp_tools:\s*$", fm)
                if block:
                    for line in fm[block.end():].splitlines():
                        item = re.match(r"\s+-\s*(\S+)", line)
                        if item:
                            tools.append(item.group(1).strip())
                        elif line.strip() and not line[0].isspace():
                            break
            for t in tools:
                # Gotchas take precedence (traps matter more than the concept);
                # otherwise first page wins deterministically (gotchas scanned first).
                if t not in mapping or (is_gotcha and not mapping[t][1]):
                    mapping[t] = (slug, is_gotcha)
    return mapping


def print_tool_table(results, catalog, wiki_map):
    if not results:
        print("No matching tools.")
        return
    print("{:<33} {:<13} {:<40} {:<36} {}".format(
        "TOOL", "SERVER", "METHOD ID", "WIKI (concept / (!) gotcha)", "DOC"))
    print("-" * 140)
    for doc in results:
        extra = ""
        cat = catalog.get(doc["method_id"]) or catalog.get(doc["name"])
        if cat and cat.get("stability"):
            extra = "  [{}]".format(cat["stability"])
        wk = wiki_map.get(doc["name"])
        if wk:
            wiki_cell = ("(!) " if wk[1] else "") + wk[0]
        else:
            wiki_cell = "-"
        print("{:<33} {:<13} {:<40} {:<36} {}{}".format(
            doc["name"][:33],
            doc["server"][:13],
            (doc["method_id"] or "-")[:40],
            wiki_cell[:36],
            doc["path"],
            extra,
        ))


def print_algo_table(results):
    if not results:
        print("No matching algorithms.")
        return
    print("\nALGORITHMS")
    print("TITLE                                              TEST UNIT ID          DOC")
    print("-" * 120)
    for doc in results:
        print("{:<50} {:<21} {}".format(
            doc["title"][:50],
            (doc["test_unit"] or "-")[:21],
            doc["path"],
        ))


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="find-tool.py",
        description="Search Seqeron MCP tool docs (and algorithm docs) for keywords.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  python3 scripts/skills/find-tool.py melting temperature\n"
            "  python3 scripts/skills/find-tool.py crispr --server moltools\n"
            "  python3 scripts/skills/find-tool.py alignment --algorithms\n"
        ),
    )
    parser.add_argument("keywords", nargs="+", help="Keywords (ALL must match, case-insensitive).")
    parser.add_argument("--server", help="Filter to one server dir ({}).".format(", ".join(SERVERS)))
    parser.add_argument("--algorithms", action="store_true", help="Also search docs/algorithms/**.")
    parser.add_argument("--limit", type=int, default=20, help="Max results per section (default 20).")
    args = parser.parse_args(argv)

    keywords = [k.lower() for k in args.keywords]
    server_filter = args.server.lower() if args.server else None
    if server_filter and server_filter not in SERVERS:
        parser.error("unknown --server '{}'. Choose from: {}".format(server_filter, ", ".join(SERVERS)))

    catalog = load_catalog()
    wiki_map = build_tool_wiki_map()

    tool_hits = []
    for path in iter_tool_docs(server_filter):
        doc = parse_tool_doc(path)
        if not doc:
            continue
        s, ok = score_tool(doc, keywords)
        if ok:
            tool_hits.append((s, doc))
    # Deterministic: score desc, then name asc, then path asc.
    tool_hits.sort(key=lambda t: (-t[0], t[1]["name"], t[1]["path"]))
    tool_results = [d for _s, d in tool_hits[:args.limit]]

    print_tool_table(tool_results, catalog, wiki_map)
    if len(tool_hits) > args.limit:
        print("... {} more tool hits (use --limit to see more).".format(len(tool_hits) - args.limit))

    if args.algorithms:
        algo_hits = []
        for path in iter_algo_docs():
            doc = parse_algo_doc(path)
            if not doc:
                continue
            s, ok = score_algo(doc, keywords)
            if ok:
                algo_hits.append((s, doc))
        algo_hits.sort(key=lambda t: (-t[0], t[1]["title"], t[1]["path"]))
        algo_results = [d for _s, d in algo_hits[:args.limit]]
        print_algo_table(algo_results)
        if len(algo_hits) > args.limit:
            print("... {} more algorithm hits (use --limit to see more).".format(len(algo_hits) - args.limit))

    return 0


if __name__ == "__main__":
    sys.exit(main())
