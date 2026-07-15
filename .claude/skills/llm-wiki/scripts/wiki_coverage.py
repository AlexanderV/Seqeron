#!/usr/bin/env python3
"""Report source files not covered by any wiki page.

Usage: python wiki_coverage.py <wiki-dir> [--docs-dir docs] [--ext .md,.rst,.txt]
       [--exclude glob ...] [--extra FILE ...] [--strict]

A source file is covered when at least one wiki page lists it in `sources:`
or `doc_path:` frontmatter. `--docs-dir` is scanned recursively; `--extra`
names individual source files that live outside the docs directory (e.g.
root-level README.md) — repeatable, and prints a warning (not an error) when
a named path does not exist, guarding against typos and renames. Default exit
code is 0 (report only); --strict exits 1 when uncovered files exist (CI gate).
Stdlib only.
"""
import argparse, fnmatch, re, sys
from pathlib import Path

FM = re.compile(r"^---\n(.*?)\n---", re.S)

def referenced_paths(wiki_dir):
    refs = set()
    for page in Path(wiki_dir).rglob("*.md"):
        m = FM.match(page.read_text(encoding="utf-8", errors="replace"))
        if not m:
            continue
        key = None
        for line in m.group(1).splitlines():
            if re.match(r"^(sources|doc_path):", line):
                key, _, val = line.partition(":")
                val = val.strip().strip("\"'")
                if val:
                    refs.add(val)
                    key = None
            elif key and line.strip().startswith("- "):
                refs.add(line.strip()[2:].strip().strip("\"'"))
            elif key and line.strip() and not line.startswith(" "):
                key = None
    return refs

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wiki_dir")
    ap.add_argument("--docs-dir", default="docs")
    ap.add_argument("--ext", default=".md")
    ap.add_argument("--exclude", action="append", default=[])
    ap.add_argument("--extra", action="append", default=[],
                    help="Individual source file outside --docs-dir (repeatable).")
    ap.add_argument("--strict", action="store_true")
    args = ap.parse_args()
    exts = {e if e.startswith(".") else "." + e for e in args.ext.split(",")}
    refs = referenced_paths(args.wiki_dir)
    uncovered = []
    candidates = list(Path(args.docs_dir).rglob("*"))
    for extra in args.extra:
        p = Path(extra)
        if not p.is_file():
            print(f"WARNING: --extra path not found (typo or renamed?): {extra}",
                  file=sys.stderr)
            continue
        candidates.append(p)
    for f in sorted(candidates):
        if not f.is_file():
            continue
        # docs-dir files are filtered by extension; --extra files are named
        # explicitly, so accept them regardless of suffix.
        rel = f.as_posix()
        is_extra = rel in args.extra
        if not is_extra and f.suffix not in exts:
            continue
        if any(fnmatch.fnmatch(rel, pat) for pat in args.exclude):
            continue
        if rel not in refs and rel not in uncovered:
            uncovered.append(rel)
    if uncovered:
        print(f"UNCOVERED sources ({len(uncovered)}): no wiki page references them")
        for rel in uncovered:
            print(f"  {rel}")
        if args.strict:
            sys.exit(1)
    else:
        print("OK: all source files are referenced by the wiki")

if __name__ == "__main__":
    main()
