#!/usr/bin/env python3
"""Flag wiki pages whose docs/** sources changed after source_commit.

Usage: python wiki_stale.py <wiki-dir> [--repo-root .]
Exit code 1 if any stale pages found (usable as a CI/pre-commit gate).
Requires: git available, run inside the repository. Stdlib only.
"""
import argparse, re, subprocess, sys
from pathlib import Path

FM = re.compile(r"^---\n(.*?)\n---", re.S)

def frontmatter(text):
    m = FM.match(text)
    if not m:
        return {}
    fm, key = {}, None
    for line in m.group(1).splitlines():
        if re.match(r"^[A-Za-z_][\w-]*:", line):
            key, _, val = line.partition(":")
            val = val.strip()
            fm[key.strip()] = val if val else []
        elif key and line.strip().startswith("- ") and isinstance(fm.get(key), list):
            fm[key].append(line.strip()[2:].strip().strip("\"'"))
    return fm

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wiki_dir")
    ap.add_argument("--repo-root", default=".")
    args = ap.parse_args()
    root = Path(args.repo_root).resolve()
    stale = []
    for page in sorted(Path(args.wiki_dir).rglob("*.md")):
        fm = frontmatter(page.read_text(encoding="utf-8", errors="replace"))
        commit = fm.get("source_commit") or ""
        paths = fm.get("sources") or ([fm["doc_path"]] if fm.get("doc_path") else [])
        if not commit or not paths or isinstance(commit, list):
            continue
        for src in paths:
            if not (root / src).exists():
                stale.append((str(page), src, "MISSING", 0))
                continue
            r = subprocess.run(
                ["git", "log", "--oneline", f"{commit}..HEAD", "--", src],
                cwd=root, capture_output=True, text=True)
            if r.returncode != 0:
                print(f"WARN {page}: git log failed for {src}: {r.stderr.strip()}", file=sys.stderr)
                continue
            if r.stdout.strip():
                n = len(r.stdout.strip().splitlines())
                stale.append((str(page), src, commit[:9], n))
    if stale:
        print("STALE pages (source changed after source_commit):")
        for page, src, commit, n in stale:
            if commit == "MISSING":
                print(f"  {page}: source {src} no longer exists (deleted or renamed)")
            else:
                print(f"  {page}: {src} has {n} commit(s) since {commit}")
        sys.exit(1)
    print("OK: no stale pages")

if __name__ == "__main__":
    main()
