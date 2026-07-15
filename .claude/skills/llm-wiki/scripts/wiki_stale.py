#!/usr/bin/env python3
"""Flag wiki pages not refreshed after docs/** sources changed.

Usage: python wiki_stale.py <wiki-dir> [--repo-root .]
Exit code 1 if any stale pages found (usable as a CI/pre-commit gate).
Requires: git available, run inside the repository. Stdlib only.
"""
import argparse
import re
import subprocess
import sys
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


def commits_after(root: Path, baseline: str, path: str) -> tuple[list[str], str | None]:
    result = subprocess.run(
        ["git", "log", "--format=%H", f"{baseline}..HEAD", "--", path],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return [], result.stderr.strip() or "git log failed"
    return result.stdout.split(), None


def source_is_covered(
    root: Path,
    baseline: str,
    source_path: str,
    page_path: str,
) -> tuple[bool, int, str | None]:
    """Return whether the latest source change has a same-or-later page refresh."""
    source_commits, error = commits_after(root, baseline, source_path)
    if error or not source_commits:
        return not error, len(source_commits), error

    page_commits, page_error = commits_after(root, baseline, page_path)
    if page_error:
        return False, len(source_commits), page_error

    latest_source = source_commits[0]
    for page_commit in page_commits:
        covered = subprocess.run(
            ["git", "merge-base", "--is-ancestor", latest_source, page_commit],
            cwd=root,
            capture_output=True,
            text=True,
        )
        if covered.returncode == 0:
            return True, len(source_commits), None
        if covered.returncode not in {0, 1}:
            return False, len(source_commits), covered.stderr.strip() or "git merge-base failed"
    return False, len(source_commits), None


def staged_paths(root: Path) -> tuple[set[str], str | None]:
    result = subprocess.run(
        ["git", "diff", "--cached", "--name-only", "--diff-filter=ACMR", "-z"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return set(), result.stderr.strip() or "git diff --cached failed"
    return {path.replace("\\", "/") for path in result.stdout.split("\0") if path}, None


def source_path_is_staged(root: Path, source_path: str, staged: set[str]) -> bool:
    normalized = source_path.replace("\\", "/").rstrip("/")
    if (root / source_path).is_dir():
        prefix = normalized + "/"
        return any(path.startswith(prefix) for path in staged)
    return normalized in staged


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("wiki_dir")
    ap.add_argument("--repo-root", default=".")
    args = ap.parse_args()
    root = Path(args.repo_root).resolve()
    staged, staged_error = staged_paths(root)
    if staged_error:
        print(f"WARN staged-change check failed: {staged_error}", file=sys.stderr)
    stale = []
    for page in sorted(Path(args.wiki_dir).rglob("*.md")):
        fm = frontmatter(page.read_text(encoding="utf-8", errors="replace"))
        commit = fm.get("source_commit") or ""
        paths = fm.get("sources") or ([fm["doc_path"]] if fm.get("doc_path") else [])
        if not commit or not paths or isinstance(commit, list):
            continue
        try:
            page_path = page.resolve().relative_to(root).as_posix()
        except ValueError:
            print(f"WARN {page}: page is outside repository root {root}", file=sys.stderr)
            continue
        page_staged = page_path in staged
        for src in paths:
            if not (root / src).exists():
                stale.append((str(page), src, "MISSING", 0))
                continue
            source_staged = source_path_is_staged(root, src, staged)
            if source_staged and not page_staged:
                stale.append((str(page), src, "STAGED", 1))
                continue
            covered, count, error = source_is_covered(root, commit, src, page_path)
            if error:
                print(f"WARN {page}: history check failed for {src}: {error}", file=sys.stderr)
                continue
            if not covered and not page_staged:
                stale.append((str(page), src, commit[:9], count))
    if stale:
        print("STALE pages (source changed after source_commit):")
        for page, src, commit, n in stale:
            if commit == "MISSING":
                print(f"  {page}: source {src} no longer exists (deleted or renamed)")
            elif commit == "STAGED":
                print(f"  {page}: source {src} is staged without a staged page refresh")
            else:
                print(f"  {page}: {src} has {n} commit(s) since {commit}")
        sys.exit(1)
    print("OK: no stale pages")

if __name__ == "__main__":
    main()
