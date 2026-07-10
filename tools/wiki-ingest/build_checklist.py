"""(Re)build WIKI_INGEST_CHECKLIST.md from human-readable docs under docs/.

Only human-information documents are listed (prose/spec files a person reads):
`.md`, `.markdown`, `.pdf`, `.txt`, `.rst`, `.adoc`. Machine-readable artifacts
— MCP tool schemas and other `.json`, configs, metadata, source files — are
deliberately excluded; they are not sources to summarize into the wiki.

Run from the repo root. Refuses to overwrite an existing checklist unless
--force is passed. On --force it PRESERVES existing `[x]` progress: any listed
file already marked done in the current checklist stays done.
"""
import os
import sys
import pathlib
import re
from collections import defaultdict

CHK = pathlib.Path("WIKI_INGEST_CHECKLIST.md")

# Human-readable document extensions only. Everything else (.json, .yaml, .toml,
# .cs, configs, metadata) is intentionally excluded.
DOC_EXTS = {".md", ".markdown", ".pdf", ".txt", ".rst", ".adoc"}


def read_done(path: pathlib.Path) -> dict[str, str]:
    """Map each done file path -> its trailing annotation (e.g. '  (done-via-concept)').

    Preserves both the [x] state and any annotation, keyed by the bare docs path so
    it matches the filesystem enumeration.
    """
    if not path.exists():
        return {}
    done: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        m = re.match(r"- \[x\] (\S.*?)(\s{2,}\([^)]*\))?\s*$", line.strip())
        if m:
            done[m.group(1)] = m.group(2) or ""
    return done


def main() -> None:
    if CHK.exists() and "--force" not in sys.argv:
        print("EXISTS: WIKI_INGEST_CHECKLIST.md already present; pass --force to rebuild "
              "([x] progress is preserved on rebuild).")
        sys.exit(1)

    done = read_done(CHK)

    root = pathlib.Path("docs")
    files = sorted(
        (pathlib.Path(dp) / fn).as_posix()
        for dp, _, fns in os.walk(root)
        for fn in fns
        if pathlib.Path(fn).suffix.lower() in DOC_EXTS
    )
    total = len(files)
    done_count = sum(1 for f in files if f in done)
    # Done paths from the old checklist that no longer exist on disk (should be none
    # once extensions match); surfaced so a rebuild never silently drops progress.
    orphaned = sorted(set(done) - set(files))

    groups = defaultdict(list)
    for f in files:
        parts = f.split("/")
        groups[parts[1] if len(parts) > 2 else "(top-level)"].append(f)

    lines = [
        "# Wiki Ingest Checklist",
        "",
        f"Checklist of all {total} human-readable documents under `docs/` "
        "(`.md`, `.pdf`, and other prose/spec files) to ingest into the LLM Wiki, "
        "one per `/wiki:ingest` run.",
        "Machine-readable artifacts (MCP `.json` tool schemas, configs, metadata, "
        "sources) are intentionally excluded — they are not human-information sources.",
        "Each processed file is marked `[x]` and committed. Do not stop until all are done.",
        "",
        f"Progress: {done_count} / {total}",
        "",
    ]
    for key in sorted(groups):
        lines.append(f"## {key}  ({len(groups[key])})")
        lines.append("")
        lines.extend(
            f"- [{'x' if f in done else ' '}] {f}{done.get(f, '')}" for f in groups[key]
        )
        lines.append("")

    CHK.write_text("\n".join(lines), encoding="utf-8")
    print(f"WROTE {total} documents ({done_count} already done, preserved)")
    if orphaned:
        print(f"WARNING: {len(orphaned)} done path(s) in the old checklist are not on disk "
              "and were dropped:")
        for f in orphaned:
            print(f"  - {f}")


if __name__ == "__main__":
    main()
