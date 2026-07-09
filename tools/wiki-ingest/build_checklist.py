"""(Re)build WIKI_INGEST_CHECKLIST.md from every file under docs/.

Run from the repo root. Refuses to overwrite an existing checklist unless
--force is passed, so a resume in a fresh context never wipes progress.
"""
import os
import sys
import pathlib
from collections import defaultdict

CHK = pathlib.Path("WIKI_INGEST_CHECKLIST.md")


def main() -> None:
    if CHK.exists() and "--force" not in sys.argv:
        print("EXISTS: WIKI_INGEST_CHECKLIST.md already present; pass --force to rebuild "
              "(this DISCARDS all [x] progress).")
        sys.exit(1)

    root = pathlib.Path("docs")
    files = sorted(
        (pathlib.Path(dp) / fn).as_posix()
        for dp, _, fns in os.walk(root)
        for fn in fns
    )
    total = len(files)

    groups = defaultdict(list)
    for f in files:
        parts = f.split("/")
        groups[parts[1] if len(parts) > 2 else "(top-level)"].append(f)

    lines = [
        "# Wiki Ingest Checklist",
        "",
        f"Checklist of all {total} files under `docs/` to ingest into the LLM Wiki, "
        "one per `/wiki:ingest` run.",
        "Each processed file is marked `[x]` and committed. Do not stop until all are done.",
        "",
        f"Progress: 0 / {total}",
        "",
    ]
    for key in sorted(groups):
        lines.append(f"## {key}  ({len(groups[key])})")
        lines.append("")
        lines.extend(f"- [ ] {f}" for f in groups[key])
        lines.append("")

    CHK.write_text("\n".join(lines), encoding="utf-8")
    print(f"WROTE {total} files")


if __name__ == "__main__":
    main()
