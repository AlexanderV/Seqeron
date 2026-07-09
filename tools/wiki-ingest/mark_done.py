"""Mark one file as processed in WIKI_INGEST_CHECKLIST.md.

Usage: python tools/wiki-ingest/mark_done.py <docs-relative-path>
Flips `- [ ] <path>` to `- [x] <path>` and refreshes the Progress line.
Exit 2 if the path is not found as an unprocessed entry.
"""
import re
import sys
import pathlib

CHK = pathlib.Path("WIKI_INGEST_CHECKLIST.md")


def main() -> None:
    target = sys.argv[1].replace("\\", "/")
    t = CHK.read_text(encoding="utf-8")
    if f"- [ ] {target}" not in t:
        if f"- [x] {target}" in t:
            print(f"ALREADY-DONE {target}")
        else:
            print(f"NOT-FOUND {target}")
            sys.exit(2)
        return
    t = t.replace(f"- [ ] {target}", f"- [x] {target}", 1)
    done = len(re.findall(r"^- \[x\] ", t, flags=re.M))
    total = len(re.findall(r"^- \[[ x]\] ", t, flags=re.M))
    t = re.sub(r"^Progress: \d+ / \d+", f"Progress: {done} / {total}", t, flags=re.M)
    CHK.write_text(t, encoding="utf-8")
    print(f"MARKED {target} ({done}/{total})")


if __name__ == "__main__":
    main()
