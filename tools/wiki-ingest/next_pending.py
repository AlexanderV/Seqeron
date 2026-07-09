"""Print the next unprocessed file from WIKI_INGEST_CHECKLIST.md.

Run from the repo root. Prints a progress comment line, then either the first
`- [ ] <path>` entry or `ALL_DONE`.
"""
import re
import pathlib

CHK = pathlib.Path("WIKI_INGEST_CHECKLIST.md")


def main() -> None:
    t = CHK.read_text(encoding="utf-8")
    pending = re.findall(r"^- \[ \] (.+)$", t, flags=re.M)
    done = len(re.findall(r"^- \[x\] ", t, flags=re.M))
    print(f"# {done}/{done + len(pending)} done, {len(pending)} pending")
    print(pending[0] if pending else "ALL_DONE")


if __name__ == "__main__":
    main()
