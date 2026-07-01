#!/usr/bin/env python3
"""Verify that every relative Markdown link in the skills + skill-docs trees resolves.

Part of the Seqeron skills anti-drift guardrail (see docs/skills/STRATEGY.md §6).
A prior template bug used ../../../ from reference/ files where ../../../../ is
required; this check makes that class of error fail fast in CI and locally.

Scans:  .claude/skills/**/*.md  and  docs/skills/**/*.md
Checks:  every `](<relative>)` target that starts with `./` or `../` exists on disk
         (anchors after `#` are ignored). Absolute URLs and in-page #anchors are skipped.

Exit 0 = all links resolve. Exit 1 = one or more broken links (each printed).
Run from the repo root.
"""
import glob
import os
import re
import sys

LINK_RE = re.compile(r'\]\((\.{1,2}/[^)]+)\)')
TREES = (".claude/skills", ".github/skills", "docs/skills")


def main() -> int:
    broken = []
    checked = 0
    files = []
    for tree in TREES:
        files.extend(glob.glob(os.path.join(tree, "**", "*.md"), recursive=True))
    for f in sorted(files):
        base = os.path.dirname(f)
        with open(f, encoding="utf-8") as fh:
            text = fh.read()
        for m in LINK_RE.finditer(text):
            raw = m.group(1).split("#", 1)[0].strip()
            if not raw:
                continue
            checked += 1
            target = os.path.normpath(os.path.join(base, raw))
            if not os.path.exists(target):
                broken.append((f, m.group(1)))
    if broken:
        print(f"[check-links] BROKEN relative links: {len(broken)} (of {checked} checked)")
        for f, link in broken:
            print(f"  {f}  ->  {link}")
        return 1
    print(f"[check-links] OK: all {checked} relative links across "
          f"{len(files)} files resolve")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
