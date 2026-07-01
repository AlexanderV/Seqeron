#!/usr/bin/env python3
"""Mirror the Seqeron skills from `.claude/skills/` into `.github/skills/` byte-identically.

The skill-development pipeline (`.github/agents/`) installs every skill into BOTH trees so it
works in Claude Code AND GitHub Copilot / VS Code; the coherence auditor checks they stay
byte-identical. This script is the mechanism + the CI guard for that invariant.

Scope: only the Seqeron-authored skills (dir name starts with `bio-` or `seqeron-`). Third-party
vendored skills (clean-code, clean-architecture) are Claude-only and are NOT mirrored.

Modes:
  (default)  copy .claude/skills/<seqeron-skill>/ -> .github/skills/<same>/, making the mirror
             an exact copy (adds/updates/removes files as needed).
  --check    verify the mirror already matches; exit 1 (listing drift) if not. Changes nothing.

Run from the repo root.
"""
import filecmp
import os
import shutil
import sys

SRC_ROOT = ".claude/skills"
DST_ROOT = ".github/skills"
PREFIXES = ("bio-", "seqeron-")


def seqeron_skills():
    return sorted(
        d for d in os.listdir(SRC_ROOT)
        if os.path.isdir(os.path.join(SRC_ROOT, d)) and d.startswith(PREFIXES)
    )


def rel_files(root):
    out = set()
    for dirpath, _dirs, files in os.walk(root):
        for fn in files:
            out.add(os.path.relpath(os.path.join(dirpath, fn), root))
    return out


def diff_skill(name):
    """Return list of drift strings for one skill (empty = in sync)."""
    src, dst = os.path.join(SRC_ROOT, name), os.path.join(DST_ROOT, name)
    drift = []
    if not os.path.isdir(dst):
        return [f"{name}: missing in {DST_ROOT}"]
    sfiles, dfiles = rel_files(src), rel_files(dst)
    for f in sorted(sfiles - dfiles):
        drift.append(f"{name}: missing in mirror -> {f}")
    for f in sorted(dfiles - sfiles):
        drift.append(f"{name}: stale in mirror (extra) -> {f}")
    for f in sorted(sfiles & dfiles):
        if not filecmp.cmp(os.path.join(src, f), os.path.join(dst, f), shallow=False):
            drift.append(f"{name}: content differs -> {f}")
    return drift


def main():
    check = "--check" in sys.argv[1:]
    skills = seqeron_skills()
    if check:
        drift = []
        for name in skills:
            drift.extend(diff_skill(name))
        # extra skill dirs in the mirror that no longer exist in source
        if os.path.isdir(DST_ROOT):
            for d in sorted(os.listdir(DST_ROOT)):
                if d.startswith(PREFIXES) and d not in skills:
                    drift.append(f"{d}: present in {DST_ROOT} but not in {SRC_ROOT}")
        if drift:
            print(f"[sync-github-skills] MIRROR DRIFT: {len(drift)} issue(s)")
            for d in drift:
                print("  " + d)
            print("  fix: run  python3 scripts/skills/sync-github-skills.py")
            return 1
        print(f"[sync-github-skills] OK: {len(skills)} skills byte-identical in both trees")
        return 0

    os.makedirs(DST_ROOT, exist_ok=True)
    for name in skills:
        src, dst = os.path.join(SRC_ROOT, name), os.path.join(DST_ROOT, name)
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
        print(f"[sync-github-skills] mirrored {name}")
    print(f"[sync-github-skills] {len(skills)} skills synced -> {DST_ROOT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
