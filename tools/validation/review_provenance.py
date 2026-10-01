#!/usr/bin/env python3
"""Regenerate docs/Validation/review-2026-09/REVIEW_PROVENANCE.md.

Purpose: every unit of the 2026-09 review campaign was reviewed by an autonomous
agent session, not by a human. When the review tooling/platform is reworked, those
units must be re-runnable *selectively*. This script derives, per unit:
  * the batch that owns it and the verdict its batch report records,
  * the commits on the campaign branch whose subject names the unit,
  * the agent session that produced them (Claude-Session trailer of each commit).
The reviewing agent itself is recorded per commit in the Co-Authored-By trailer,
so `git log --format='%h %s %(trailers:key=Co-Authored-By,valueonly)'` recovers it.

Usage:  python3 tools/validation/review_provenance.py [--since <rev>] > out.md
"""
import argparse, collections, re, subprocess, sys

CAMPAIGN = "review-2026-09"
README = f"docs/Validation/{CAMPAIGN}/README.md"


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=True).stdout


def batch_units() -> list[tuple[str, str]]:
    text = open(README).read()
    rows = re.findall(r"^\| (B\d\d) \| (.+?) \| (.+?) \|$", text, re.M)
    return [(b, u) for b, units, _ in rows for u in units.split()]


def report_verdicts(batch: str) -> dict[str, tuple[str, str, str]]:
    """unit -> (Stage A, Stage B, State) from the batch report's unit table."""
    try:
        text = open(f"docs/Validation/{CAMPAIGN}/{batch}.md").read()
    except FileNotFoundError:
        return {}
    out = {}
    for line in text.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 4 and re.fullmatch(r"[A-Z0-9-]+-\d{3}", cells[0]):
            out[cells[0]] = (cells[1], cells[2], cells[3])
    return out


def commits(since: str) -> list[tuple[str, str, str]]:
    fmt = "%h\x01%s\x01%(trailers:key=Claude-Session,valueonly)"
    log = git("log", "--format=" + fmt, f"{since}..HEAD")
    return [tuple(l.split("\x01")) for l in log.strip().split("\n") if l]  # type: ignore[misc]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--since", default="", help="rev to start from (default: campaign start, auto-detected)")
    args = ap.parse_args()
    since = args.since or git("log", "--format=%h", "-1", "--diff-filter=A", "--", README).strip() + "~1"

    units = batch_units()
    log = commits(since)
    by_unit: dict[str, list[tuple[str, str, str]]] = collections.defaultdict(list)
    for sha, subj, sess in log:
        for _, unit in units:
            if unit in subj:
                by_unit[unit].append((sha, subj, sess.strip()))

    print(f"# Review provenance & re-review queue — campaign {CAMPAIGN}\n")
    print("> **Generated** by `tools/validation/review_provenance.py` — do not hand-edit.\n>")
    print("> Every unit below was reviewed by an autonomous agent session (Stage A + Stage B per")
    print("> `docs/Validation/VALIDATION_PROTOCOL.md`), not by a human reviewer. This table exists so the")
    print("> whole campaign can be **re-run selectively** once the review tooling is reworked: take the")
    print("> units whose `Re-review` box is unchecked (or all of them, per batch) and replay the protocol.")
    print("> The reviewing agent and its session are recorded *per commit* on the campaign branch:")
    print("> `git log --format='%h %s %(trailers:key=Co-Authored-By,valueonly) %(trailers:key=Claude-Session,valueonly)'`.\n")
    print(f"Campaign start: `{since}`. Commits scanned: {len(log)}. Units: {len(units)}.\n")
    print("| Unit | Batch | Stage A | Stage B | State | Review commits | Re-review |")
    print("|---|---|---|---|---|---|---|")
    pending_batches: dict[str, int] = collections.Counter()
    for batch, unit in units:
        a, b, st = report_verdicts(batch).get(unit, ("—", "—", "not reviewed"))
        shas = " ".join(f"`{s}`" for s, _, _ in by_unit.get(unit, [])[:6]) or "—"
        if st == "not reviewed":
            pending_batches[batch] += 1
        print(f"| {unit} | {batch} | {a} | {b} | {st} | {shas} | ☐ |")
    print("\n## Batches with units not yet recorded as reviewed\n")
    for batch, n in sorted(pending_batches.items()):
        print(f"- {batch}: {n} unit(s) without a verdict row in `{batch}.md`")
    return 0


if __name__ == "__main__":
    sys.exit(main())
