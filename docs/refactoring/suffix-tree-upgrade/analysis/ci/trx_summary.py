#!/usr/bin/env python3
"""The job summary of a test run for the deep workflow (work package S2.2): reads a .trx and prints Markdown (the step appends it to $GITHUB_STEP_SUMMARY).

usage: trx_summary.py [--title TEXT] [--expected-red FILE] [--grep PREFIX] FILE.trx
  --expected-red FILE  the tests that are red by design (large_expected_red.txt: `test<TAB>package`): the red tests are split into the expected ones, the UNEXPECTED ones
                       (a failure that is not in the list: a real finding) and the ones of the list that are green now (the row is to be removed)
  --grep PREFIX        also print the lines of the tests' output that start with PREFIX (the LARGE-REPORT lines of the Large category)
Always exits 0: the verdict of the run is the exit code of the test step; a trx that is missing (the build failed) is said in the summary and is not an error here.
Only the Python standard library.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trx_lib  # noqa: E402

LIMIT = 60


def read_expected(path):
    expected = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            cells = line.rstrip("\n").split("\t")
            expected[cells[0]] = cells[1] if len(cells) > 1 else ""
    return expected


def first_line(text, limit=200):
    text = " ".join((text or "").split())
    return text if len(text) <= limit else text[:limit - 3] + "..."


def minutes(seconds):
    return f"{int(seconds // 60)} min {int(seconds % 60):02d} s"


def build(results, wall, title, expected=None, grep=None):
    from collections import Counter
    counts = Counter(r.outcome for r in results)
    failed = [r for r in results if r.outcome in ("Failed", "Error", "Timeout", "Aborted")]
    lines = [f"### {title}", ""]
    ran = sum(v for k, v in counts.items() if k in trx_lib.EXECUTED)
    lines.append(f"{counts['Passed']} passed, {len(failed)} failed, {sum(v for k, v in counts.items() if k not in trx_lib.EXECUTED)} not executed; "
                 f"{ran} ran in {minutes(wall)} (the span of the run, from the trx).")
    lines.append("")
    if expected is not None:
        red = {r.key: r for r in failed}
        unexpected = [r for r in failed if r.key not in expected]
        as_expected = [r for r in failed if r.key in expected]
        green = [k for k in expected if k not in red and any(r.key == k for r in results)]
        absent = [k for k in expected if not any(r.key == k for r in results)]
        lines.append(f"**Red as expected ({len(as_expected)} of {len(expected)} listed):**")
        lines.extend(f"- `{r.key}` (green at {expected[r.key]}): {first_line(r.message, 140)}" for r in as_expected)
        lines.append("")
        lines.append(f"**UNEXPECTED failures ({len(unexpected)}):** " + ("none" if not unexpected else ""))
        lines.extend(f"- `{r.key}`: {first_line(r.message)}" for r in unexpected[:LIMIT])
        if green:
            lines.append("")
            lines.append(f"**Listed as red but green now ({len(green)}): remove the row from large_expected_red.txt:**")
            lines.extend(f"- `{k}`" for k in green)
        if absent:
            lines.append("")
            lines.append(f"**Listed but not in the run ({len(absent)}):** " + ", ".join(f"`{k}`" for k in absent))
        lines.append("")
    elif failed:
        lines.append(f"**Failed ({len(failed)}):**")
        lines.extend(f"- `{r.key}`: {first_line(r.message)}" for r in failed[:LIMIT])
        if len(failed) > LIMIT:
            lines.append(f"- ... and {len(failed) - LIMIT} more (the .trx artifact has them all)")
        lines.append("")
    if grep:
        found = [ln.strip() for r in results for ln in r.stdout.splitlines() if ln.strip().startswith(grep)]
        if found:
            lines.append(f"`{grep}` lines of the tests:")
            lines.append("")
            lines.append("```")
            lines.extend(found[:LIMIT])
            lines.append("```")
            lines.append("")
    return "\n".join(lines)


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument("trx")
    parser.add_argument("--title", default="Test run")
    parser.add_argument("--expected-red")
    parser.add_argument("--grep")
    args = parser.parse_args(argv)
    try:
        results = trx_lib.read_trx(args.trx)
    except trx_lib.TrxError as e:
        print(f"### {args.title}\n\nno usable test result file: {e}\n")
        return 0
    expected = read_expected(args.expected_red) if args.expected_red else None
    print(build(results, trx_lib.wall_seconds(args.trx), args.title, expected, args.grep))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
