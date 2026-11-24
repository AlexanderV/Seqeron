#!/usr/bin/env python3
"""The skip accounting of work package S2.2 (review finding R04: "a skipped test is green").

A test cannot see the outcome of the other tests of its own run, so the check reads the .trx of a finished run and compares it with a checked-in list,
tests/SuffixTree/SuffixTree.Esa.Tests/skips.tsv (the columns are described in its header). It fails (exit 1) on

  * UNLISTED SKIP               a test that was not executed (NotExecuted: an Assert.Ignore, an Assume, an Explicit test that was not selected; Inconclusive, ...)
                                and is not in the list: `Assert.Ignore("x")` added to any test is this failure;
  * STALE                       a row that says "not executed" whose test ran (the list cannot rot: a test that is no longer skipped must leave the list);
  * CHECKED-BUILD TEST SKIPPED  (mode checked) one of the tests that exist only for the checked accessor was skipped: the build was not checked after all;
  * LISTED BUT ABSENT           a row whose test is not in the run at all (a renamed or deleted test; a filtered run; not reported with --partial).

--mode default   a Release build of the default set (what dotnet.yml runs): `ignore` rows of in_builds `default` and `both` are skipped, `host` rows are free.
--mode checked   the checked build (-p:EsaBounds=true): only `both` is skipped; `default` rows must have RUN; `host` rows are free.
Explicit tests (kind explicit, explicit-class) are never skips of the run: the console of `dotnet test` does not count them (the default run says 37 skipped while its trx
holds 897 not-executed results), and the list names them so that a test cannot be made Explicit, and leave the run, unnoticed. The trx cannot tell an ignored test from
an Explicit one (both are NotExecuted with a message), so the list is the only judge.

usage: skip_accounting.py --trx FILE --list FILE --mode {default,checked} [--partial] [--propose]
exit 0: the run is accounted for; 1: at least one failure (each is printed); 2: the list or the trx cannot be read (or a usage error).
--partial  the run is a filtered one (a hand run of one class): rows without a result are not reported (everything else is).
--propose  after the failures, print the rows of the unlisted skips in the format of the list (to be completed with a real reason).
Only the Python standard library.
"""
import argparse
import os
import sys
from collections import Counter, defaultdict, namedtuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import trx_lib  # noqa: E402

HEADER = ["kind", "test", "in_builds", "reason", "how_to_unskip"]
KINDS = ("ignore", "explicit", "explicit-class")
IN_BUILDS = ("default", "both", "host")
MODES = ("default", "checked")

Row = namedtuple("Row", "kind test in_builds reason how")
Report = namedtuple("Report", "failures counts proposals lines")


class ListError(Exception):
    """The list is not well formed (the message names the line)."""


def load_list(path):
    try:
        with open(path, encoding="utf-8") as f:
            lines = f.read().split("\n")
    except OSError as e:
        raise ListError(f"cannot read {path}: {e}") from e
    rows, seen, header_seen, problems = [], {}, False, []
    for number, line in enumerate(lines, 1):
        if not line.strip() or line.startswith("#"):
            continue
        cells = line.split("\t")
        if not header_seen:
            if cells != HEADER:
                raise ListError(f"{path}:{number}: the header must be exactly {'<TAB>'.join(HEADER)}")
            header_seen = True
            continue
        if len(cells) != len(HEADER):
            problems.append(f"{path}:{number}: {len(cells)} columns, expected {len(HEADER)}")
            continue
        row = Row(*cells)
        if row.kind not in KINDS:
            problems.append(f"{path}:{number}: kind {row.kind!r} is not one of {KINDS}")
        if row.in_builds not in IN_BUILDS:
            problems.append(f"{path}:{number}: in_builds {row.in_builds!r} is not one of {IN_BUILDS}")
        if row.kind != "ignore" and row.in_builds != "both":
            problems.append(f"{path}:{number}: an {row.kind} row says in_builds `both` (it is not run in either build)")
        for name, value in zip(HEADER, row):
            if not value.strip() or value != value.strip():
                problems.append(f"{path}:{number}: column {name} is empty or has blanks around it")
        if row.test in seen:
            problems.append(f"{path}:{number}: duplicate test {row.test!r} (first at line {seen[row.test]})")
        seen[row.test] = number
        rows.append(row)
    if not header_seen:
        raise ListError(f"{path}: no header and no rows")
    if problems:
        raise ListError("\n".join(problems))
    return rows


def _short(text, limit=160):
    text = " ".join((text or "").split())
    return text if len(text) <= limit else text[:limit - 3] + "..."


def account(results, rows, mode, partial=False):
    """Compares the results of a run with the list. Returns a Report (failures, counts, proposals, lines of the summary)."""
    if mode not in MODES:
        raise ValueError(f"mode {mode!r} is not one of {MODES}")
    by_key = defaultdict(list)
    by_class = defaultdict(list)
    for r in results:
        by_key[r.key].append(r)
        by_class[r.class_name].append(r)
    ignore = {r.test: r for r in rows if r.kind == "ignore"}
    explicit = {r.test: r for r in rows if r.kind == "explicit"}
    explicit_class = {r.test: r for r in rows if r.kind == "explicit-class"}

    failures, proposals, proposed_classes = [], [], set()
    counts = Counter()
    for r in results:
        if r.outcome in trx_lib.EXECUTED:
            counts["executed"] += 1
            continue
        # a result that did not run: it is a skip, an Explicit test, or an unlisted skip
        if r.key in explicit or r.class_name in explicit_class:
            counts["explicit"] += 1
        elif r.key in ignore:
            row = ignore[r.key]
            counts["skipped"] += 1
            if mode == "checked" and row.in_builds == "default":
                failures.append(f"CHECKED-BUILD TEST SKIPPED {r.key}: {_short(r.message)} (a test that needs the checked accessor was skipped: this run is not a checked build)")
        else:
            counts["unlisted"] += 1
            failures.append(f"UNLISTED SKIP {r.key} [{r.outcome}]: {_short(r.message)}")
            if r.message.startswith("OneTimeSetUp:"):
                if r.class_name not in proposed_classes:
                    proposed_classes.add(r.class_name)
                    proposals.append("\t".join(["explicit-class", r.class_name, "both", _short(r.message[len("OneTimeSetUp:"):]), "TODO"]))
            else:
                in_builds = "default" if "checked build" in r.message or "EsaBounds" in r.message else "both"
                proposals.append("\t".join(["ignore", r.key, in_builds, _short(r.message), "TODO"]))

    # the rows: every one must be what it says
    for row in rows:
        if row.kind == "explicit-class":
            group = by_class.get(row.test, [])
            where = f"the fixture {row.test}"
        else:
            group = by_key.get(row.test, [])
            where = f"{row.test}"
        if not group:
            if not partial:
                failures.append(f"LISTED BUT ABSENT {where}: the run has no such test (renamed, deleted or a filtered run; --partial for a filtered run)")
            continue
        ran = [r for r in group if r.outcome in trx_lib.EXECUTED]
        if row.kind in ("explicit", "explicit-class"):
            if ran:
                failures.append(f"STALE {where}: listed as {row.kind} (never run) but {len(ran)} of {len(group)} result(s) ran: remove the row")
        elif row.in_builds == "both" or (row.in_builds == "default" and mode == "default"):
            if ran:
                failures.append(f"STALE {where}: listed as skipped in a {mode} run but {len(ran)} of {len(group)} result(s) ran ({ran[0].outcome}): remove the row")
        # a `default` row in a checked run is judged in the loop above (skipped = CHECKED-BUILD TEST SKIPPED); a `host` row is free in both modes

    lines = [f"results in the run: {len(results)}",
             f"executed: {counts['executed']}",
             f"skipped (ignore rows of the list): {counts['skipped']}",
             f"explicit, not selected: {counts['explicit']}",
             f"unlisted skips: {counts['unlisted']}"]
    return Report(failures, counts, proposals, lines)


def main(argv):
    parser = argparse.ArgumentParser(description="skip accounting of a .trx against skips.tsv", add_help=True)
    parser.add_argument("--trx", required=True)
    parser.add_argument("--list", required=True, dest="list_path")
    parser.add_argument("--mode", required=True)
    parser.add_argument("--partial", action="store_true")
    parser.add_argument("--propose", action="store_true")
    try:
        args = parser.parse_args(argv)
    except SystemExit as e:
        return 2 if e.code else 0
    if args.mode not in MODES:
        print(f"--mode must be one of {MODES}, not {args.mode!r}", file=sys.stderr)
        return 2
    try:
        rows = load_list(args.list_path)
        results = trx_lib.read_trx(args.trx)
    except (ListError, trx_lib.TrxError) as e:
        print(f"CANNOT ACCOUNT: {e}", file=sys.stderr)
        return 2
    report = account(results, rows, args.mode, args.partial)
    print(f"### Skip accounting ({args.mode} build{', partial run' if args.partial else ''})")
    print()
    for line in report.lines:
        print(f"- {line}")
    print()
    if report.failures:
        for failure in report.failures:
            print(f"FAIL: {failure}")
        if args.propose and report.proposals:
            print()
            print("rows to add (complete the reason and how_to_unskip; add a row only for a skip that is meant to stay):")
            for proposal in report.proposals:
                print(proposal)
        print()
        print(f"{len(report.failures)} discrepancy(ies) between the run and {os.path.basename(args.list_path)}")
        return 1
    print(f"the skips of the run are exactly the skips of the list: holds ({report.counts['skipped']} skipped, {report.counts['explicit']} explicit not selected)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
