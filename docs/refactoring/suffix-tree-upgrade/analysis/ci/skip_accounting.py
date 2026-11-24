#!/usr/bin/env python3
"""STUB of the tests-first commit of work package S2.2: it accounts nothing (every run holds), so that test_skip_accounting.py is red against it."""
import sys
from collections import namedtuple

Row = namedtuple("Row", "kind test in_builds reason how")
Report = namedtuple("Report", "failures counts proposals lines")


class ListError(Exception):
    pass


def load_list(path):
    return []


def account(results, rows, mode, partial=False):
    return Report([], {"skipped": 0, "explicit": 0, "executed": 0}, [], [])


def main(argv):
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
