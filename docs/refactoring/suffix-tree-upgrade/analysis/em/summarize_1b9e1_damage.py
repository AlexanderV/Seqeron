#!/usr/bin/env python3
"""Tables of run_1b9e1_damage.py from its rows: the orchestrator's cases before and after the fix (bytes that differ from the uninterrupted build), and the matrix (every marker x every kind of damage, each resumed in a NEW process) per size.
usage: summarize_1b9e1_damage.py <repro head rows> <repro fixed rows> <matrix rows> [<matrix rows> ...]"""
import collections
import sys


def read(path):
    return [l.rstrip("\n").split("\t") for l in open(path) if l.strip()]


def main():
    head, fixed = read(sys.argv[1]), read(sys.argv[2])
    print("| marker | damage | skipped (before / after) | bytes that differ from the uninterrupted build, before the fix | after the fix |")
    print("| --- | --- | --- | --- | --- |")
    for h, f in zip(head, fixed):
        print(f"| {h[1]} | {h[2]} | {h[4]} / {f[4]} | {h[6] if h[5] != 'eq' else 0} | {f[6] if f[5] != 'eq' else 0} |")
    print()
    print("| file | resumes (each in a new process) | markers | kinds of damage | exit 0 | SHA-256 equal | nothing but the target left |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for path in sys.argv[3:]:
        rows = read(path)
        kinds = collections.Counter(r[2].split("@")[0] for r in rows)
        print(f"| {path.split('rows_1b9e1_damage_')[-1]} | {len(rows)} | {len({r[1] for r in rows})} | {len(kinds)} | {sum(r[3] == '0' for r in rows)} | {sum(r[5] == 'eq' for r in rows)} | {sum(r[7] == 'clean' for r in rows)} |")


if __name__ == "__main__":
    main()
