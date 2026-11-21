#!/usr/bin/env python3
"""Tables of the kill campaign of run_1b9e_kill.py from its rows: kills per kind and phase, identical or not, and the time that a resume saves per phase.
usage: summarize_1b9e_kill.py <rows file> [<rows file> ...]"""
import collections
import statistics
import sys

COLS = "kind n id mode aimed last phase killed target skipped resume_s fresh_s saved_s hash_ok leftover_ok extra".split()
ORDER = ["ph0", "ph1", "ph23", "ph4", "bwt", "psi", "ph7", "seal"]


def main():
    rows = []
    for path in sys.argv[1:]:
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) >= 15:
                rows.append(dict(zip(COLS, p)))
    by = collections.defaultdict(list)
    for r in rows:
        by[(r["kind"], r["n"])].append(r)
    print("| input | kills (signal landed) | boundary | random | chained | finished before the signal | target after the kill: none / complete / damaged | hash equal | nothing else left |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    total = collections.Counter()
    for (kind, n), rs in sorted(by.items()):
        k = [r for r in rs if r["killed"] == "1"]
        t = collections.Counter(r["target"] for r in rs)
        print(f"| {kind} (n = {n}) | {len(k)} | {sum(r['mode'] == 'boundary' for r in k)} | {sum(r['mode'] == 'random' for r in k)} | {sum(r['extra'].startswith('chain-killed') for r in rs)} | {sum(r['killed'] == '0' for r in rs)} | {t['none']} / {t['complete']} / {t['damaged']} | {sum(r['hash_ok'] == '1' for r in rs)} of {len(rs)} | {sum(r['leftover_ok'] == '1' for r in rs)} of {len(rs)} |")
        total["kills"] += len(k)
        total["rows"] += len(rs)
        total["hash"] += sum(r["hash_ok"] == "1" for r in rs)
        total["left"] += sum(r["leftover_ok"] == "1" for r in rs)
    print(f"\nTotal: {total['kills']} kills landed in {total['rows']} cases; hash equal in {total['hash']}, nothing else left in {total['left']}.\n")
    print("| input | phase of the kill (last marker seen) | kills | stages skipped by the resume (most common) | resume s (median [min-max]) | uninterrupted build s | saved s (median [min-max]) |")
    print("| --- | --- | --- | --- | --- | --- | --- |")
    for (kind, n), rs in sorted(by.items()):
        for ph in ORDER:
            sel = [r for r in rs if r["phase"] == ph and r["killed"] == "1" and r["resume_s"]]
            if not sel:
                continue
            res = sorted(float(r["resume_s"]) for r in sel)
            sav = sorted(float(r["saved_s"]) for r in sel)
            common = collections.Counter(r["skipped"] or "(none)" for r in sel).most_common(1)[0]
            print(f"| {kind} | {ph} | {len(sel)} | {common[0]} ({common[1]}) | {statistics.median(res):.2f} [{res[0]:.2f}-{res[-1]:.2f}] | {sel[0]['fresh_s']} | {statistics.median(sav):.2f} [{sav[0]:.2f}-{sav[-1]:.2f}] |")
    bad = [r for r in rows if r["hash_ok"] != "1" or r["leftover_ok"] != "1"]
    for r in bad:
        print("FAILED:", r)


if __name__ == "__main__":
    main()
