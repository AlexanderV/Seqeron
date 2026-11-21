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
    rss = [float(x.split("rssanon_mb=")[1].split()[0]) for r in rows for x in [r["extra"]] if "rssanon_mb=" in x]
    print(f"\nTotal: {total['kills']} kills landed in {total['rows']} cases; hash equal in {total['hash']}, nothing else left in {total['left']}." + (f" RssAnon of the finishing runs: max {max(rss):.1f} MB over {len(rss)} runs (the gate is 320)." if rss else "") + "\n")
    print("Kills by the phase that the process was in (the last marker that it printed before the signal) and what the resume skipped (the most common set of stages). The seconds of the resumes are in the raw rows but are not a measurement of the saving: the resumes of these campaigns run with 1, 2 or 4 threads at random and right after the kill; the saving is measured by run_1b9e_timing.py.\n")
    print("| input | " + " | ".join(ORDER) + " |")
    print("| --- | " + " | ".join("---" for _ in ORDER) + " |")
    for (kind, n), rs in sorted(by.items()):
        cells = []
        for ph in ORDER:
            sel = [r for r in rs if r["phase"] == ph and r["killed"] == "1"]
            if not sel:
                cells.append("0")
                continue
            common = collections.Counter(r["skipped"] or "(none)" for r in sel).most_common(1)[0]
            cells.append(f"{len(sel)} ({common[0]})")
        print(f"| {kind} | " + " | ".join(cells) + " |")
    bad = [r for r in rows if r["hash_ok"] != "1" or r["leftover_ok"] != "1"]
    for r in bad:
        print("FAILED:", r)


if __name__ == "__main__":
    main()
