#!/usr/bin/env python3
"""Work package 1b.9, part C: the table of the file-backed production chain from the logs of run_1b9_session3.sh: one line per (kind, n) with the builds of the log(s) (the chain's own and the --rebuild repeats):
the time inside SesaBuilder.Build, the RssAnon peak of the independent sampler, the container size, the verdict of the checker and the hash comparison.
usage: chain_table_1b9.py <log>...   (the logs are read in order and their builds are pooled per (kind, n))"""
import re
import statistics
import sys


def main(paths):
    d = {}
    for path in paths:
        cur = None
        for line in open(path):
            m = re.match(r"== \S+ (chain|rebuild) (\S+) (\S+)(?: #(\d))?", line)
            if m:
                cur = (m.group(2), m.group(3))
                d.setdefault(cur, {"total": [], "rss": [], "ph": {}, "size": None, "sum": [], "struct": [], "ref": []})
                continue
            if cur is None:
                continue
            e = d[cur]
            for key, pat in (("total", r"\btotal_s\s+(\S+) s"), ("rss", r"independent sampler: RssAnon peak (\S+) MB")):
                mm = re.search(pat, line)
                if mm:
                    e[key].append(float(mm.group(1)))
            mm = re.search(r"file_size\s+(\d+) bytes", line)
            if mm:
                e["size"] = int(mm.group(1))
            mm = re.search(r"SUMMARY: (\d+) PASS, (\d+) FAIL", line)
            if mm:
                e["sum"].append(f"{mm.group(1)} PASS, {mm.group(2)} FAIL")
            mm = re.search(r"S1 Structure level of the reader: pass in (\S+) s", line)
            if mm:
                e["struct"].append(float(mm.group(1)))
            mm = re.search(r"I1 container SHA-256 (\S+)\.\.\. equals the reference", line)
            if mm:
                e["ref"].append(mm.group(1))
    print("| kind | n | builds | time inside SesaBuilder.Build s: median [min-max] | RssAnon peak MB (independent sampler): median [min-max] | container bytes (B/char) | Structure s | checker | = reference (I1) |")
    print("|---|---|---|---|---|---|---|---|---|")
    for (kind, n), e in d.items():
        nn = {"10m": 10**7, "100m": 10**8}[n]
        t, r = e["total"], e["rss"]
        st = f"{statistics.median(e['struct']):.1f}" if e["struct"] else "-"
        print(f"| {kind} | {n} | {len(t)} | {statistics.median(t):.2f} [{min(t):.2f}-{max(t):.2f}] | {statistics.median(r):.1f} [{min(r):.1f}-{max(r):.1f}] | {e['size']:,} ({e['size'] / nn:.2f}) | {st} | "
              f"{', '.join(sorted(set(e['sum'])))} ({len(e['sum'])} chain{'s' if len(e['sum']) != 1 else ''}) | {'yes' if e['ref'] else '-'} |")


main(sys.argv[1:])
