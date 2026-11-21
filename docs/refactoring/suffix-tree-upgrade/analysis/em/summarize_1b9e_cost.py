#!/usr/bin/env python3
"""Tables of run_1b9e_cost.sh: per label the median [min-max] of the pair ratios and of the seconds, from the raw rows (no number is typed by hand).
usage: summarize_1b9e_cost.py <rows file> [<rows file> ...]"""
import collections
import statistics
import sys


def main():
    series = collections.defaultdict(lambda: collections.defaultdict(list))
    for path in sys.argv[1:]:
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if p[0] != "R" or len(p) < 6 or not p[1].startswith("resume-"):
                continue
            op, metric, value, notes = p[1], p[2], p[3], p[5]
            label = notes.split(" ")[0]
            if op == "resume-ab" and metric in ("total_s", "dev_sectors_written_synced", "flushes", "ckpt_total_s"):
                base, pair, which = label.rsplit(".", 2)
                if pair == "0":
                    continue  # the warm-up pair of the process
                series[(op, base)][f"{metric}.{which}"].append(float(value))
            elif op == "resume-abhead" and metric == "total_s":
                base, pair, which = label.rsplit(".", 2)
                if pair == "0":
                    continue
                series[(op, base)][f"total_s.{which}"].append(float(value))
    def stat(xs):
        xs = sorted(xs)
        return f"{statistics.median(xs):.3f} [{xs[0]:.3f}-{xs[-1]:.3f}]"
    print("| experiment | label | n pairs | seconds A (median [min-max]) | seconds B | B / A per pair (median [min-max]) | device MB written A | device MB written B | fsync calls B | recording s (B) |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for (op, base), d in sorted(series.items()):
        if op == "resume-ab":
            a, b = d["total_s.off"], d["total_s.on"]
            ratios = [y / x for x, y in zip(a, b)]
            sa = d.get("dev_sectors_written_synced.off", [0])
            sb = d.get("dev_sectors_written_synced.on", [0])
            print(f"| Resume off (A) vs on (B) | {base} | {len(a)} | {stat(a)} | {stat(b)} | {stat(ratios)} | {statistics.median(sa) * 512 / 1e6:.0f} | {statistics.median(sb) * 512 / 1e6:.0f} | {statistics.median(d.get('flushes.on', [0])):.0f} | {stat(d['ckpt_total_s.on']) if d.get('ckpt_total_s.on') else ''} |")
        else:
            a, b = d["total_s.before"], d["total_s.after"]
            ratios = [y / x for x, y in zip(a, b)]
            print(f"| the build before the package (A) vs this commit, Resume off (B) | {base} | {len(a)} | {stat(a)} | {stat(b)} | {stat(ratios)} | | | | |")


if __name__ == "__main__":
    main()
