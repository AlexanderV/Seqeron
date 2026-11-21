#!/usr/bin/env python3
"""Tables of run_1b9e1_timing.py: per input and kill point the median [min-max] of the seconds of the restart, of the resume of the commit before the package (old), of the resume of this package (new) and of the resume of this
package on a file whose every byte outside the kept ranges was overwritten with 0xA5 (garbage, the worst case), the seconds of the normalisation itself, the bytes it read and wrote, and the device writes; from the raw rows (no number is typed by hand).
usage: summarize_1b9e1_timing.py <rows file> [<rows file> ...]"""
import collections
import statistics
import sys


def main():
    cells = collections.defaultdict(list)
    order = []
    for path in sys.argv[1:]:
        for line in open(path):
            p = line.rstrip("\n").split("\t")
            if len(p) < 13:
                continue
            kind, marker, mode = p[0], p[1], p[2]
            if (kind, marker) not in order:
                order.append((kind, marker))
            cells[(kind, marker, mode)].append(p)

    def stat(rows, i, fmt="{:.2f}"):
        v = sorted(float(r[i]) for r in rows if r[i] != "")
        if not v:
            return "-"
        return (fmt + " [" + fmt + "-" + fmt + "]").format(statistics.median(v), v[0], v[-1])

    def med(rows, i):
        v = [float(r[i]) for r in rows if r[i] != ""]
        return statistics.median(v) if v else float("nan")

    print("| input | killed at | kept (skipped by the resume) | restart s | resume, before s | resume, this package s | resume, garbage s | normalisation s | read MB | written MB | written MB (garbage) | device MB: restart / before / this / garbage | saved by this resume s |")
    print("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
    for kind, marker in order:
        rs, ol, nw, gb = (cells.get((kind, marker, m), []) for m in ("restart", "old", "new", "garbage"))
        if not nw:
            continue
        kept = " or ".join(sorted({r[12] or "(none)" for r in nw}, key=len))  # the kills of one kill point do not always land at the same side of a boundary
        saved = med(rs, 4) - med(nw, 4)
        print(f"| {kind} | {marker} | {kept} | {stat(rs, 4)} | {stat(ol, 4)} | {stat(nw, 4)} | {stat(gb, 4)} | {stat(nw, 8, '{:.3f}')} | {med(nw, 9) / 1e6:.0f} | {med(nw, 10) / 1e6:.1f} | {med(gb, 10) / 1e6:.0f} | "
              f"{med(rs, 6):.0f} / {med(ol, 6):.0f} / {med(nw, 6):.0f} / {med(gb, 6):.0f} | {saved:.2f} |")
    bad = [r for rows in cells.values() for r in rows if r[11] != "1"]
    print(f"\nruns: {sum(len(v) for v in cells.values())}, every one finished with exit 0, the bytes of the uninterrupted build and nothing else left: {'yes' if not bad else 'NO: ' + str(len(bad))}")


if __name__ == "__main__":
    main()
