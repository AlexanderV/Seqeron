#!/usr/bin/env python3
"""Paired, interleaved A/B of the DCS cost on dna_10m (E2 gate: slowdown < 25 %), fresh process per run, order alternates inside the pairs.
A = e1 (direct comparator, no rank table), B = e2 --dcs <X> (rank table built, then the same build); the metric is sa_plus_dcs_s of B against sa_s of A.
A/A = the same A configuration on both sides (noise of the session). Usage: python3 -I ab_dcs.py <dll> <data dir> <pairs> [X]
"""
import statistics
import subprocess
import sys


def run(dll, data, args):
    out = subprocess.run(["dotnet", dll] + args + ["--data", data], capture_output=True, text=True, check=True).stdout
    vals = {}
    for line in out.splitlines():
        c = line.split("\t")
        if len(c) >= 6 and c[0] == "R":
            vals[c[3]] = c[4]
    return vals


def main():
    dll, data, pairs = sys.argv[1], sys.argv[2], int(sys.argv[3])
    x = sys.argv[4] if len(sys.argv) > 4 else "1024"
    a_args = ["e1", "--family", "dna_10m", "--no-verify"]
    b_args = ["e2", "--family", "dna_10m", "--no-verify", "--dcs", x]
    metric = {"e1": "sa_s", "e2": "sa_plus_dcs_s"}
    aa, ab = [], []
    for i in range(pairs):
        first = float(run(dll, data, a_args)[metric["e1"]])
        second = float(run(dll, data, a_args)[metric["e1"]])
        aa.append(second / first if i % 2 == 0 else first / second)
    for i in range(pairs):
        if i % 2 == 0:
            a = float(run(dll, data, a_args)[metric["e1"]])
            b = float(run(dll, data, b_args)[metric["e2"]])
        else:
            b = float(run(dll, data, b_args)[metric["e2"]])
            a = float(run(dll, data, a_args)[metric["e1"]])
        ab.append(b / a)
    print(f"A/A ratios ({pairs} pairs): median {statistics.median(aa):.4f} range [{min(aa):.4f}, {max(aa):.4f}]")
    print(f"B/A ratios ({pairs} pairs, X={x}): median {statistics.median(ab):.4f} range [{min(ab):.4f}, {max(ab):.4f}]  -> DCS slowdown {100 * (statistics.median(ab) - 1):.1f} % (gate < 25 %)")
    print("pairs:", " ".join(f"{r:.3f}" for r in ab))


if __name__ == "__main__":
    main()
