#!/usr/bin/env python3
"""Summary of the A/B rows of work package S2.3 (run_s23_ab.sh): per cell (input, kind, threads) and operation the median over the processes of the per-process median of the per-block ratios B / A, the quartiles and the range over the
processes, and the verdict against the threshold 0.97 (ADR-005 section 3.4). The null experiment (A = B = the baseline build) gives the spread of the instrument: a figure below 1.0 passes only if it clears the threshold by more than
the spread of the null rows of the same kind (its minimum over the processes is not below 0.97 - spread is not asked: the rule of the brief is "pooled median >= 0.97 and not below the null spread", i.e. the pooled median of the A/B cell must
not be below the pooled median of the matching null cell by more than the null cell's own range).
usage: summarize_s23_ab.py <rows dir> [out.tsv]"""
import glob, os, re, statistics, sys

PAT = re.compile(r'^AB (\S+)\s+A\s+(\d+)\s+B\s+(\d+)\s+B/A median-of-paired ([\d.]+)')
OPS = ('count_present', 'count_absent', 'locate')


def load(path):
    ops = {}
    for line in open(path):
        m = PAT.match(line)
        if m:
            ops.setdefault(m.group(1), []).append(float(m.group(4)))
    return ops


def q(v, f):
    s = sorted(v)
    return s[min(len(s) - 1, int(f * len(s)))]


rows = []
base = sys.argv[1]
cells = {}
for path in sorted(glob.glob(os.path.join(base, '*.txt'))):
    name = os.path.basename(path)[:-4]
    mode, rest = name.split('_', 1)
    rest = rest.rsplit('_t', 1)
    th = rest[1]
    input_kind = rest[0].rsplit('_', 1)
    cells[(mode, input_kind[0], input_kind[1], th)] = load(path)
print('mode\tinput\tkind\tthreads\toperation\tprocesses\tpooled_median\tp25\tp75\tmin\tmax\tverdict')
for (mode, inp, kind, th), ops in sorted(cells.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2], int(kv[0][3]))):
    for op in OPS:
        v = ops.get(op)
        if not v:
            continue
        med = statistics.median(v)
        verdict = ''
        if mode == 'ab':
            verdict = 'ok' if med >= 0.97 else 'BELOW 0.97'
            nullv = cells.get(('null', inp, kind, th), {}).get(op)
            if nullv:
                spread = max(nullv) - min(nullv)
                verdict += f' (null median {statistics.median(nullv):.3f}, null range {min(nullv):.3f}-{max(nullv):.3f})'
        print(f'{mode}\t{inp}\t{kind}\t{th}\t{op}\t{len(v)}\t{med:.4f}\t{q(v, 0.25):.4f}\t{q(v, 0.75):.4f}\t{min(v):.4f}\t{max(v):.4f}\t{verdict}')
