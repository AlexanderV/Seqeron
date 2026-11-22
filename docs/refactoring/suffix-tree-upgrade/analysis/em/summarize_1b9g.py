#!/usr/bin/env python3
"""Tables of work package 1b.9g from the rows of the in-process A/B client (analysis/bench/abapp) and the four-thread rows of ReadBench (`mt:`).
   usage: summarize_1b9g.py <directory with the abapp row files> [<directory with the ReadBench rows>]
   Row files (one process = three lines "AB <op> A <n> B <n> B/A median-of-paired <r> ..."; the median over the processes of <r> is the figure, the range is the spread):
     null_*                  the same build twice (the noise floor and the bias of the second spec)
     before_idx_vs_file4_*   four threads, A = the index behind the facade, B = the facade over the container (order 1), the head before the change
     before_file_vs_idx4_*   the same with the order of the two specs exchanged (order 2)
     after_*                 the same for the final build
     single_base_vs_final_*  one thread, A = the head, B = the final build, both the facade over the container
     {before,after}_large_mem_vs_file_*  one thread, large pattern sets, A = the in-memory tree, B = the file-backed tree of the same build
   The balanced figure of the gate is sqrt(r1 / r2) with r1 the median of order 1 (B = facade) and r2 the median of order 2 (B = index): the bias of the second spec (it favours B) cancels."""
import glob, math, os, re, statistics, sys

OPS = ('count_present', 'locate', 'count_absent')
ROW = re.compile(r'^AB (\S+)\s+A\s+(\d+)\s+B\s+(\d+)\s+B/A median-of-paired ([\d.]+)')


def load(path):
    ops = {}
    for line in open(path, errors='replace'):
        m = ROW.match(line)
        if m:
            ops.setdefault(m.group(1), []).append((float(m.group(4)), int(m.group(2)), int(m.group(3))))
    return ops


def cell(ops, op):
    r = ops.get(op)
    if not r:
        return None
    v = [x[0] for x in r]
    return statistics.median(v), min(v), max(v), len(v)


def fmt(c):
    return '-' if c is None else '%.3f [%.3f-%.3f] n=%d' % c


def find(d, pattern):
    return sorted(glob.glob(os.path.join(d, pattern)))


def main():
    d = sys.argv[1]
    print('### Null experiment (the same build twice; B / A, median of the processes [range], n)')
    print('| file | count present | locate | count absent |')
    print('|---|---|---|---|')
    for p in find(d, 'null_*.txt'):
        o = load(p)
        print('| %s | %s |' % (os.path.basename(p), ' | '.join(fmt(cell(o, op)) for op in OPS)))
    print()
    print('### Gate: four threads, the facade over the container against the index behind it')
    print('Order 1: A = index, B = facade (B / A = facade / index). Order 2: A = facade, B = index. Balanced = sqrt(order 1 / order 2).')
    print('| build | text | order 1 (count / locate / absent) | order 2 (count / locate / absent) | balanced (count / locate / absent) |')
    print('|---|---|---|---|---|')
    for build in ('before', 'after'):
        for text in ('dna', 'protu', 'protr', 'protO'):
            f1 = find(d, '%s_idx_vs_file4_%s.txt' % (build, text))
            f2 = find(d, '%s_file_vs_idx4_%s.txt' % (build, text))
            c1 = [cell(load(f1[0]), op) if f1 else None for op in OPS]
            c2 = [cell(load(f2[0]), op) if f2 else None for op in OPS]
            if not any(c1) and not any(c2):
                continue
            bal = []
            for a, b in zip(c1, c2):
                bal.append('-' if a is None or b is None else '%.3f' % math.sqrt(a[0] / b[0]))
            short = lambda cs: ' / '.join('-' if c is None else '%.3f (%.3f-%.3f, n=%d)' % (c[0], c[1], c[2], c[3]) for c in cs)
            print('| %s | %s | %s | %s | %s |' % (build, text, short(c1), short(c2), ' / '.join(bal)))
    print()
    print('### One thread: the facade over the container, head (A) against the final build (B)')
    print('| text | count present | locate | count absent |')
    print('|---|---|---|---|')
    for p in find(d, 'single_base_vs_final_*.txt'):
        o = load(p)
        print('| %s | %s |' % (os.path.basename(p)[len('single_base_vs_final_'):-4], ' | '.join(fmt(cell(o, op)) for op in OPS)))
    print()
    print('### Large pattern sets (200,000 / 20,000 / 100,000), one thread: file-backed / in-memory of the same build')
    print('| build | text | count present | locate | count absent |')
    print('|---|---|---|---|---|')
    for build in ('before', 'after'):
        for p in find(d, '%s_large_mem_vs_file_*.txt' % build):
            o = load(p)
            print('| %s | %s | %s |' % (build, os.path.basename(p)[len(build + '_large_mem_vs_file_'):-4], ' | '.join(fmt(cell(o, op)) for op in OPS)))
    if len(sys.argv) > 2:
        print()
        print('### ReadBench mt: four threads on one tree, checked (facade) / ix-file (index behind it)')
        print('| rows file | processes | count | locate | absent |')
        print('|---|---|---|---|---|')
        for p in find(sys.argv[2], 'threads_*.tsv'):
            ratios = {1: [], 2: [], 3: []}
            for line in open(p):
                f = line.rstrip('\n').split('\t')
                if len(f) >= 11 and f[0] == 'mt' and f[3] == 'checked':
                    for k in (1, 2, 3):
                        ratios[k].append(float(f[7 + k]))
            if ratios[1]:
                print('| %s | %d | %s |' % (os.path.basename(p), len(ratios[1]), ' | '.join('%.3f [%.3f-%.3f]' % (statistics.median(ratios[k]), min(ratios[k]), max(ratios[k])) for k in (1, 2, 3))))


if __name__ == '__main__':
    main()
