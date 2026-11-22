#!/usr/bin/env python3
"""Tables of work package 1b.9g.1 (the flags of the gate in an object of their own): the rows of the in-process A/B client `abapp` and of `ReadBench mt:` taken after the change, with the rows of 1b.9g (the flags in the gate) for comparison.
   usage: summarize_1b9g1.py <directory of the rows of 1b.9g (res3)> <directory of the rows of 1b.9g.1 (res6)>
   Row files of the second directory: null_file_<text>_new.txt (the final build twice), after_idx_vs_file4_<text>.txt (order 1: A = index, B = facade) and after_file_vs_idx4_<text>.txt (order 2), single_{base,prev}_vs_new_<text>.txt (one thread,
   the facade over the container, A = the head of 1b.9g (base) or the head before this change (prev), B = this build), large_cross_base_vs_new_<text>.txt, and rbalt_<text>_{prev,new}.tsv (ReadBench mt: one row per process, the facade over the
   index behind it, `ratio_to_first`), rbdiag_dna_{prev,new}.tsv (the diagnostic variant `ix-gate`), large4_*_dna.txt (the large sets at four threads). The balanced figure is sqrt(order 1 / order 2) as in summarize_1b9g.py.
   In the repository the rows carry the prefixes `rows_1b9g_ab_` (1b.9g) and `rows_1b9g1_` (this package); strip them into two directories to run this."""
import math, os, re, statistics, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize_1b9g import OPS, cell, find, fmt, load  # noqa: E402


def short(cs):
    return ' / '.join('-' if c is None else '%.3f (%.3f-%.3f, n=%d)' % (c[0], c[1], c[2], c[3]) for c in cs)


def balanced(c1, c2):
    return ' / '.join('-' if a is None or b is None else '%.3f' % math.sqrt(a[0] / b[0]) for a, b in zip(c1, c2))


def gate_rows(d, build_prefix, text):
    f1 = find(d, '%s_idx_vs_file4_%s.txt' % (build_prefix, text)) + find(d, '%s_idx_vs_file4_%s_more.txt' % (build_prefix, text))
    f2 = find(d, '%s_file_vs_idx4_%s.txt' % (build_prefix, text)) + find(d, '%s_file_vs_idx4_%s_more.txt' % (build_prefix, text))
    c1 = [cell(load(f1), op) if f1 else None for op in OPS]
    c2 = [cell(load(f2), op) if f2 else None for op in OPS]
    return c1, c2


def readbench(*paths):
    """The `checked` rows of one or several files (the processes of all of them count)."""
    ratios = {1: [], 2: [], 3: []}
    for path in paths:
        if not os.path.exists(path):
            continue
        for line in open(path):
            f = line.rstrip('\n').split('\t')
            if len(f) >= 11 and f[0] == 'mt' and f[3] == 'checked':
                for k in (1, 2, 3):
                    ratios[k].append(float(f[7 + k]))
    return ratios


def main():
    old, new = sys.argv[1], sys.argv[2]
    print('### Null experiment (the final build twice, four threads; B / A, median of the processes [range], n)')
    print('| file | count present | locate | count absent |')
    print('|---|---|---|---|')
    for p in find(new, 'null_*.txt'):
        print('| %s | %s |' % (os.path.basename(p), ' | '.join(fmt(cell(load(p), op)) for op in OPS)))
    print()
    print('### Gate: four threads, the facade over the container against the index behind it (abapp)')
    print('Order 1: A = index, B = facade (B / A = facade / index). Order 2: A = facade, B = index. Balanced = sqrt(order 1 / order 2). "1b.9g" = the flags in the gate object (head ea11f848), "1b.9g.1" = the flags in an object of their own.')
    print('| build | text | order 1 (count / locate / absent; median, range, n) | order 2 | balanced |')
    print('|---|---|---|---|---|')
    for text in ('dna', 'protu', 'protr', 'protO'):
        for label, d, prefix in (('1b.9g', old, 'after'), ('1b.9g.1', new, 'after')):
            c1, c2 = gate_rows(d, prefix, text)
            if not any(c1) and not any(c2):
                continue
            print('| %s | %s | %s | %s | %s |' % (label, text, short(c1), short(c2), balanced(c1, c2)))
    print()
    print('### One thread, the facade over the container (B / A above 1: this build is faster)')
    print('| A | text | count present | locate | count absent |')
    print('|---|---|---|---|---|')
    def pooled(pattern, prefix):
        """{text: [files]} for the files `<prefix><text>.txt` and `<prefix><text>_more.txt`."""
        groups = {}
        for p in find(new, pattern):
            name = os.path.basename(p)[len(prefix):-4]
            groups.setdefault(name[:-5] if name.endswith('_more') else name, []).append(p)
        return groups

    for tag, name in (('base', 'the head of 1b.9g (fd9bd0d8)'), ('prev', 'the head before this change (ea11f848)')):
        for text, fs in pooled('single_%s_vs_new_*.txt' % tag, 'single_%s_vs_new_' % tag).items():
            print('| %s | %s | %s |' % (name, text, ' | '.join(fmt(cell(load(fs), op)) for op in OPS)))
    print()
    print('### Large pattern sets (200,000 / 20,000 / 100,000), one thread: the file-backed tree of the head of 1b.9g (A) against this build (B)')
    print('| text | count present | locate | count absent |')
    print('|---|---|---|---|')
    for text, fs in pooled('large_cross_base_vs_new_*.txt', 'large_cross_base_vs_new_').items():
        print('| %s | %s |' % (text, ' | '.join(fmt(cell(load(fs), op)) for op in OPS)))
    print()
    print('### ReadBench mt: four threads on one tree, checked (facade) / ix-file (index behind it); one row per process, the two builds alternate')
    print('| text | build | processes | count | locate | absent |')
    print('|---|---|---|---|---|---|')
    for text in ('dna', 'protu', 'protr'):
        for build, label in (('prev', '1b.9g (flags in the gate)'), ('new', '1b.9g.1 (flags in their own object)')):
            r = readbench(os.path.join(new, 'rbalt_%s_%s.tsv' % (text, build)), os.path.join(new, 'rbalt_%s_b_%s.tsv' % (text, build)))
            if r[1]:
                print('| %s | %s | %d | %s |' % (text, label, len(r[1]), ' | '.join('%.3f [%.3f-%.3f]' % (statistics.median(r[k]), min(r[k]), max(r[k])) for k in (1, 2, 3))))
    print()
    print('### Diagnostic: ReadBench mt: three variants in one process (index, index inside the gate of the build, facade), DNA; B / index')
    print('| build | variant | processes | count | locate | absent |')
    print('|---|---|---|---|---|---|')
    for build, label in (('prev', '1b.9g (flags in the gate)'), ('new', '1b.9g.1 (flags in their own object)')):
        rows = {}
        p = os.path.join(new, 'rbdiag_dna_%s.tsv' % build)
        if os.path.exists(p):
            for line in open(p):
                f = line.rstrip('\n').split('\t')
                if len(f) >= 11 and f[0] == 'mt' and f[3] != 'ix-file':
                    rows.setdefault(f[3], []).append([float(f[8]), float(f[9]), float(f[10])])
        for v, r in rows.items():
            print('| %s | %s | %d | %s |' % (label, v, len(r), ' | '.join('%.3f [%.3f-%.3f]' % (statistics.median(x[k] for x in r), min(x[k] for x in r), max(x[k] for x in r)) for k in range(3))))
    print()
    print('### The large pattern sets at four threads, the facade over the container against the index behind it (abapp, memory-bound regime)')
    print('| order | count present | locate | count absent |')
    print('|---|---|---|---|')
    o1 = find(new, 'large4_idx_vs_file_*.txt')
    o2 = find(new, 'large4_file_vs_idx_*.txt')
    for lbl, fs in (('order 1 (A = index, B = facade)', o1), ('order 2 (A = facade, B = index)', o2)):
        if fs:
            print('| %s | %s |' % (lbl, ' | '.join(fmt(cell(load(fs), op)) for op in OPS)))
    if o1 and o2:
        print('| balanced | %s |' % balanced([cell(load(o1), op) for op in OPS], [cell(load(o2), op) for op in OPS]).replace(' / ', ' | '))
    print()
    print('Per process (count / locate / absent), in the order of the run:')
    for text in ('dna', 'protu', 'protr'):
        for build in ('prev', 'new'):
            r = readbench(os.path.join(new, 'rbalt_%s_%s.tsv' % (text, build)), os.path.join(new, 'rbalt_%s_b_%s.tsv' % (text, build)))
            if r[1]:
                print('- %s %s: %s' % (text, build, '; '.join('%.2f %.2f %.2f' % (a, b, c) for a, b, c in zip(r[1], r[2], r[3]))))


if __name__ == '__main__':
    main()
