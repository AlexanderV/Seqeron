#!/usr/bin/env python3
"""GENOMIC-TANDEM-001 oracle (A2-34): MUMmer `exact-tandems` (repeat-match -t + tandem-repeat.awk), MUMmer 3.23.

`exact-tandems <file> <min-match-len>` lists, for every start position, the exact tandem repeat with the
smallest unit length whose maximal match (repeat-match -t -n L) is at least L long:
    Start (1-based)  Extent  UnitLen  Copies(= Extent/UnitLen)
The Debian wrapper /usr/bin/exact-tandems is broken (a stray line break in its `awk -f $scriptdir/...`
line), so this script runs the two documented stages itself, exactly as the wrapper intends:
    repeat-match -t -n L f.fa | tail -n +3 | sort -k1n -k2n | awk -f /usr/lib/mummer/tandem-repeat.awk

Output TSV (tests/Seqeron/Seqeron.Genomics.Tests/TestData/Analysis/exact_tandems_oracle.tsv), one row per
(sequence, L):  sequence  L  rows
  rows = start0:unitLen:extent:canonicalUnit;...   or "-"
  start0 = Start-1 (0-based), canonicalUnit = lexicographically least rotation of seq[start0:start0+unitLen]
  (computed here by brute force, independent of Seqeron).

Usage: python3 scripts/analysis/generate_exact_tandems_oracle.py > exact_tandems_oracle.tsv
Never aborts: a failing external call yields the row "ERR" (the test then fails loudly).
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261119)
AWK = '/usr/lib/mummer/tandem-repeat.awk'


def rand(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


def mutate(s, alpha, rate):
    return ''.join(rnd.choice(alpha) if rnd.random() < rate else c for c in s)


def planted(alpha):
    parts = [rand(rnd.randint(0, 12), alpha)]
    for _ in range(rnd.randint(1, 4)):
        unit = rand(rnd.randint(1, 8), alpha)
        run = unit * rnd.randint(2, 6) + unit[:rnd.randint(0, len(unit))]
        parts.append(mutate(run, alpha, 0.04) if rnd.random() < 0.3 else run)
        parts.append(rand(rnd.randint(0, 12), alpha))
    return ''.join(parts)


def least_rotation(u):
    return min(u[i:] + u[:i] for i in range(len(u)))


def exact_tandems(seq, n):
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 's.fa')
        open(f, 'w').write('>s\n' + seq + '\n')
        try:
            rm = subprocess.run(['repeat-match', '-t', '-n', str(n), f], capture_output=True, text=True, timeout=60)
            lines = rm.stdout.splitlines()[2:]
            text = '\n'.join(lines) + ('\n' if lines else '')
            srt = subprocess.run(['sort', '-k1n', '-k2n'], input=text, capture_output=True, text=True, timeout=60)
            aw = subprocess.run(['awk', '-f', AWK], input=srt.stdout, capture_output=True, text=True, timeout=60)
        except Exception as e:  # never abort
            sys.stderr.write('ERR %s %s\n' % (seq, e))
            return None
    rows = []
    for ln in aw.stdout.splitlines()[1:]:
        start, extent, unit, _ = ln.split()
        start0, extent, unit = int(start) - 1, int(extent), int(unit)
        rows.append('%d:%d:%d:%s' % (start0, unit, extent, least_rotation(seq[start0:start0 + unit])))
    return rows


def main():
    seqs = []
    for alpha in ('ACGT', 'AC', 'A', 'AT', 'ACGT'):
        for _ in range(240):
            seqs.append(planted(alpha) if rnd.random() < 0.8 else rand(rnd.randint(4, 60), alpha))
    seqs += ['ATATATAT', 'AAAAA', 'AAAGAAAG', 'ACGACGACGTT', 'GGGGGGGGGG', 'ACACACGTGTGT', 'AATAATAATAAT', 'TTAGGGTTAGGGTTAGGG', 'ACGT', 'A']
    count = 0
    for seq in seqs:
        if not seq:
            continue
        for n in (1, 2, 4):
            rows = exact_tandems(seq, n)
            print('%s\t%d\t%s' % (seq, n, 'ERR' if rows is None else (';'.join(rows) if rows else '-')))
            count += 1
    sys.stderr.write('rows=%d\n' % count)


if __name__ == '__main__':
    main()
