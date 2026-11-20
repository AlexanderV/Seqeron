#!/usr/bin/env python3
"""GENOMIC-REPEAT-001 oracle (A3-13): EMBOSS palindrome 6.6.0 with -gaplimit 0 -nummismatches 0 (contiguous reverse-complement palindromes) vs GenomicAnalyzer.FindPalindromes.

Random DNA (30-150 nt, planted palindromes with arm 3-9, optional lower case), `palindrome -minpallen m -maxpallen M -gaplimit 0 -nummismatches 0 -overlap`, 2 <= m <= 5, M in {6, 12, 20}.
The report lists each hit as two arms `start1 arm1 end1 / start2 arm2 end2` (1-based, the second arm read on the other strand); with no loop the hit is the palindromic substring start1 .. start2.
TSV (tests/.../TestData/Analysis/palindrome_gap0_oracle.tsv):  sequence \t m \t M \t hits ('start:end:arm' 1-based, sorted, space separated; '-' none)
Usage: python3 scripts/analysis/generate_palindrome_gap0_oracle.py [rows] > palindrome_gap0_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261126)
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 330
    done = 0
    with tempfile.TemporaryDirectory() as tmp:
        for _ in range(n):
            s = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(30, 150)))
            for _ in range(rnd.randint(0, 4)):
                arm = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(3, 9)))
                p = rnd.randrange(len(s)); s = s[:p] + arm + rc(arm) + s[p:]
            m, M = rnd.choice([2, 3, 4, 5]), rnd.choice([6, 12, 20])
            fa, out = os.path.join(tmp, 'p.fa'), os.path.join(tmp, 'p.out')
            open(fa, 'w').write('>p\n' + s + '\n')
            r = subprocess.run(['palindrome', '-sequence', fa, '-minpallen', str(m), '-maxpallen', str(M), '-gaplimit', '0', '-nummismatches', '0', '-overlap', '-outfile', out, '-auto'], capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(out):
                sys.stderr.write('FAIL %s: %s\n' % (s, r.stderr.strip()[:200])); continue
            hits = [(int(a), int(d), len(b)) for a, b, c, d, e, f in re.findall(r'^(\d+)\s+([acgt]+)\s+(\d+)\n\s+\|+\n(\d+)\s+([acgt]+)\s+(\d+)', open(out).read(), re.M)]
            hits.sort()
            if rnd.random() < 0.2: s = s.lower()
            print('%s\t%d\t%d\t%s' % (s, m, M, ' '.join('%d:%d:%d' % h for h in hits) or '-'))
            done += 1
    sys.stderr.write('rows: %d\n' % done)


if __name__ == '__main__':
    main()
