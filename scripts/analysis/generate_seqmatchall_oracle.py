#!/usr/bin/env python3
"""GENOMIC-COMMON-001 oracle (A3-13): EMBOSS seqmatchall 6.6.0 on random sequence sets.

Sets of 2-5 sequences (nucleotide ACGT / AC / with N, and protein), related by point mutations, shared segments and unrelated flanks (20-200 residues), word sizes 2..12;
`seqmatchall -sequence s.fa -wordsize K -aformat match`: per pair (i < j) the lines `length name1 strand start1..end1 name2 strand start2..end2` (1-based).
TSV (tests/.../TestData/Analysis/seqmatchall_oracle.tsv):  sequences (comma separated) \t wordsize \t matches
matches = space separated `i:j:start1:start2:length` (0-based, sorted), '-' when none.
Usage: python3 scripts/analysis/generate_seqmatchall_oracle.py [rows] > seqmatchall_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261123)


def rnd_seq(n, al):
    return ''.join(rnd.choice(al) for _ in range(n))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 330
    done = 0
    with tempfile.TemporaryDirectory() as tmp:
        for _ in range(n):
            protein = rnd.random() < 0.25
            al = 'ACDEFGHIKLMNPQRSTVWY' if protein else rnd.choice(['ACGT', 'ACGT', 'AC', 'ACGTN'])
            base = rnd_seq(rnd.randint(20, 160), al)
            seqs = [base]
            for _ in range(rnd.randint(1, 4)):
                s = ''.join(c if rnd.random() > 0.08 else rnd.choice(al) for c in base)
                if rnd.random() < 0.6:
                    s = rnd_seq(rnd.randint(0, 25), al) + s[rnd.randint(0, 15):] + rnd_seq(rnd.randint(0, 15), al)
                if rnd.random() < 0.2:
                    s = rnd_seq(rnd.randint(2, 60), al)
                seqs.append(s)
            if rnd.random() < 0.3:
                seqs = [s.lower() if rnd.random() < 0.5 else s for s in seqs]
            k = rnd.choice([2, 3, 4, 5, 6, 8, 10, 12]) if not protein else rnd.choice([2, 3, 4, 5, 6])
            fa = os.path.join(tmp, 's.fa')
            with open(fa, 'w') as f:
                for i, s in enumerate(seqs):
                    f.write('>s%d\n%s\n' % (i, s))
            out = os.path.join(tmp, 's.out')
            r = subprocess.run(['seqmatchall', '-sequence', fa, '-wordsize', str(k), '-outfile', out, '-aformat', 'match', '-auto'], capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(out):
                sys.stderr.write('FAIL %s k=%d: %s\n' % (seqs, k, r.stderr.strip()[:200]))
                continue
            matches = []
            cur = None
            a = b = None
            for line in open(out):
                if line.startswith('# 1: s'): a = int(line.split('s')[-1])
                elif line.startswith('# 2: s'):
                    b = int(line.split('s')[-1]); cur = (a, b)
                m = re.match(r'\s*(\d+) s(\d+)\s+\+\s+(\d+)\.\.(\d+)\s+s(\d+)\s+\+\s+(\d+)\.\.(\d+)', line)
                if m and cur:
                    matches.append((cur[0], cur[1], int(m.group(3)) - 1, int(m.group(6)) - 1, int(m.group(1))))
            matches.sort()
            print('%s\t%d\t%s' % (','.join(seqs), k, ' '.join('%d:%d:%d:%d:%d' % t for t in matches) or '-'))
            done += 1
    sys.stderr.write('rows: %d\n' % done)


if __name__ == '__main__':
    main()
