#!/usr/bin/env python3
"""GENOMIC-MOTIFS-001 oracle (A3-12): EMBOSS fuzztran 6.6.0 (protein pattern search after translation) on random DNA.

DNA 24-300 nt (ACGT, some N, lower case); a pattern drawn from the protein language and, in most rows, built around a window of a translated frame (transeq `-frame 6`
of the same table) so that hits occur; `-frame` one of 1 2 3 F -1 -2 -3 R 6, `-table` one of 0 1 2 3 4 5 9 11 (EMBOSS menu: 0 = standard).
`fuzztran -sequence s.fa -pattern P -frame F -table T -rformat gff` -> per hit: GFF start / end (as printed), strand, frame notes: *frame, *start, *end (protein), *translation.
TSV (tests/.../TestData/Analysis/fuzztran_oracle.tsv):  sequence \t pattern \t frame \t table \t hits (`frame:pstart:pend:translation:gffstart:gffend` space separated, sorted; '-' none; 'ERR' invalid pattern)
Usage: python3 scripts/analysis/generate_fuzztran_oracle.py [rows] > fuzztran_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

rnd = random.Random(20261128)
AA = 'ACDEFGHIKLMNPQRSTVWY'


def window_pattern(prot):
    stretch = ''.join(c for c in prot if c in AA)
    if len(stretch) < 6:
        return None
    L = rnd.randint(2, 5); p = rnd.randrange(0, len(stretch) - L)
    w = stretch[p:p + L]
    parts = []
    for c in w:
        k = rnd.random()
        if k < 0.55: parts.append(c)
        elif k < 0.75: parts.append('x')
        elif k < 0.9: parts.append('[' + c + rnd.choice(AA) + ']')
        else: parts.append('{' + rnd.choice([a for a in AA if a != c]) + '}')
    if rnd.random() < 0.3:
        parts.insert(rnd.randrange(len(parts) + 1), 'x(1,2)')
    return ''.join(parts)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 760
    done = 0
    with tempfile.TemporaryDirectory() as tmp:
        for _ in range(n):
            s = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(24, 300)))
            if rnd.random() < 0.2:
                s = ''.join('N' if rnd.random() < 0.03 else c for c in s)  # N only: EMBOSS resolves the other IUPAC codons with its own table (trnComplete), which differs from the Biopython-style GeneticCode.Translate for 13 codon patterns (TRANS-SIXFRAME-001)
            if rnd.random() < 0.2:
                s = s.lower() if rnd.random() < 0.5 else ''.join(c.lower() if rnd.random() < 0.5 else c for c in s)
            table = rnd.choice([0, 0, 0, 1, 2, 3, 4, 5, 9, 11])
            frame = rnd.choice(['1', '2', '3', 'F', '-1', '-2', '-3', 'R', '6', '6', 'F'])
            fa = os.path.join(tmp, 's.fa'); open(fa, 'w').write('>s\n' + s + '\n')
            pat = None
            if rnd.random() < 0.85:
                tr = os.path.join(tmp, 't.out')
                subprocess.run(['transeq', '-sequence', fa, '-frame', '6', '-table', str(table), '-outseq', tr, '-auto'], capture_output=True)
                prots = [''.join(l.strip() for l in blk.split('\n')[1:]) for blk in open(tr).read().split('>')[1:]]
                pat = window_pattern(rnd.choice(prots)) if prots else None
            if pat is None:
                pat = ''.join(rnd.choice(AA + 'x') for _ in range(rnd.randint(1, 3)))
            out = os.path.join(tmp, 'o.gff')
            if os.path.exists(out): os.remove(out)
            r = subprocess.run(['fuzztran', '-sequence', fa, '-pattern', pat, '-frame', frame, '-table', str(table), '-outfile', out, '-rformat', 'gff', '-auto'], capture_output=True, text=True)
            if r.returncode != 0:
                sys.stderr.write('FAIL %s %s: %s\n' % (s, pat, r.stderr.strip()[:200])); continue
            if 'Illegal pattern' in (r.stdout + r.stderr):
                print('%s\t%s\t%s\t%d\tERR' % (s, pat, frame, table)); done += 1; continue
            hits = []
            if os.path.exists(out):
                for line in open(out):
                    if line.startswith('#') or not line.strip(): continue
                    f = line.rstrip('\n').split('\t')
                    note = f[8]
                    fr = int(re.search(r'\*frame (\d+)', note).group(1)); strand = f[6]
                    ps = int(re.search(r'\*start (\d+)', note).group(1)); pe = int(re.search(r'\*end (\d+)', note).group(1))
                    tl = re.search(r'\*translation (\S+)', note).group(1)
                    hits.append('%d:%d:%d:%s:%s:%s' % (fr if strand == '+' else -fr, ps, pe, tl, f[3], f[4]))
            hits.sort()
            print('%s\t%s\t%s\t%d\t%s' % (s, pat, frame, table, ' '.join(hits) or '-'))
            done += 1
    sys.stderr.write('rows: %d\n' % done)


if __name__ == '__main__':
    main()
