#!/usr/bin/env python3
"""ALIGN-LOCAL-001 (audit A2-28, B09 F122) oracle: EMBOSS 6.6.0 `water` (Smith-Waterman, affine gaps): the rows, the coordinates and the score of the program.

Seqeron API: SequenceAligner.WaterAlign(seq1, seq2, matrix, gapOpen, gapExtend).
Reference (run here): the installed `water -asequence a.fa -bsequence b.fa -auto -stdout -datafile M -gapopen G -gapextend E -aformat3 fasta` (rows), `-aformat3 markx10` (coordinates, sq_len / al_start / al_stop)
and `-aformat3 pair` (score) on random DNA (ACGT, AC, IUPAC) and protein pairs of 2 - 200 residues, related and unrelated, fractional and integer gap penalties (including non-dyadic ones).
Output TSV (tests/.../TestData/Alignment/water_oracle.tsv):  id  matrix  gapopen  gapextend  seq1  seq2  row1  row2  start1  start2  score      (start = 0-based first aligned residue; rows empty when water reports no alignment)
Usage: python3 scripts/alignment/generate_water_oracle.py OUTFILE [cases] [seed]      (needs EMBOSS water on PATH)
"""
import os, random, re, subprocess, sys, tempfile

PROT = ['EBLOSUM62', 'EBLOSUM62', 'EBLOSUM45', 'EBLOSUM80', 'EPAM250', 'EBLOSUM30', 'EPAM30']
DNA = ['EDNAFULL', 'EDNAFULL', 'EDNAMAT', 'EDNASIMPLE']
GAPOPEN = ['10', '10', '1', '2', '5', '8.5', '12', '15', '20', '50', '100', '3.7']
GAPEXT = ['0.5', '0.5', '0', '0.1', '0.3', '1', '2', '2.5', '5', '10', '0.25', '0.75', '1.5']


def mutate(rnd, s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate: out.append(rnd.choice(alpha))
        elif r < rate * 1.4: continue
        elif r < rate * 1.8: out.append(c); out.extend(rnd.choice(alpha) for _ in range(rnd.randint(1, 4)))
        else: out.append(c)
    return ''.join(out)


def main(outfile, cases, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        for i in range(cases):
            typ = 'n' if i % 2 == 0 else 'p'
            if typ == 'n':
                alpha = rnd.choice(['ACGT', 'ACGT', 'AC', 'ACGTRYN'])
            else:
                alpha = rnd.choice(['ACDEFGHIKLMNPQRSTVWY', 'ACDEFGHIKLMNPQRSTVWY', 'AG', 'ACDEFGHIKLMNPQRSTVWYBZX'])
            n = rnd.choice([2, 3, 8, 20, 45, 90, 150, 200])
            a = ''.join(rnd.choice(alpha) for _ in range(n))
            if rnd.random() < 0.85:
                core = a[rnd.randint(0, n // 3):]
                b = rnd.choice(alpha) * rnd.randint(0, 6) + mutate(rnd, core, alpha, rnd.choice([0.0, 0.05, 0.2, 0.5])) + ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 8)))
            else:
                b = ''.join(rnd.choice(alpha) for _ in range(rnd.choice([2, 5, 30, 120])))
            if len(b) < 2: b = b + rnd.choice(alpha)
            if rnd.random() < 0.2: a = a.lower()
            matrix = rnd.choice(DNA if typ == 'n' else PROT)
            go, ge = rnd.choice(GAPOPEN), rnd.choice(GAPEXT)
            open(os.path.join(d, 'a.fa'), 'w').write(f'>a\n{a}\n'); open(os.path.join(d, 'b.fa'), 'w').write(f'>b\n{b}\n')
            args = ['water', '-asequence', 'a.fa', '-bsequence', 'b.fa', '-auto', '-stdout', '-datafile', matrix, '-gapopen', go, '-gapextend', ge]
            f = subprocess.run(args + ['-aformat3', 'fasta'], cwd=d, capture_output=True, text=True, errors='replace')
            m = subprocess.run(args + ['-aformat3', 'markx10'], cwd=d, capture_output=True, text=True, errors='replace')
            p = subprocess.run(args + ['-aformat3', 'pair'], cwd=d, capture_output=True, text=True, errors='replace')
            recs = re.findall(r'>\w+\n([^>]*)', f.stdout)
            sc = re.search(r'^# Score: (-?[\d.]+)$', p.stdout, re.M)
            st = re.findall(r'; al_start: (\d+)', m.stdout)
            if f.returncode or len(recs) != 2 or not sc or len(st) != 2:
                continue
            r1, r2 = (x.replace('\n', '') for x in recs)
            if r1.replace('-', '').upper() not in a.upper() or r2.replace('-', '').upper() not in b.upper():
                continue    # the program converted symbols while reading (e.g. X -> N in a nucleotide-looking sequence)
            rows.append('\t'.join([str(len(rows)), matrix, go, ge, a, b, r1, r2, str(int(st[0]) - 1), str(int(st[1]) - 1), sc.group(1)]))
    with open(outfile, 'w') as fh:
        fh.write('id\tmatrix\tgapopen\tgapextend\tseq1\tseq2\trow1\trow2\tstart1\tstart2\tscore\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 900, int(sys.argv[3]) if len(sys.argv) > 3 else 20261123)
