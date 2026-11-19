#!/usr/bin/env python3
"""ALIGN-GLOBAL-001 (audit A2-28, B09 F121) oracle: EMBOSS 6.6.0 `stretcher` (Myers-Miller linear-space global alignment): the rows and the score of the program.

Seqeron API: SequenceAligner.StretcherAlign(seq1, seq2, matrix, gapOpen, gapExtend).
Reference (run here): the installed `stretcher -asequence a.fa -bsequence b.fa -datafile M -gapopen G -gapextend E -aformat3 fasta` (rows) and `-aformat3 pair` (score).  Random DNA (ACGT, low-complexity AC,
IUPAC) and protein pairs of 1 - 260 residues, related and unrelated, integer gap penalties 0 - 20 / 0 - 6, the nucleotide and protein matrices of EMBOSS.
Output TSV (tests/.../TestData/Alignment/stretcher_oracle.tsv):  id  matrix  gapopen  gapextend  seq1  seq2  row1  row2  score
Usage: python3 scripts/alignment/generate_stretcher_oracle.py OUTFILE [cases] [seed]      (needs EMBOSS stretcher on PATH)
"""
import os, random, re, subprocess, sys, tempfile

PROT = ['EBLOSUM62', 'EBLOSUM62', 'EBLOSUM45', 'EBLOSUM80', 'EPAM250', 'EBLOSUM30', 'EPAM30']
DNA = ['EDNAFULL', 'EDNAFULL', 'EDNAMAT', 'EDNASIMPLE']


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
            n = rnd.choice([1, 2, 3, 8, 20, 45, 90, 150, 260])
            a = ''.join(rnd.choice(alpha) for _ in range(n))
            b = mutate(rnd, a, alpha, rnd.choice([0.0, 0.05, 0.2, 0.5])) if rnd.random() < 0.85 else ''.join(rnd.choice(alpha) for _ in range(rnd.choice([1, 5, 30, 120])))
            if not b: b = rnd.choice(alpha)
            if rnd.random() < 0.2: a = a.lower()
            matrix = rnd.choice(DNA if typ == 'n' else PROT)
            go = rnd.choice([0, 1, 4, 8, 10, 12, 16, 20]); ge = rnd.choice([0, 1, 2, 4, 6])
            open(os.path.join(d, 'a.fa'), 'w').write(f'>a\n{a}\n'); open(os.path.join(d, 'b.fa'), 'w').write(f'>b\n{b}\n')
            args = ['stretcher', '-asequence', 'a.fa', '-bsequence', 'b.fa', '-auto', '-stdout', '-datafile', matrix, '-gapopen', str(go), '-gapextend', str(ge)]
            f = subprocess.run(args + ['-aformat3', 'fasta'], cwd=d, capture_output=True, text=True, errors='replace')
            p = subprocess.run(args + ['-aformat3', 'pair'], cwd=d, capture_output=True, text=True, errors='replace')
            recs = re.findall(r'>\w+\n([^>]*)', f.stdout)
            m = re.search(r'^# Score: (-?\d+)$', p.stdout, re.M)
            if f.returncode or len(recs) != 2 or not m:
                continue
            r1, r2 = (x.replace('\n', '') for x in recs)
            if r1.replace('-', '').upper() != a.upper() or r2.replace('-', '').upper() != b.upper():
                continue    # the program converted symbols while reading (e.g. X -> N in a nucleotide-looking sequence)
            rows.append('\t'.join([str(len(rows)), matrix, str(go), str(ge), a, b, r1, r2, m.group(1)]))
    with open(outfile, 'w') as fh:
        fh.write('id\tmatrix\tgapopen\tgapextend\tseq1\tseq2\trow1\trow2\tscore\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 700, int(sys.argv[3]) if len(sys.argv) > 3 else 20261121)
