#!/usr/bin/env python3
"""B09 audit A2-15: oracle for MultipleAlignmentStatistics.CalculateDistanceMatrix(rows, MsaDistanceOptions) = EMBOSS 6.6.0 `distmat` with every correction and option
(-nucmethod 0-5 / -protmethod 0-2, -gapweight, -ambiguous, -position 123/12/1/2/3, -calculatea, -parametera), on random gapped alignments (DNA with IUPAC codes, protein with X/B/Z/J).

    python3 scripts/alignment/generate_distmat_oracle.py [n=900] > tests/.../TestData/Alignment/distmat_oracle.tsv

Columns: rows(comma-joined)  type(dna|protein)  method  gapweight  position  ambiguous(0|1)  calculatea(0|1)  parametera  expected
  expected = the upper triangle of the printed matrix (x100, two decimals; 'nan' / '-nan' / 'inf' / '-inf' as the program prints them) joined by ',',
             or ERROR (the program exits non-zero: e.g. a position other than 1 2 3 12 123) or TIMEOUT (the program does not terminate within 4 s: -ambiguous on nucleotides with a gap).
"""
import random, re, subprocess, sys, tempfile

rnd = random.Random(20261119)
DNA = ['A', 'C', 'G', 'T'] * 6 + ['N', 'R', 'Y', 'M', 'K', 'S', 'W', 'B', 'D', 'H', 'V']
PROT = list('ACDEFGHIKLMNPQRSTVWY') * 2 + ['X', 'B', 'Z', 'J']


def make_msa(kind, gaps=True, ambig=True):
    L = rnd.randint(9, 60)
    n = rnd.randint(2, 7)
    alpha = (DNA if ambig else list('ACGT')) if kind == 'dna' else (PROT if ambig else list('ACDEFGHIKLMNPQRSTVWY'))
    base = [rnd.choice(alpha) for _ in range(L)]
    rows = []
    for _ in range(n):
        rate = rnd.choice([0.0, 0.05, 0.15, 0.3, 0.5, 0.7])
        row = [c if rnd.random() > rate else rnd.choice(alpha) for c in base]
        if gaps:
            for _ in range(rnd.randint(0, 3)):
                s = rnd.randrange(L); e = min(L, s + rnd.randint(1, 5))
                for k in range(s, e): row[k] = '-'
        rows.append(row)
    for k in range(L):
        if all(r[k] == '-' for r in rows):
            rows[rnd.randrange(n)][k] = rnd.choice(alpha)
    rows = [''.join(r) for r in rows]
    if kind == 'dna' and rnd.random() < 0.5:     # EMBOSS decides the type from the first row only: keep it clean DNA (ACGT + ambiguity + gaps)
        rows[0] = re.sub('[^ACGTNRYMKSWBDHV-]', 'A', rows[0])
    if any(set(r) == {'-'} for r in rows): return make_msa(kind, gaps, ambig)
    return rows


def run(d, rows, kind, method, gapw, pos, amb, calca, parama):
    open(d + '/a.fa', 'w').write(''.join(f'>s{i + 1}\n{r}\n' for i, r in enumerate(rows)))
    cmd = ['distmat', '-sequence', d + '/a.fa', '-outfile', d + '/o.txt', '-auto', '-gapweight', str(gapw), '-position', str(pos),
           '-nucmethod', str(method), '-protmethod', str(method if method <= 2 else 0), '-parametera', str(parama)]    # EMBOSS decides the type from the first row: give both menus the method
    if amb: cmd.append('-ambiguous')
    if calca: cmd.append('-calculatea')
    try:
        p = subprocess.run(cmd, capture_output=True, timeout=4)
    except subprocess.TimeoutExpired:
        return 'TIMEOUT'
    if p.returncode != 0:
        return 'ERROR'
    n = len(rows)
    lines = open(d + '/o.txt').read().split('\n')
    start = [i for i, l in enumerate(lines) if re.match(r'^\s+1\s+2', l)][0] + 1
    vals = []
    for j in range(n - 1):
        f = lines[start + j].split('\t')
        cols = f[1:n + 1]
        vals += [cols[i].strip() for i in range(j + 1, n)]
    return ','.join(vals)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 900
    print('rows\ttype\tmethod\tgapweight\tposition\tambiguous\tcalculatea\tparametera\texpected')
    with tempfile.TemporaryDirectory() as d:
        for t in range(n):
            kind = rnd.choice(['dna', 'dna', 'protein'])
            method = rnd.choice([0, 1, 2, 3, 4, 5]) if kind == 'dna' else rnd.choice([0, 1, 2])
            amb = rnd.random() < (0.2 if kind == 'dna' else 0.3) and method in (0, 1, 4)
            gaps = not (amb and kind == 'dna' and rnd.random() < 0.8)    # -ambiguous with a gap in the first sequence of a pair does not terminate
            rows = make_msa(kind, gaps)
            gapw = rnd.choice([0, 0, 0.5, 1, 2.5])
            pos = rnd.choice([123, 123, 123, 12, 1, 2, 3]) if kind == 'dna' else rnd.choice([123, 123, 123, 12])
            if rnd.random() < 0.03: pos = rnd.choice([23, 13, 4, 0]) if kind == 'dna' else rnd.choice([1, 2, 3, 23])
            calca = method == 5 and rnd.random() < 0.5
            parama = rnd.choice([1.0, 1.0, 0.5, 2.0, 3.5]) if method == 5 else 1.0
            exp = run(d, rows, kind, method, gapw, pos, amb, calca, parama)
            print('\t'.join([','.join(rows), kind, str(method), repr(float(gapw)), str(pos), '1' if amb else '0', '1' if calca else '0', repr(float(parama)), exp]))


if __name__ == '__main__':
    main()
