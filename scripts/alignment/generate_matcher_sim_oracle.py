#!/usr/bin/env python3
"""EMBOSS matcher (Huang-Miller SIM) oracle for SequenceAligner.LocalAlignTopN(..., TopNTieOrder.Sim) (audit A2-28d, B09 F145).

Every row is one `matcher -alternatives K` run (EMBOSS 6.6.0, installed): the full list of reported alignments with the program's own tie order.
Columns: kind(dna|prot) matrix gapopen gapextend seq1 seq2 K alignments
  alignments = 'ABORT' when the program aborts (fewer than K local alignments exist: "The number of alignments computed is too large"),
               else ';'-joined  score|start1|start2|end1|end2|row1|row2  (0-based start, 0-based inclusive end, as AlignmentResult).
Row sources: (1) every matcher row of topn_oracle.tsv (397; includes the 7 whose score list differs from Waterman-Eggert), (2) N new tie-rich cases
(tandem repeats, repeated domains, short alphabets, the sequence generators of generate_topn_oracle.py with a different seed and longer inputs).
Usage: python3 scripts/alignment/generate_matcher_sim_oracle.py [n_new=700] > matcher_sim_oracle.tsv
"""
import os, random, re, subprocess, sys, tempfile

here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, 'generate_topn_oracle.py')).read()
assert src.rstrip().endswith('main()')
exec(compile(src.rstrip()[:-len('main()')], 'generate_topn_oracle_nomain', 'exec'), globals())
rnd = random.Random(20261120)
TOPN = os.path.join(here, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'topn_oracle.tsv')


def run_matcher(a, b, k, G, E, matrix, kind, d):
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind == 'dna' else ['-sprotein1', '-sprotein2']
    r = subprocess.run(['matcher', d + '/a.fa', d + '/b.fa', '-gapopen', str(G), '-gapextend', str(E), '-datafile', matrix, '-alternatives', str(k),
                        '-aformat3', 'srspair', '-outfile', d + '/o.txt', '-auto'] + typ, capture_output=True, text=True)
    if r.returncode != 0:
        return 'ABORT'
    txt = open(d + '/o.txt').read()
    out = []
    for part in re.split(r'#=======================================\n#\n# Aligned_sequences', txt)[1:]:
        sc = int(re.search(r'# Score: (\d+)', part).group(1))
        ra = rb = ''; sa = sb = None
        for line in part.split('\n'):
            m_ = re.match(r'^a\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sa = sa or int(m_.group(1)); ra += m_.group(2); continue
            m_ = re.match(r'^b\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sb = sb or int(m_.group(1)); rb += m_.group(2)
        ea = sa - 1 + sum(1 for x in ra if x != '-'); eb = sb - 1 + sum(1 for x in rb if x != '-')
        out.append(f'{sc}|{sa - 1}|{sb - 1}|{ea - 1}|{eb - 1}|{ra}|{rb}')
    return ';'.join(out)


def main():
    n_new = int(sys.argv[1]) if len(sys.argv) > 1 else 700
    print('kind\tmatrix\tgapopen\tgapextend\tseq1\tseq2\tK\talignments')
    seen = set(); counts = {'rows': 0, 'abort': 0}
    with tempfile.TemporaryDirectory() as d:
        def emit(kind, G, E, a, b, k):
            key = (kind, G, E, a, b, k)
            if key in seen: return
            seen.add(key)
            res = run_matcher(a, b, k, G, E, 'EDNAFULL' if kind == 'dna' else 'EBLOSUM62', kind, d)
            counts['rows'] += 1; counts['abort'] += res == 'ABORT'
            print('\t'.join([kind, 'EDNAFULL' if kind == 'dna' else 'EBLOSUM62', str(G), str(E), a, b, str(k), res]))
        for f in [l.rstrip('\n').split('\t') for l in open(TOPN)][1:]:
            if f[0].startswith('#') or f[12] == 'NA': continue
            kind = 'dna' if f[0] == 'scalar' else 'prot'
            go, ge = float(f[5]), float(f[6])
            emit(kind, int(-(go + ge)), int(-ge), f[7], f[8], len(f[12].split(',')))
        for idx in range(n_new):
            if idx % 3 < 2:
                alpha = rnd.choice(['ACGT', 'ACGT', 'AC', 'ACGT'])
                a, b = gen_pair(alpha, rnd.choice([12, 18, 30, 45]), idx % 11 == 0)
                go, ge = rnd.choice(EMBOSS_DNA)[2:]
                go = rnd.choice([go, -12, -10, -6, -4, -8, -2]); ge = rnd.choice([ge, -1, -2, -4, -1])
                if go > ge: go = ge
                emit('dna', int(-go), int(-ge), a, b, rnd.randint(1, 10))
            else:
                a, b = gen_pair(rnd.choice([PROT, 'ACDEF', 'LIVM', 'AG']), rnd.choice([15, 25, 40]))
                go, ge = rnd.choice(EMBOSS_PROT_GAPS)
                emit('prot', int(-(go + ge)), int(-ge), a, b, rnd.randint(1, 8))
    print(f'# rows {counts["rows"]}, abort {counts["abort"]}', file=sys.stderr)


main()
