#!/usr/bin/env python3
"""EMBOSS matcher `-gapopen` / `-gapextend` / `-alternatives` oracle (audit A6-3, B09 F178): SequenceAligner.MatcherProgram on every integer pair the program accepts, in particular -gapopen BELOW -gapextend
(`matcher a b -gapopen 0 -gapextend 4`: the first gap residue costs less than a further one, SIM's q = gapopen - gapextend is negative), 0, 1 and 100, and on the values the ACD resets or refuses.
Columns: kind(dna|prot) gapopen gapextend alternatives seq1 seq2 status alignments
  gapopen / gapextend / alternatives are the TEXTS given on the command line;  status = 'OK' | 'DIED' (the program died: invalid integer / fewer alignments than -alternatives) ;
  alignments = ';'-joined score|start1|start2|end1|end2|row1|row2 (0-based start, 0-based inclusive end) as in matcher_sim_oracle.tsv, then '~' and the program's stderr warning lines joined by '|'.
Usage: python3 scripts/alignment/generate_matcher_gap_oracle.py [n=800] > matcher_gap_oracle.tsv      (needs EMBOSS 6.6.0 matcher on PATH)
"""
import os, random, re, subprocess, sys, tempfile

here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, 'generate_topn_oracle.py')).read()
assert src.rstrip().endswith('main()')
exec(compile(src.rstrip()[:-len('main()')], 'generate_topn_oracle_nomain', 'exec'), globals())
rnd = random.Random(20261222)


def run(a, b, G, E, K, kind, d):
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind == 'dna' else ['-sprotein1', '-sprotein2']
    matrix = 'EDNAFULL' if kind == 'dna' else 'EBLOSUM62'
    if os.path.exists(d + '/o.txt'): os.remove(d + '/o.txt')
    r = subprocess.run(['matcher', d + '/a.fa', d + '/b.fa', '-gapopen', G, '-gapextend', E, '-datafile', matrix, '-alternatives', K, '-aformat3', 'srspair', '-outfile', d + '/o.txt', '-auto'] + typ,
                       capture_output=True, text=True)
    warn = '|'.join(l for l in r.stderr.split('\n') if l.startswith('Warning: integer value out of range'))
    if r.returncode != 0 or not os.path.exists(d + '/o.txt'):
        return 'DIED', ''
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
    return 'OK', ';'.join(out) + '~' + warn


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 800
    print('kind\tgapopen\tgapextend\talternatives\tseq1\tseq2\tstatus\talignments')
    seen = set(); counts = {}
    odd = ['-1', '-5', '1.5', '0.5', 'abc', '1e2', '+3', '007', '', '2147483647', '99999999999']
    with tempfile.TemporaryDirectory() as d:
        for idx in range(n):
            kind = 'dna' if idx % 2 == 0 else 'prot'
            alpha = rnd.choice(['ACGT', 'ACGT', 'AC']) if kind == 'dna' else rnd.choice([PROT, 'ACDEF', 'LIVM'])
            a, b = gen_pair(alpha, rnd.choice([12, 18, 30, 45]), idx % 11 == 0)
            if idx % 5 == 4:
                G = rnd.choice(odd); E = rnd.choice(['4', '2', '0', '1', '-3', '1.5'] + odd[:6]); K = rnd.choice(['1', '1', '2', '0', '-2'])
            else:
                G = str(rnd.choice([0, 1, 2, 3, 4, 5, 8, 10, 14, 16, 20, 100, 0, 1, 2, 3])); E = str(rnd.choice([0, 1, 2, 3, 4, 4, 5, 8, 10])); K = str(rnd.choice([1, 1, 2, 3]))
            if G in ('2147483647', '99999999999') or E in ('2147483647', '99999999999'): continue        # the program overflows its 32-bit scores (undefined results)
            key = (kind, G, E, K, a, b)
            if key in seen: continue
            seen.add(key)
            status, res = run(a, b, G, E, K, kind, d)
            counts[status] = counts.get(status, 0) + 1
            print('\t'.join([kind, G, E, K, a, b, status, res]))
    print('# rows', counts, file=sys.stderr)


main()
