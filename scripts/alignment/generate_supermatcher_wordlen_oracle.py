#!/usr/bin/env python3
"""EMBOSS supermatcher `-wordlen` above the targets (audit A6-5, B09 F180): `supermatcher -wordlen W` with W around and above the length of every target.
The program builds its word table from the targets (embWordGetTable); when the table is empty (W above every target, or only N / X runs) it prints `Error: no kmers found` on stderr and goes on (no alignment; the
error file lists "No wordmatch start points for <query> vs <target>. No alignment" for every pair; exit status 0 up to about W = 100 000, a segmentation fault from about 1 000 000 on, which this oracle does not use).
TSV: id  typ  fastaA  fastaB  wordlen  stderr  errfile  alignments    (fasta / stderr / errfile escaped as in emboss_multi_oracle.tsv; alignments = ';'-joined  query|target|score  of the -aformat pair file)
Usage: python3 scripts/alignment/generate_supermatcher_wordlen_oracle.py OUTFILE [n=240] [seed]      (needs EMBOSS supermatcher on PATH)
"""
import os, random, re, sys, tempfile

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import generate_emboss_multi_oracle as base  # noqa: E402


def main(outfile, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            typ, qs, ts = base.make_case(rnd, 'supermatcher')
            if any(set(x.upper()) <= set('ACGTUBDHKMNRSVWXY') and 'X' in x.upper() for x in qs + ts): continue
            longest = max(len(t) for t in ts)
            c = rnd.random()
            if c < 0.45: w = rnd.randint(3, longest + 3)
            elif c < 0.65: w = longest + rnd.randint(1, 40)
            elif c < 0.80: w = rnd.choice([longest, longest - 1, longest + 1])
            else: w = rnd.choice([100, 1000, 100000, 999999, longest * 3 + 7])
            w = max(w, 3)
            qn = ['q%d' % (i + 1) for i in range(len(qs))]; tn = ['t%d' % (i + 1) for i in range(len(ts))]
            fa = base.fasta(qn, qs, rnd); fb = base.fasta(tn, ts, rnd)
            open(os.path.join(d, 'a.fa'), 'w').write(fa); open(os.path.join(d, 'b.fa'), 'w').write(fb)
            err = os.path.join(d, 'err.txt')
            if os.path.exists(err): os.remove(err)
            p = base.run('supermatcher', ['-asequence', 'a.fa', '-bsequence', 'b.fa', '-wordlen', str(w), '-auto', '-stdout', '-errfile', 'err.txt', '-aformat', 'pair'], d)
            if p.returncode != 0: continue
            errfile = open(err).read() if os.path.exists(err) else ''
            aligns = []
            for part in re.split(r'#=======================================\n#\n# Aligned_sequences', p.stdout)[1:]:
                m1 = re.search(r'# 1: (\S+)', part); m2 = re.search(r'# 2: (\S+)', part); sc = re.search(r'# Score: (\S+)', part)
                aligns.append('%s|%s|%s' % (m1.group(1), m2.group(1), sc.group(1)))
            rows.append('\t'.join([str(len(rows)), typ, base.esc(fa), base.esc(fb), str(w), base.esc(p.stderr), base.esc(errfile), ';'.join(aligns)]))
    with open(outfile, 'w', newline='\n') as fh:
        fh.write('id\ttyp\tfastaA\tfastaB\twordlen\tstderr\terrfile\talignments\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if 'no kmers' in r), 'with "no kmers found"')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 240, int(sys.argv[3]) if len(sys.argv) > 3 else 20261223)
