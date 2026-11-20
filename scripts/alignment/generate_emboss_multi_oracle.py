#!/usr/bin/env python3
"""ALIGN-LOCAL-001 / ALIGN-GLOBAL-001 (audit A2-28c, B09 F124) oracle: EMBOSS 6.6.0 `supermatcher` (Rabin-Karp word seeds + banded Smith-Waterman) and `needleall` (all-vs-all needle), whole output files.

Seqeron API: SequenceAligner.Supermatcher(queries, targets, ...) / SequenceAligner.NeedleAll(setA, setB, ...) formatted with EmbossAlignmentWriter.
Reference (run here): the installed programs with `-auto -stdout -errfile err.txt -aformat <fmt>` on random related / unrelated DNA, IUPAC-DNA and protein sets (1-3 query and 1-4 target sequences, 20-600 residues, word length,
band width, minimum score, matrix, gap penalties, end-gap options, upper / lower case, N / X runs).
Output TSV (tests/.../TestData/Alignment/emboss_multi_oracle.tsv):  id  program  args  fastaA  fastaB  format  date  file  errfile
  args = the command line words after the program name (space separated; file names a.fa / b.fa), date = the Rundate of the output as ISO, file / errfile = the output (\\n \\t \\\\ escaped).
Usage: python3 scripts/alignment/generate_emboss_multi_oracle.py OUTFILE [casesPerProgram] [seed]      (needs EMBOSS supermatcher and needleall on PATH)
"""
import datetime, os, random, re, subprocess, sys, tempfile

PROT = ['EBLOSUM62', 'EBLOSUM62', 'EBLOSUM45', 'EBLOSUM80', 'EPAM250', 'EBLOSUM30']
DNA = ['EDNAFULL', 'EDNAFULL', 'EDNAMAT', 'EDNASIMPLE']
FORMATS = ['pair', 'srspair', 'fasta', 'markx0', 'markx3', 'markx10', 'simple', 'score', 'match', 'sam', 'tcoffee', 'msf']


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')


def mutate(rnd, s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate: out.append(rnd.choice(alpha))
        elif r < rate * 1.4: continue
        elif r < rate * 1.8: out.append(c); out.extend(rnd.choice(alpha) for _ in range(rnd.randint(1, 6)))
        else: out.append(c)
    return ''.join(out)


def fasta(names, seqs, rnd):
    out = []
    for n, s in zip(names, seqs):
        d = '' if rnd.random() < 0.6 else ' ' + ' '.join(rnd.choice(['alpha', 'beta', 'gene', 'kinase', 'x1']) for _ in range(rnd.randint(1, 3)))
        out.append('>' + n + d + '\n' + '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
    return ''.join(out)


def make_case(rnd, program):
    typ = 'n' if rnd.random() < 0.55 else 'p'
    alpha = rnd.choice(['ACGT', 'ACGT', 'AC', 'ACGTRYN']) if typ == 'n' else rnd.choice(['ACDEFGHIKLMNPQRSTVWY'] * 3 + ['ACDEFGHIKLMNPQRSTVWYBZX'])
    base_len = rnd.choice([20, 40, 80, 150, 300, 600] if program == 'supermatcher' else [10, 30, 60, 120])
    base = ''.join(rnd.choice(alpha) for _ in range(base_len))
    if rnd.random() < 0.25:                      # runs of N / X
        k = rnd.randint(0, max(0, len(base) - 10)); base = base[:k] + ('N' if typ == 'n' else 'X') * rnd.randint(1, 12) + base[k + 12:]
    nq, nt = rnd.choice([1, 1, 1, 2, 3]), rnd.choice([1, 1, 2, 3, 4])
    def derive():
        s = mutate(rnd, base[rnd.randint(0, base_len // 3):len(base) - rnd.randint(0, base_len // 3)], alpha, rnd.choice([0.0, 0.03, 0.1, 0.25]))
        if rnd.random() < 0.2: s = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(8, 60))) + s
        if rnd.random() < 0.15: s = ''.join(rnd.choice(alpha) for _ in range(rnd.choice([4, 15, 50])))
        if rnd.random() < 0.15: s = s.lower()
        return s or alpha[0]
    qs = [derive() for _ in range(nq)]
    ts = [derive() for _ in range(nt)]
    return typ, qs, ts


def run(prog, args, d):
    return subprocess.run([prog] + args, cwd=d, capture_output=True, text=True, errors='replace')


def main(outfile, per, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        for program in ['supermatcher', 'needleall']:
            made = 0
            attempts = 0
            while made < per and attempts < per * 4:
                attempts += 1
                typ, qs, ts = make_case(rnd, program)
                if any(set(x.upper()) <= set('ACGTUBDHKMNRSVWXY') and 'X' in x.upper() for x in qs + ts):
                    continue            # EMBOSS reads such a sequence as nucleotide and converts X to N
                qn = ['q%d' % (i + 1) for i in range(len(qs))]; tn = ['t%d' % (i + 1) for i in range(len(ts))]
                a_seqs, b_seqs = (qs, ts)
                open(os.path.join(d, 'a.fa'), 'w').write(fasta(qn, qs, rnd)); open(os.path.join(d, 'b.fa'), 'w').write(fasta(tn, ts, rnd))
                args = ['-asequence', 'a.fa', '-bsequence', 'b.fa', '-auto', '-stdout', '-errfile', 'err.txt']
                if rnd.random() < 0.8: args += ['-datafile', rnd.choice(DNA if typ == 'n' else PROT)]
                if rnd.random() < 0.7: args += ['-gapopen', rnd.choice(['10', '1', '3', '5', '8.5', '12', '20', '50', '3.7']) if program == 'needleall' else rnd.choice(['10', '0', '3', '5', '8.5', '12', '20', '50', '3.7'])]
                if rnd.random() < 0.7: args += ['-gapextend', rnd.choice(['0.5', '0', '0.1', '1', '2', '5', '0.25', '1.5'])]
                if program == 'supermatcher':
                    if rnd.random() < 0.7: args += ['-wordlen', str(rnd.choice([3, 3, 4, 5, 6, 6, 7, 8, 10, 15]))]
                    if rnd.random() < 0.7: args += ['-width', str(rnd.choice([1, 2, 3, 4, 7, 8, 16, 16, 20, 33, 50, 100]))]
                    if rnd.random() < 0.4: args += ['-minscore', str(rnd.choice([0, 5, 20, 50, 100, 300]))]
                else:
                    if rnd.random() < 0.5: args += ['-endweight', '-endopen', rnd.choice(['10', '1', '5', '20']), '-endextend', rnd.choice(['0.5', '0', '1', '2'])]
                    if rnd.random() < 0.5: args += ['-minscore', str(rnd.choice([-10, 0, 1, 5, 20, 50, 100]))]
                    if rnd.random() < 0.4: args += ['-nobrief']
                fmt = rnd.choice(FORMATS)
                args += ['-aformat', fmt]
                for f in ['err.txt']:
                    if os.path.exists(os.path.join(d, f)): os.remove(os.path.join(d, f))
                p = run(program, args, d)
                if p.returncode or 'Error' in p.stderr or 'Fatal' in p.stderr or 'Warning' in p.stderr:
                    continue
                m = re.search(r'^# Rundate: (.*)$', p.stdout, re.M)
                if not m: continue
                date = datetime.datetime.strptime(m.group(1).strip(), '%a %d %b %Y %H:%M:%S').strftime('%Y-%m-%dT%H:%M:%S')
                err = open(os.path.join(d, 'err.txt')).read() if os.path.exists(os.path.join(d, 'err.txt')) else ''
                rows.append('\t'.join([str(len(rows)), program, ' '.join(args), esc(open(os.path.join(d, 'a.fa')).read()), esc(open(os.path.join(d, 'b.fa')).read()), fmt, date, esc(p.stdout), esc(err)]))
                made += 1
    with open(outfile, 'w') as fh:
        fh.write('id\tprogram\targs\tfastaA\tfastaB\tformat\tdate\tfile\terrfile\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 200, int(sys.argv[3]) if len(sys.argv) > 3 else 20261125)
