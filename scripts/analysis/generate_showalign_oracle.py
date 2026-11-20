#!/usr/bin/env python3
"""ALIGN-STATS-001 oracle (audit A3-13c): EMBOSS showalign 6.6.0 (pretty display of a multiple alignment) vs MultipleAlignmentDisplay.ShowAlign.

Random alignments (1-8 rows, 4-120 columns; nucleotide with ambiguity letters or protein; upper / lower case letters, gaps - . ~, ragged rows) and random options:
-refseq (0 / number / name), -[no]bottom, -show A|I|N|S|D, -order I|A|S, -[no]similarcase, -[no]consensus, -[no]number, -[no]ruler, -width, -margin, -uppercase ranges,
-matrix, -plurality, -setcase, -identity, -[no]gaps. Output = the text file of the program (never -html / -highlight: not ported).
JSON lines: {"names": [...], "rows": [...], "args": {...program options...}, "out": "<text>"}
Never aborts: a failing showalign call is reported on stderr and the case skipped.
Usage: python3 scripts/analysis/generate_showalign_oracle.py [n=500] > showalign_oracle.jsonl
"""
import json, os, random, subprocess, sys, tempfile

rnd = random.Random(20261126)
NUC = 'ACGT'; AMB = 'RYSWKMBDHVN'; AA = 'ACDEFGHIKLMNPQRSTVWY'


def row(n, kind, base):
    out = []
    for i in range(n):
        r = rnd.random()
        if r < 0.12: c = rnd.choice('-.~') if rnd.random() < 0.3 else '-'
        elif base and r < 0.82: c = base[i]
        elif kind == 'nuc': c = rnd.choice(AMB) if rnd.random() < 0.12 else rnd.choice(NUC)
        else: c = rnd.choice('BZX') if rnd.random() < 0.05 else rnd.choice(rnd.choice([AA, 'ACDEFG', 'LIVM']))
        out.append(c)
    return ''.join(out)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    rows_out = bad = 0
    with tempfile.TemporaryDirectory() as d:
        for _ in range(n):
            kind = 'prot' if rnd.random() < 0.4 else 'nuc'
            nseq = rnd.randint(1, 8) if rnd.random() < 0.9 else rnd.randint(2, 4)
            mlen = rnd.choice([4, 9, 10, 11, 20, 37, 60, 61, 95, 120])
            alpha = AA if kind == 'prot' else NUC
            base = ''.join(rnd.choice(alpha) for _ in range(mlen))
            seqs = [row(mlen, kind, base) for _ in range(nseq)]
            if rnd.random() < 0.15 and nseq > 1:
                i = rnd.randrange(nseq); seqs[i] = seqs[i][:rnd.randint(1, mlen)]
            if rnd.random() < 0.3:
                seqs = [s.lower() if rnd.random() < 0.5 else s for s in seqs]
            names = ['seq%s' % rnd.choice(['A', 'b', 'Cc', '_x', 'LongerName%d' % j]) + str(j) for j in range(nseq)]
            args = {}
            ref = rnd.choice(['0', '0', '0', str(rnd.randint(1, nseq)), names[rnd.randrange(nseq)]])
            args['refseq'] = ref
            if rnd.random() < 0.5: args['bottom'] = rnd.choice(['Y', 'N'])
            args['show'] = rnd.choice('ANNIISD')
            args['order'] = rnd.choice('IIAS')
            if rnd.random() < 0.4: args['similarcase'] = rnd.choice(['Y', 'N'])
            if rnd.random() < 0.3: args['consensus'] = rnd.choice(['Y', 'N'])
            if rnd.random() < 0.3: args['number'] = rnd.choice(['Y', 'N'])
            if rnd.random() < 0.3: args['ruler'] = rnd.choice(['Y', 'N'])
            args['width'] = rnd.choice([60, 60, 10, 7, 25, 1, 40, 80])
            args['margin'] = rnd.choice([-1, -1, -1, 0, 3, 12])
            if rnd.random() < 0.3:
                a = rnd.randint(1, mlen); b = rnd.randint(a, mlen); args['uppercase'] = '%d-%d' % (a, b)
            if rnd.random() < 0.25: args['plurality'] = rnd.choice([0, 25, 50, 75, 100])
            if rnd.random() < 0.15: args['setcase'] = rnd.choice([0, 1, 2, 3])
            if rnd.random() < 0.15: args['identity'] = rnd.choice([30, 60, 100])
            if rnd.random() < 0.3: args['gaps'] = rnd.choice(['Y', 'N'])
            if rnd.random() < 0.15: args['matrix'] = rnd.choice(['EBLOSUM62', 'EPAM250', 'EDNAFULL', 'EBLOSUM80'] if kind == 'prot' else ['EDNAFULL', 'EDNASIMPLE'])
            f = os.path.join(d, 's.fa'); o = os.path.join(d, 'o.txt')
            open(f, 'w').write(''.join('>%s\n%s\n' % (nm, s) for nm, s in zip(names, seqs)))
            cmd = ['showalign', f, '-outfile', o, '-auto']
            for k, v in args.items():
                if k in ('bottom', 'similarcase', 'consensus', 'number', 'ruler', 'gaps'):
                    cmd.append(('-' if v == 'Y' else '-no') + k)
                else:
                    cmd += ['-' + k, str(v)]
            r = subprocess.run(cmd, capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(o):
                sys.stderr.write('REJECTED %s: %s\n' % (cmd[1:], r.stderr.strip()[-100:])); bad += 1
                continue
            out = open(o, encoding='latin-1').read(); os.remove(o)
            print(json.dumps({'names': names, 'rows': seqs, 'args': args, 'out': out}))
            rows_out += 1
    sys.stderr.write('rows=%d rejected=%d\n' % (rows_out, bad))


main()
