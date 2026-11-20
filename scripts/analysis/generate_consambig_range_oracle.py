#!/usr/bin/env python3
"""MOTIF-CONS-001 oracle (audit A4-9): EMBOSS consambig 6.6.0 with -sbegin / -send vs MotifFinder.GenerateEmbossAmbiguousConsensus (which has no range options).

Measured on the installed program: consambig IGNORES -sbegin and -send (the consensus covers every column of the alignment whatever the range, positive, negative or beyond the alignment). The same random alignments as
generate_consambig_oracle.py (A3-13d), every row with a random -sbegin and -send; the recorded consensus is the program's output WITH the range, which the test compares with the plain call of the port.

Random alignments (2-8 rows, 1-30 columns; ragged rows included - the program pads them with trailing gaps): nucleotides (A C G T, the 11 IUPAC ambiguity letters, U, X, ?, gaps - . ~, lower case)
and proteins (20 amino acids, B Z J X U O, *, gaps, lower case); type = auto (decided from the first row by the program), or forced with -snucleotide / -sprotein.
TSV: type(auto|nuc|prot)  rows(comma-separated)  (empty)  sbegin  send  consensus   (the NUL character the program writes for a column without a scored amino acid is written '@')
Never aborts: a failing consambig call is reported on stderr and the case skipped.
Usage: python3 scripts/analysis/generate_consambig_range_oracle.py [n=300] > consambig_range_oracle.tsv
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261204)
NUC = 'ACGT'; AMB = 'RYSWKMBDHVN'; GAPS = '-.~'
AA = 'ACDEFGHIKLMNPQRSTVWY'


def nuc_row(n, rate_amb, rate_gap, rate_odd):
    out = []
    for _ in range(n):
        r = rnd.random()
        if r < rate_gap: c = rnd.choice(GAPS)
        elif r < rate_gap + rate_amb: c = rnd.choice(AMB)
        elif r < rate_gap + rate_amb + rate_odd: c = rnd.choice('UX?')
        else: c = rnd.choice(NUC)
        out.append(c.lower() if rnd.random() < 0.08 else c)
    return ''.join(out)


def prot_row(n, rate_amb, rate_gap, rate_odd):
    out = []
    for _ in range(n):
        r = rnd.random()
        if r < rate_gap: c = rnd.choice(GAPS)
        elif r < rate_gap + rate_amb: c = rnd.choice('BZJXUO')
        elif r < rate_gap + rate_amb + rate_odd: c = rnd.choice('*?')
        else: c = rnd.choice(rnd.choice([AA, 'ACDE', 'DNEQIL']))
        out.append(c.lower() if rnd.random() < 0.08 else c)
    return ''.join(out)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    rows_out = bad = 0
    with tempfile.TemporaryDirectory() as d:
        for _ in range(n):
            protein = rnd.random() < 0.4
            nseq = rnd.randint(2, 8); mlen = rnd.randint(1, 30)
            amb, gap, odd = rnd.choice([(0.0, 0.0, 0.0), (0.1, 0.05, 0.0), (0.25, 0.1, 0.03), (0.05, 0.3, 0.02)])
            gen = prot_row if protein else nuc_row
            seqs = [gen(mlen, amb, gap, odd) for _ in range(nseq)]
            if rnd.random() < 0.12:   # ragged
                i = rnd.randrange(nseq); seqs[i] = seqs[i][:rnd.randint(1, mlen)]
            if rnd.random() < 0.05:   # a row of the other kind
                i = rnd.randrange(nseq); seqs[i] = (nuc_row if protein else prot_row)(mlen, 0.0, 0.0, 0.0)
            typ = rnd.choice(['auto', 'auto', 'prot' if protein else 'nuc'])
            f = os.path.join(d, 's.fa')
            open(f, 'w').write(''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs)))
            args = ['consambig', f, '-outseq', os.path.join(d, 'o.fa'), '-osformat', 'fasta', '-auto']
            sb = rnd.choice([0, 1, 2, 3, mlen, mlen + 4, -1, -2, -mlen, -(mlen + 3), rnd.randint(1, mlen)]); se = rnd.choice([0, 1, 2, mlen, mlen + 4, -1, -3, -(mlen + 3), rnd.randint(1, mlen)])
            args += ['-sbegin', str(sb), '-send', str(se)]
            if typ == 'nuc': args.append('-snucleotide')
            if typ == 'prot': args.append('-sprotein')
            r = subprocess.run(args, capture_output=True, text=True)
            if r.returncode != 0 or not os.path.exists(os.path.join(d, 'o.fa')):
                sys.stderr.write('REJECTED %s %s: %s\n' % (typ, seqs, r.stderr.strip()[-120:])); bad += 1; continue
            raw = open(os.path.join(d, 'o.fa'), 'rb').read().decode('latin-1').split('\n')
            cons = ''.join(l for l in raw[1:] if not l.startswith('>')).replace('\x00', '@')
            os.remove(os.path.join(d, 'o.fa'))
            print('%s\t%s\t%s\t%d\t%d\t%s' % (typ, ','.join(seqs), '', sb, se, cons)); rows_out += 1
    sys.stderr.write('rows=%d rejected=%d\n' % (rows_out, bad))


main()
