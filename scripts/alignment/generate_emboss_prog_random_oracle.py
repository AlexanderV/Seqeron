#!/usr/bin/env python3
"""B09 WP62 (F281) oracle of the random sequence programs of EMBOSS 6.6.0: shuffleseq, msbar, makenucseq, makeprotseq.
ajRandomSeed takes the microseconds of gettimeofday as the seed (tv_usec % 9999 + 1) unless the EMBOSS variable timetoday is defined (then 0), so a run is deterministic when the clock is fixed:
every row fixes it by the LD_PRELOAD shim faketime.c (FAKE_USEC, row field 'usec') or by `SET timetoday 1` in ~/.embossrc (row field 'timetoday'). Every row is run twice (a row whose two runs differ is dropped) and the
rows that can touch memory outside a string (very short sequences, deletions) are run under valgrind: a reference run with an invalid access has no defined result and is stored as a refusal (rc -11).
Usage: python3 scripts/alignment/generate_emboss_prog_random_oracle.py PROGRAM OUT.jsonl.gz [n=400] [seed]      (needs EMBOSS 6.6.0)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402
import generate_emboss_prog_seq_oracle as G  # noqa: E402


def clock(rnd):
    return {'timetoday': rnd.choice(['2012-03-04 05:06:07', '1999-12-31 23:59:59', '2030-01-01 00:00:00', '1', '2012-03-04', '12-3-4 5:6:7', '2012-13-45 25:61:61', '2100-02-29 12:00:00', 'abc'])} if rnd.random() < 0.12 else {'usec': rnd.choice([0, 1, 7, 99, 4242, 9998, 9999, 10000, 500000, 999999, rnd.randint(0, 999999), rnd.randint(0, 999999)])}


def gen_shuffleseq(rnd, i, d):
    kind = rnd.choice(['n', 'n', 'p'])
    text, fmt, args = G.common(rnd, d, 'shuffleseq', kind, nseq=rnd.randint(1, 4), lens=(1, 120), gaps=rnd.random() < 0.1)
    if rnd.random() < 0.6: args += ['-shuffle', str(rnd.choice([1, 2, 3, 5, 0, -1, 10]))]
    if rnd.random() < 0.03: args += ['-shuffle', rnd.choice(['x', '1.5', ''])]
    ex = clock(rnd)
    return args + G.out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}, ex


LIST_VALUES = ['0', '1', '2', '3', '4', '5', '6', '2,3', '3,4,5', '2,3,4,5', '6,5', '1', '*', 'Insertions', 'Del', 'None', 'Any', 'ins,del', 'Moves,Changes', '2,3,4,5,6', '7', 'x', '', 'Ch', 'D', '2,,3', ',2', '2,', '0,2', '2,0']


def msb_lists(rnd, args):
    for q in ('-point', '-block', '-codon'):
        if rnd.random() < 0.55:
            args += [q, rnd.choice(LIST_VALUES[:12] if rnd.random() < 0.85 else LIST_VALUES)]


def gen_msbar(rnd, i, d):
    kind = rnd.choice(['n', 'n', 'p'])
    short = rnd.random() < 0.12
    items = G.seqs_for(rnd, 'msbar', kind, rnd.choice([1, 1, 1, 2, 3]), (1, 5) if short else (1, 150))
    text, fmt = G.seq_file(rnd, d, items)
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt)
    if rnd.random() < 0.7: args += ['-count', str(rnd.choice([0, 1, 1, 2, 3, 5, 10, -1]))]
    msb_lists(rnd, args)
    if rnd.random() < 0.3: args += ['-inframe']
    if rnd.random() < 0.4: args += ['-minimum', str(rnd.choice([0, 1, 2, 3, 5, 10, 30, -2]))]
    if rnd.random() < 0.4: args += ['-maximum', str(rnd.choice([0, 1, 3, 5, 10, 20, 60]))]
    files = {'in.dat': text}
    r = rnd.random()
    if r < 0.12:
        args += ['-othersequence', 'asis:' + rnd.choice([items[0][2][:rnd.randint(1, 6)], 'N', items[0][2]])]
    elif r < 0.22:
        files['other.dat'] = L.fasta(rnd, [(n, '', s) for n, _, s in items] + [('z', '', 'ACGT')])
        args += ['-othersequence', 'other.dat']
    elif r < 0.25:
        args += ['-othersequence', 'missing.dat']
    if rnd.random() < 0.03: args += [rnd.choice(['-count', '-minimum', '-maximum']), rnd.choice(['x', '1.5'])]
    ex = clock(rnd)
    if short or any(len(s) <= 6 for _, _, s in items) or '3' in ''.join(a for a in args if a in ('3', '4', '3,4,5', '2,3,4,5', 'Del', 'ins,del', '*', 'Insertions')) or True:
        ex['vg'] = True
    return args + G.out_quals(rnd) + ['-outseq', 'stdout'], files, ex


def pepstats_file(rnd):
    aas = 'ARNDCQEGHILKMFPSTWYV'
    names = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}
    cnt = [rnd.randint(0, 30) for _ in aas]
    tot = max(1, sum(cnt))
    out = 'PEPSTATS of SEQ from 1 to %d\n\nMolecular weight = 1234.5\t\tResidues = %d\nAverage Weight  = 110\t\tCharge   = 1\n\n\nResidue\t\tNumber\t\tMole%%\t\tDayhoffStat\n' % (tot, tot)
    for a, c in zip(aas, cnt):
        out += '%s = %s\t\t%d\t\t%.3f\t\t%.3f\n' % (a, names[a], c, 100.0 * c / tot, rnd.random() * 2)
    out += '\nProperty\tResidues\t\tNumber\t\tMole%\n'
    if rnd.random() < 0.15:
        out = out.replace('PEPSTATS', rnd.choice(['pepstats', 'Something', '']))
    if rnd.random() < 0.1:
        out = out.replace('Residue\t', 'Residues\t', 1) if rnd.random() < 0.5 else out.split('Residue')[0]
    if rnd.random() < 0.1:
        out = out.replace('\n\n', '\n')
    if rnd.random() < 0.1:
        out = out.replace('\t', ' ')
    if rnd.random() < 0.05:
        out += '\n'.join(['Z\t1'])
    return out


def cusp_file(rnd):
    codons = [a + b + c for a in 'ACGT' for b in 'ACGT' for c in 'ACGT']
    code = 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF'
    cnt = {c: rnd.choice([0, 1, 5, 20, 100, 400]) for c in codons}
    tot = max(1, sum(cnt.values()))
    out = '#CdsCount: 100\n#Species: Test species\n#Codon AA Fraction Frequency Number\n'
    for c, aa in zip(codons, code):
        n = cnt[c]
        out += '%s %s %.3f %.3f %d\n' % (c, aa, rnd.random(), 1000.0 * n / tot, n)
    if rnd.random() < 0.1: out = out.replace('.', ',', 3)
    if rnd.random() < 0.05: out += 'AAA K x y z\n'
    return out


def make_common(rnd, nuc):
    args = []
    if rnd.random() < 0.85: args += ['-amount', str(rnd.choice([1, 1, 2, 3, 5, 10, 50, 0, -1]))]
    if rnd.random() < 0.85: args += ['-length', str(rnd.choice([1, 2, 5, 10, 30, 60, 100, 250, 1000, 0]))]
    files = {}
    if rnd.random() < 0.3:
        ins = ''.join(rnd.choice('ACGT' if nuc else 'ACDEFGHIKLMNPQRSTVWY') for _ in range(rnd.randint(1, 8)))
        if rnd.random() < 0.05: ins = rnd.choice(['xyz*', 'a b', 'NNNN', 'uuu', 'U', '123'])
        args += ['-insert', ins]
        if rnd.random() < 0.6: args += ['-start', str(rnd.choice([1, 2, 3, 5, 10, 100, 0, -3]))]
        if rnd.random() < 0.3: args += ['-useinsert']
    if rnd.random() < 0.03: args += [rnd.choice(['-amount', '-length', '-start']), rnd.choice(['x', '1.5'])]
    return args, files


def gen_makenucseq(rnd, i, d):
    args, files = make_common(rnd, True)
    if rnd.random() < 0.25:
        files['usage.cut'] = cusp_file(rnd)
        args += ['-codonfile', 'usage.cut' if rnd.random() < 0.9 else 'nothere.cut']
    return args + G.out_quals(rnd) + ['-outseq', 'stdout'], files, clock(rnd)


def gen_makeprotseq(rnd, i, d):
    args, files = make_common(rnd, False)
    if rnd.random() < 0.3:
        files['stats.txt'] = pepstats_file(rnd) if rnd.random() < 0.9 else 'just text\n'
        args += ['-pepstatsfile', 'stats.txt' if rnd.random() < 0.95 else 'nothere.txt']
    return args + G.out_quals(rnd) + ['-outseq', 'stdout'], files, clock(rnd)


GENS = {'shuffleseq': gen_shuffleseq, 'msbar': gen_msbar, 'makenucseq': gen_makenucseq, 'makeprotseq': gen_makeprotseq}

if __name__ == '__main__':
    prog = sys.argv.pop(1)
    L.main(prog, GENS[prog], 400)
