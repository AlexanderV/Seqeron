#!/usr/bin/env python3
"""GENOMIC-TANDEM-001 oracle (A3-8): EMBOSS 6.6.0 `etandem` and `equicktandem` run on random and planted nucleotide sequences.

The C# ports (GenomicAnalyzer.FindTandemRepeatsEtandem / FindTandemRepeatsEquicktandem) are statement-by-statement ports of etandem.c / equicktandem.c
(raw.githubusercontent.com/kimrutherford/EMBOSS); this file records what the installed binaries print so the tests can lock it.

Sequences: random flanks over ACGT, planted tandem arrays (unit 1-14, 2-14 copies, point mutations, one-base indels), optional N runs, IUPAC symbols (R Y K M S W) and
lower case. Parameters: etandem -minrepeat 2..6 -maxrepeat min..min+9 -threshold 2..30 (a few <= 0) -mismatch -uniform; equicktandem -maxrepeat 1..14 -threshold 1..25.
Each (program, parameters, sequence) is run with `-origfile` (the program's original output); the rows are
  etandem      : score start end size count identity consensus
  equicktandem : score start end size count
TSV (tests/.../TestData/Analysis/emboss_tandem_oracle.tsv):  program \t params \t sequence \t rows ('|' separated, fields ',' separated; '-' = none)
params: etandem  min,max,threshold,mismatch(0/1),uniform(0/1)   equicktandem  max,threshold
Never aborts: a failing program call is reported on stderr and the row skipped.
Usage: python3 scripts/analysis/generate_emboss_tandem_oracle.py [n_etandem n_equick] > emboss_tandem_oracle.tsv
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261120)


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate:
            out.append(rnd.choice('ACGT'))
        elif r < rate * 1.2:
            continue
        elif r < rate * 1.4:
            out.append(c)
            out.append(rnd.choice('ACGT'))
        else:
            out.append(c)
    return ''.join(out)


def gen_sequence():
    parts = [''.join(rnd.choice('ACGT') for _ in range(rnd.randint(0, 40)))]
    for _ in range(rnd.randint(0, 3)):
        unit = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(1, 14)))
        arr = unit * rnd.randint(2, 14)
        parts.append(mutate(arr, rnd.choice([0.0, 0.0, 0.03, 0.08, 0.15])))
        parts.append(''.join(rnd.choice('ACGT') for _ in range(rnd.randint(0, 30))))
    s = ''.join(parts)
    if not s:
        s = 'ACGT'
    r = rnd.random()
    if r < 0.25:
        s = ''.join('N' if rnd.random() < 0.04 else c for c in s)
    elif r < 0.33:
        s = ''.join(rnd.choice('RYKMSW') if rnd.random() < 0.03 else c for c in s)
    r = rnd.random()
    if r < 0.15:
        s = s.lower()
    elif r < 0.25:
        s = ''.join(c.lower() if rnd.random() < 0.5 else c for c in s)
    return s


def parse(path, ncol):
    rows = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) >= ncol:
                rows.append(','.join(t[:ncol]))
    return '|'.join(rows) if rows else '-'


def run_one(tmp, prog, args, seq, ncol):
    fa = os.path.join(tmp, 's.fa')
    with open(fa, 'w') as f:
        f.write('>s\n' + seq + '\n')
    rep = os.path.join(tmp, 'o.rep')
    orig = os.path.join(tmp, 'o.orig')
    for p in (rep, orig):
        if os.path.exists(p):
            os.remove(p)
    r = subprocess.run([prog, '-sequence', fa, '-outfile', rep, '-origfile', orig, '-auto'] + args,
                       capture_output=True, text=True)
    if r.returncode != 0 or not os.path.exists(orig):
        sys.stderr.write('FAIL %s %s %s: %s\n' % (prog, args, seq, r.stderr.strip()[:200]))
        return None
    return parse(orig, ncol)


def main():
    n_e = int(sys.argv[1]) if len(sys.argv) > 1 else 1700
    n_q = int(sys.argv[2]) if len(sys.argv) > 2 else 1300
    with tempfile.TemporaryDirectory() as tmp:
        for i in range(n_e):
            seq = gen_sequence()
            mn = rnd.randint(2, 6)
            mx = mn + rnd.randint(0, 9)
            th = rnd.choice([rnd.randint(2, 30)] * 9 + [rnd.randint(-3, 1)])
            mm = int(rnd.random() < 0.3)
            un = int(rnd.random() < 0.25)
            args = ['-minrepeat', str(mn), '-maxrepeat', str(mx), '-threshold', str(th)]
            if mm: args += ['-mismatch']
            if un: args += ['-uniform']
            rows = run_one(tmp, 'etandem', args, seq, 7)
            if rows is not None:
                print('etandem\t%d,%d,%d,%d,%d\t%s\t%s' % (mn, mx, th, mm, un, seq, rows))
        for i in range(n_q):
            seq = gen_sequence()
            mx = rnd.randint(1, 14)
            th = rnd.randint(1, 25)
            rows = run_one(tmp, 'equicktandem', ['-maxrepeat', str(mx), '-threshold', str(th)], seq, 5)
            if rows is not None:
                print('equicktandem\t%d,%d\t%s\t%s' % (mx, th, seq, rows))


if __name__ == '__main__':
    main()
