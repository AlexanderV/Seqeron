#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/spoa_oracle.tsv (B09 audit A2-26c) with the installed SPOA 4.0.8 (`spoa`, Debian 4.0.8-1): the whole command line.

Every row = one run `spoa ARGS` in an empty directory with the input files of the row. Columns: id, args (JSON list), files (JSON {name: base64}), exit status (the shell's: 134 for SIGABRT),
standard output, standard error, written files (JSON {name: text}); text columns are escaped (backslash, newline, tab, carriage return).
The sequences are random families (substitutions / indels, low-complexity and repeat sets that force ties, long sets that need the int32 engine) in FASTA / FASTQ / gzip, plus hand-made inputs.
Run: python3 scripts/alignment/generate_spoa_oracle.py
"""
import base64
import gzip
import json
import os
import random
import subprocess
import tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'spoa_oracle.tsv')
BIN = os.environ.get('SPOA_BIN', 'spoa')
rng = random.Random(20261321)


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def mutate(s, p, alpha='ACGT'):
    out = []
    for c in s:
        r = rng.random()
        if r < p / 3:
            continue
        if r < 2 * p / 3:
            out.append(rng.choice(alpha))
            continue
        out.append(c)
        if rng.random() < p / 3:
            out.append(rng.choice(alpha))
    return ''.join(out) or s[:1]


def family(n, length, p, alpha='ACGT'):
    base = ''.join(rng.choice(alpha) for _ in range(length))
    return [mutate(base, p * rng.random() * 2, alpha) for _ in range(n)]


def fasta(seqs, wrap=0, crlf=False, names=None):
    nl = '\r\n' if crlf else '\n'
    out = []
    for i, s in enumerate(seqs):
        name = names[i] if names else 'read%d' % i
        out.append('>' + name + nl)
        if wrap:
            for k in range(0, len(s), wrap):
                out.append(s[k:k + wrap] + nl)
        else:
            out.append(s + nl)
    return ''.join(out)


def fastq(seqs, names=None, quals=None):
    out = []
    for i, s in enumerate(seqs):
        name = names[i] if names else 'r%d' % i
        q = quals[i] if quals else ''.join(chr(33 + rng.randrange(0, 60)) for _ in s)
        out.append('@%s desc\n%s\n+\n%s\n' % (name, s, q))
    return ''.join(out)


rows = []


def run(args, files):
    with tempfile.TemporaryDirectory() as d:
        for name, data in files.items():
            with open(os.path.join(d, name), 'wb') as fh:
                fh.write(data)
        r = subprocess.run([BIN] + args, cwd=d, capture_output=True, timeout=120)
        rc = r.returncode
        if rc < 0:
            rc = 128 - rc
        written = {}
        for name in os.listdir(d):
            if name not in files:
                with open(os.path.join(d, name), 'rb') as fh:
                    written[name] = fh.read().decode('latin-1')
        rows.append((json.dumps(args), json.dumps({k: base64.b64encode(v).decode() for k, v in files.items()}), str(rc),
                     r.stdout.decode('latin-1'), r.stderr.decode('latin-1'), json.dumps(written)))


def params():
    sub = rng.choice(['lin', 'aff', 'cvx', 'cvx', 'def'])
    a = []
    if rng.random() < 0.8:
        a += ['-m', str(rng.randrange(1, 9))]
    if rng.random() < 0.8:
        a += ['-n', str(-rng.randrange(1, 9))]
    if sub == 'lin':
        gg = -rng.randrange(1, 12)
        a += ['-g', str(gg), '-e', str(gg - rng.randrange(0, 1) if False else gg)]
    elif sub == 'aff':
        a += ['-g', str(-rng.randrange(6, 14)), '-e', str(-rng.randrange(1, 5)), '-q', str(-rng.randrange(6, 14)), '-c', str(-rng.randrange(1, 5))]
    elif sub == 'cvx':
        a += ['-g', str(-rng.randrange(4, 10)), '-e', str(-rng.randrange(3, 7)), '-q', str(-rng.randrange(8, 20)), '-c', str(-rng.randrange(1, 3))]
    return a


def add(args, text, name='in.fa'):
    run(args, {name: text.encode('latin-1') if isinstance(text, str) else text})


# ---- random families: all alignment types x gap models x result modes
for i in range(150):
    n = rng.choice([2, 3, 4, 5, 8, 12])
    length = rng.choice([12, 30, 60, 100, 150, 200])
    seqs = family(n, length, rng.choice([0.05, 0.15, 0.3, 0.5]))
    args = ['-l', str(rng.choice([0, 1, 2]))] + params()
    modes = rng.choice([[], ['-r', '1'], ['-r', '2'], ['-r', '3'], ['-r', '4'], ['-r', '0', '-r', '1'], ['-r', '1', '-r', '3']])
    args += modes
    if rng.random() < 0.2:
        args += ['-s']
    if rng.random() < 0.15:
        args += ['-d', 'g.dot']
    if rng.random() < 0.3:
        add(args + ['in.fq'], fastq(seqs), 'in.fq')
    else:
        add(args + ['in.fa'], fasta(seqs, wrap=rng.choice([0, 0, 60, 70, 10])))

# ---- ties: low complexity, repeats, identical sequences
for i in range(50):
    unit = ''.join(rng.choice('ACGT') for _ in range(rng.choice([1, 2, 3, 4])))
    seqs = [unit * rng.randrange(3, 40) for _ in range(rng.choice([2, 3, 5, 8]))]
    if rng.random() < 0.5:
        seqs = [mutate(s, 0.1) for s in seqs]
    args = ['-l', str(i % 3)] + params() + ['-r', rng.choice(['1', '2', '3'])]
    add(args + ['in.fa'], fasta(seqs))

# ---- long / large scores: the int32 engine and the int16 limits
for i in range(20):
    length = rng.choice([260, 400, 700])
    seqs = family(rng.choice([2, 3]), length, 0.1)
    args = ['-l', str(i % 3), '-m', str(rng.choice([20, 60, 127])), '-n', str(-rng.choice([20, 60, 127])), '-g', str(-rng.choice([20, 100])), '-e', str(-rng.choice([10, 50, 90])),
            '-q', str(-rng.choice([60, 120])), '-c', str(-rng.choice([5, 30])), '-r', '1']
    add(args + ['in.fa'], fasta(seqs))

# ---- protein / odd alphabets / case / IUPAC / strand
for i in range(30):
    alpha = rng.choice(['ACDEFGHIKLMNPQRSTVWY', 'ACGTN', 'acgt', 'ACGTRYKM', 'AC'])
    seqs = family(rng.choice([2, 3, 4]), rng.choice([20, 50]), 0.2, alpha)
    args = ['-l', str(i % 3), '-r', rng.choice(['1', '2', '4']), '-s'] if i % 2 else ['-l', str(i % 3), '-r', '2']
    add(args + ['in.fa'], fasta(seqs))
for i in range(20):
    base = ''.join(rng.choice('ACGT') for _ in range(80))
    comp = base[::-1].translate(str.maketrans('ACGT', 'TGCA'))
    seqs = [mutate(base, 0.1), mutate(comp, 0.1), mutate(base, 0.1), mutate(comp, 0.1)]
    add(['-l', str(i % 3), '-s', '-r', rng.choice(['1', '3', '4']), 'in.fa'], fasta(seqs))

# ---- FASTQ with quality weights ('!' = weight 0), gzip, CRLF, long names
for i in range(20):
    seqs = family(rng.choice([2, 3, 5]), rng.choice([30, 80]), 0.15)
    quals = [''.join(rng.choice('!"#I5+?') for _ in s) for s in seqs]
    add(['-l', str(i % 3), '-r', rng.choice(['1', '3', '4']), 'in.fq'], fastq(seqs, quals=quals), 'in.fq')
for i in range(8):
    seqs = family(3, 50, 0.2)
    run(['-r', '2', 'in.fa.gz'], {'in.fa.gz': gzip.compress(fasta(seqs).encode(), mtime=0)})
    run(['-r', '1', 'in.fq.gz'], {'in.fq.gz': gzip.compress(fastq(seqs).encode(), mtime=0)})
    run(['-r', '2', 'in.fasta'], {'in.fasta': fasta(seqs, crlf=True).encode()})
    run(['-r', '4', 'in.fna'], {'in.fna': fasta(seqs, wrap=7).encode()})
    run(['-r', '3', 'in.faa'], {'in.faa': fasta(seqs, names=['seq %d with description' % k for k in range(len(seqs))]).encode()})

# ---- command line variants
f = fasta(family(4, 40, 0.2)).encode()
variants = [
    [], ['-h'], ['--help'], ['--version'], ['--vers'], ['-h', 'nonexistent.fa'], ['-x'], ['--bogus'], ['--algorithm'], ['-m'], ['--result=3', 'in.fa'], ['--res', '3', 'in.fa'], ['--alg=2', '--res=1', 'in.fa'],
    ['--dot=a.dot', 'in.fa'], ['--dot', 'a.dot', '-r', '3', 'in.fa'], ['--strand', 'in.fa'], ['--strand-ambiguous=1', 'in.fa'], ['--help=1'], ['in.fa', '-r', '1'], ['-r1', 'in.fa'], ['-l1', '-m3', '-n-2', 'in.fa'],
    ['-sr1', 'in.fa'], ['-r', '1', '--', 'in.fa'], ['--', '-r', 'in.fa'], ['in.fa', 'in2.fa'], ['in.txt'], ['missing.fa'], ['-l', '3', 'in.fa'], ['-l', '300', 'in.fa'], ['-g', '1', 'in.fa'], ['-q', '2', 'in.fa'],
    ['-e', '5', 'in.fa'], ['-c', '5', 'in.fa'], ['-g', '-200', 'in.fa'], ['-m', '200', 'in.fa'], ['-m', 'abc', 'in.fa'], ['-r', '9', 'in.fa'], ['-r', '5', '-r', '0', 'in.fa'], ['-r', '1', '-r', '1', 'in.fa'],
    ['-r', '2', '-r', '4', '-d', 'x.dot', 'in.fa'], ['-d', '', 'in.fa'], ['-l', '1', '-g', '-4', '-e', '-4', 'in.fa'], ['-g', '-6', '-e', '-6', '-q', '-10', '-c', '-4', 'in.fa'], ['-g', '-8', '-e', '-8', 'in.fa'],
    ['-g', '-8', '-e', '-2', '-q', '-8', '-c', '-1', 'in.fa'], ['-g', '-8', '-e', '-2', '-q', '-12', '-c', '-2', 'in.fa'], ['-m', '0', '-n', '0', 'in.fa'], ['-m', '-5', 'in.fa'],
]
for v in variants:
    run(v, {'in.fa': f, 'in2.fa': f, 'in.txt': f})

# ---- hand-made inputs (valid and invalid)
hand = [
    ('in.fa', ''), ('in.fa', 'ACGT\n'), ('in.fa', '>a\n>b\nACGT\n'), ('in.fa', '>a\nACGT\n>b\n'), ('in.fa', '>a\nAC GT\n>b\nACGA\n'), ('in.fa', '>\nACGT\n'), ('in.fa', '>a\nACGT\n\n\n>b\nAC\n'),
    ('in.fa', '>a\nACGT\n>b\nACGA'), ('in.fa', '>a\nACGTA>b\nACGA\n'), ('in.fa', '>a\nAAAA\n'), ('in.fa', '>a b c\nACGT\n>d\tx\nACGG\n'), ('in.fa', '>a\n  ACGT  \n>b\nACGT\n'),
    ('in.fq', '@a\nACGT\n+\nIIII\n'), ('in.fq', '@a\nACGT\n+\nIII\n'), ('in.fq', '@a\nACGT\n+\nIIII\n@b\nACGA\n+b\n!!!!\n'), ('in.fq', 'a\nACGT\n+\nIIII\n'), ('in.fq', '@a\nACGT\n'), ('in.fq', ''),
    ('in.fq', '@a\nACGT\n+\nIIII'), ('in.fq', '@a\nACGT\n+\n    \n@b\nACGT\n+\nIIII\n'), ('in.fq', '@a\nACGT\n+\n\x01\x01\x01\x01\n@b\nACGT\n+\nIIII\n'),
]
for name, text in hand:
    for r in (['-r', '2'], ['-r', '3'], ['-l', '1', '-r', '1']):
        run(r + [name], {name: text.encode('latin-1')})

os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='\n') as fh:
    fh.write('id\targs\tfiles\trc\tstdout\tstderr\toutfiles\n')
    for i, r in enumerate(rows):
        fh.write('S%04d\t%s\n' % (i + 1, '\t'.join(esc(c) for c in r)))
print(len(rows), 'rows;', sum(1 for r in rows if r[2] != '0'), 'with a non-zero exit status')
