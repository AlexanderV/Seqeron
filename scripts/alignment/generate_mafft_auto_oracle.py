#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/mafft_auto_oracle.tsv.gz (B09 audit A2-19c, stage 6: `mafft --auto`) with the installed MAFFT 7.505 (`mafft`, Debian 7.505-1).

Every row = one run `mafft [--quiet] --auto <options> in.fa` on a deterministic set (python `random.Random`, fixed seed): fewer than 100 sequences of 1 - 120 residues (L-INS-i: related by substitutions and indels, unrelated, tiny, with unusual symbols, '-' in the input),
100 - 130 sequences (L-INS-2 iterations), 200 - 215 sequences (FFT-NS-i, two iterations) and 500 - 510 sequences (FFT-NS-2): the rows of the alignment (unwrapped), the return code and the SHA-256 of the program's standard error (the full text, `\\n` `\\r` `\\b`
escaped, for every 25th row; with --thread no hash).  Groups: the default and the options --op / --ep, --bl, --jtt, --unweight, --reorder, --leavegappyregion, --fmodel, --nuc / --amino, --thread 4 --threadit 0.
Columns: id, strategy, options, seqs ('|' separated), rc, rows ('|' separated), stderr_sha256, stderr (escaped; every 25th row).  Run: python3 scripts/alignment/generate_mafft_auto_oracle.py
"""
import gzip, hashlib, os, random, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'mafft_auto_oracle.tsv.gz')
rng = random.Random(75051)
P = 'ACDEFGHIKLMNPQRSTVWY'
N = 'ACGT'


def mut(s, alpha, r, ind):
    out = []
    for c in s:
        x = rng.random()
        if x < ind / 2:
            continue
        if x < ind:
            out.append(c)
            out.append(rng.choice(alpha))
            continue
        out.append(rng.choice(alpha) if rng.random() < r else c)
    return ''.join(out) or s[:1] or alpha[0]


def family(alpha, n, k):
    base = ''.join(rng.choice(alpha) for _ in range(n))
    r = rng.choice([0.05, 0.15, 0.3])
    ind = rng.choice([0.02, 0.05, 0.1])
    return [mut(base, alpha, r, ind) for _ in range(k)]


def protein_set():
    # the refinement of a pair library re-aligns every branch of the tree many times: the sets are smaller than in the oracle of the progressive pair-library methods
    return family(P, rng.choice([8, 20, 40, 80]), rng.randint(3, 8))


def dna_set():
    return family(N, rng.choice([12, 30, 60, 120]), rng.randint(3, 8))


def mid_set():
    # 100 - 199 sequences of fewer than 1 000 letters: the second branch of --auto
    return family(P if rng.random() < 0.5 else N, rng.choice([20, 40, 60]), rng.randint(100, 130))


def large_set():
    # 200 - 499 sequences: the third branch (FFT-NS-i, two iterations, two trees)
    return family(P if rng.random() < 0.5 else N, rng.choice([10, 20]), rng.randint(200, 215))


def xl_set():
    # 500 and more sequences: the fourth branch (FFT-NS-2)
    return family(P if rng.random() < 0.5 else N, rng.choice([8, 12]), rng.randint(500, 510))


def tiny_set():
    alpha = P if rng.random() < 0.5 else N
    return [''.join(rng.choice(alpha) for _ in range(rng.randint(1, 14))) for _ in range(rng.randint(3, 8))]


def pair_set():
    return family(P if rng.random() < 0.5 else N, rng.choice([20, 60, 150, 250]), 2)


def odd_set():
    alpha = P if rng.random() < 0.5 else N
    seqs = family(alpha, rng.choice([20, 40, 80]), rng.randint(3, 9))
    extra = 'XBZJOU*#.@x' if alpha == P else 'NRYKMSWBDHVXnryx'
    out = []
    for s in seqs:
        s = list(s)
        for _ in range(rng.randrange(0, 4)):
            s[rng.randrange(len(s))] = rng.choice(extra)
        if rng.random() < 0.4:
            s = [c.lower() if rng.random() < 0.5 else c for c in s]
        for _ in range(rng.randrange(0, 3)):
            s.insert(rng.randrange(len(s) + 1), '-')
        out.append(''.join(s))
    return out


def unrelated_set():
    alpha = P if rng.random() < 0.5 else N
    return [''.join(rng.choice(alpha) for _ in range(rng.randint(5, 150))) for _ in range(rng.randint(3, 9))]


# (strategy name, mafft options, generator, number of sets)

SCALE = float(os.environ.get('AUTO_SCALE', '1'))
OUT = os.environ.get('AUTO_OUT')
A = ['--auto']


def build():
    out = [
        ('auto', A, protein_set, 80), ('auto', A, dna_set, 80), ('auto', A, tiny_set, 30), ('auto', ['--anysymbol'] + A, odd_set, 30), ('auto', A, unrelated_set, 20), ('auto', A, pair_set, 20),
        ('auto', A, mid_set, 12), ('auto', A, large_set, 10), ('auto', A, xl_set, 6),
    ]
    for name, opts in [('op-ep', ['--op', '2.5', '--ep', '0.3']), ('bl', ['--bl', '80']), ('jtt', ['--jtt', '150']), ('unweight', ['--unweight']), ('reorder', ['--reorder']),
                       ('leavegappyregion', ['--leavegappyregion']), ('fmodel', ['--fmodel']), ('threadit0', ['--thread', '4', '--threadit', '0'])]:
        out += [('auto-' + name, A + opts, protein_set, 110), ('auto-' + name, A + opts, dna_set, 110)]
    out += [('auto-nuc', A + ['--nuc'], dna_set, 220), ('auto-amino', A + ['--amino'], protein_set, 220)]
    return out


strategies = [(n, o, g, max(1, int(c * SCALE))) for n, o, g, c in build()]
if OUT:
    out_path = OUT


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\r', '\\r').replace('\b', '\\b').replace('\t', '\\t')


def unwrap(txt):
    rows = []
    for l in txt.splitlines():
        if l.startswith('>'):
            rows.append('')
        elif l.strip():
            rows[-1] += l.strip()
    return rows


cases = []
k = 0
for name, opts, gen, count in strategies:
    for _ in range(count):
        k += 1
        quiet = (k % 10 == 9)
        cases.append((k, name, (['--quiet'] if quiet else []) + opts, gen()))
if True:
    # one extra block: the plain command with letters outside the alphabet (the program stops with a message, rc 1)
    for _ in range(max(1, int(20 * SCALE))):
        k += 1
        cases.append((k, 'auto', ['--auto'], odd_set()))


def run(case):
    k, name, opts, seqs = case
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa')
        open(fa, 'w').write(''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs)))
        p = subprocess.run(['mafft'] + opts + [fa], capture_output=True)  # bytes: the progress text holds carriage returns and backspaces
    threaded = '--thread' in opts
    digest = '-' if threaded else hashlib.sha256(p.stderr).hexdigest()
    full = esc(p.stderr.decode('latin-1')) if (k % 25 == 0 and not threaded) else ''
    return '\t'.join([str(k), name, ' '.join(opts), '|'.join(seqs), str(p.returncode), '|'.join(unwrap(p.stdout.decode('latin-1'))), digest, full])


with ThreadPoolExecutor(max_workers=4) as ex:
    rows_out = list(ex.map(run, cases))
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with gzip.open(out_path, 'wt', encoding='utf-8') as f:
    f.write('\n'.join(rows_out) + '\n')
print('wrote', len(rows_out), 'rows to', out_path)
