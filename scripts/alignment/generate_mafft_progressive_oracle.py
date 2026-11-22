#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/mafft_progressive_oracle.tsv.gz (B09 audit A2-19b, stage 1: three and more sequences, the progressive strategies) with the installed MAFFT 7.505 (`mafft`, Debian 7.505-1).

Every row = one run `mafft [--quiet] <options> in.fa` on a deterministic set (python `random.Random`, fixed seed) of 3 - 45 sequences of 1 - 900 residues (related by substitutions and indels, unrelated, with
IUPAC / unusual symbols, '-' in the input, lower case): the rows of the alignment (unwrapped, in the order of the output), the return code and the SHA-256 of the program's standard error of the same run (the progress and strategy text; the full text, `\\n` `\\r` `\\b` escaped,
is kept for every 25th row; with --thread the text depends on the order the threads run: no hash).  Strategies: the default FFT-NS-2, --retree 1 (FFT-NS-1), --nofft (NW-NS-2), --retree 1 --nofft (NW-NS-1), --retree 3, --6merpair, --10merpair, --op / --ep / --exp, --bl, --jtt, --tm, --kimura,
--kappa, --fmodel, --nwildcard, --leavegappyregion, --termgappenalty, --unweight, --averagelinkage, --minimumlinkage, --mixedlinkage, --opdist, --reorder, --thread 4, --anysymbol, --nuc / --amino.
Columns: id, strategy, options, seqs ('|' separated), rc, rows ('|' separated), stderr_sha256, stderr (escaped; every 25th row).  Run: python3 scripts/alignment/generate_mafft_progressive_oracle.py
"""
import gzip, hashlib, os, random, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'mafft_progressive_oracle.tsv.gz')
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
    return family(P, rng.choice([8, 20, 40, 80, 150, 250]), rng.randint(3, 12))


def dna_set():
    return family(N, rng.choice([12, 30, 60, 120, 250, 400, 600]), rng.randint(3, 12))


def big_set():
    return family(P if rng.random() < 0.5 else N, rng.choice([30, 60, 120]), rng.randint(21, 45))


def tiny_set():
    alpha = P if rng.random() < 0.5 else N
    return [''.join(rng.choice(alpha) for _ in range(rng.randint(1, 14))) for _ in range(rng.randint(3, 8))]


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
strategies = [
    ('FFT-NS-2', [], protein_set, 150), ('FFT-NS-2', [], dna_set, 150), ('FFT-NS-2', [], big_set, 60), ('FFT-NS-2', [], tiny_set, 60), ('FFT-NS-2', [], unrelated_set, 40),
    ('FFT-NS-2', ['--anysymbol'], odd_set, 100), ('FFT-NS-2', [], odd_set, 20),
    ('FFT-NS-1', ['--retree', '1'], protein_set, 110), ('FFT-NS-1', ['--retree', '1'], dna_set, 110),
    ('NW-NS-2', ['--nofft'], protein_set, 110), ('NW-NS-2', ['--nofft'], dna_set, 110),
    ('NW-NS-1', ['--retree', '1', '--nofft'], protein_set, 110), ('NW-NS-1', ['--retree', '1', '--nofft'], dna_set, 110),
    ('retree3', ['--retree', '3'], protein_set, 110), ('retree3', ['--retree', '3'], dna_set, 110),
    ('6merpair', ['--6merpair'], protein_set, 110), ('6merpair', ['--6merpair'], dna_set, 110),
    ('10merpair', ['--10merpair'], dna_set, 220),
    ('op-ep', ['--op', '2.5', '--ep', '0.3'], protein_set, 110), ('op-ep', ['--op', '0.8', '--ep', '0.1'], dna_set, 110),
    ('exp', ['--exp', '0.2'], protein_set, 110), ('exp', ['--exp', '0.5', '--op', '3'], dna_set, 110),
    ('opdist', ['--opdist', '3.0'], protein_set, 110), ('opdist', ['--opdist', '0.5'], dna_set, 110),
    ('bl', ['--bl', '80'], protein_set, 70), ('bl', ['--bl', '45'], protein_set, 70), ('bl', ['--bl', '30'], protein_set, 70),
    ('jtt', ['--jtt', '100'], protein_set, 110), ('jtt', ['--jtt', '250'], dna_set, 110),
    ('tm', ['--tm', '150'], protein_set, 220),
    ('kimura', ['--kimura', '100'], dna_set, 220), ('kappa', ['--kappa', '4'], dna_set, 220),
    ('fmodel', ['--fmodel'], protein_set, 110), ('fmodel', ['--fmodel'], dna_set, 110),
    ('nwildcard', ['--nwildcard'], dna_set, 220),
    ('leavegappyregion', ['--leavegappyregion'], protein_set, 110), ('leavegappyregion', ['--leavegappyregion'], dna_set, 110),
    ('termgappenalty', ['--termgappenalty'], protein_set, 110), ('termgappenalty', ['--termgappenalty'], dna_set, 110),
    ('unweight', ['--unweight'], protein_set, 110), ('unweight', ['--unweight'], dna_set, 110),
    ('averagelinkage', ['--averagelinkage'], protein_set, 110), ('averagelinkage', ['--averagelinkage'], dna_set, 110), ('minimumlinkage', ['--minimumlinkage'], protein_set, 110), ('minimumlinkage', ['--minimumlinkage'], dna_set, 110),
    ('mixedlinkage', ['--mixedlinkage', '0.4'], protein_set, 110), ('mixedlinkage', ['--mixedlinkage', '0.7'], dna_set, 110),
    ('reorder', ['--reorder'], protein_set, 110), ('reorder', ['--reorder'], dna_set, 110),
    ('fft', ['--fft'], protein_set, 110), ('fft', ['--fft'], dna_set, 110),
    ('thread', ['--thread', '4'], protein_set, 110), ('thread', ['--thread', '4'], big_set, 110),
    ('nuc', ['--nuc'], dna_set, 220), ('amino', ['--amino'], protein_set, 220),
]


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
    for _ in range(20):
        k += 1
        cases.append((k, 'FFT-NS-2', [], odd_set()))


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
