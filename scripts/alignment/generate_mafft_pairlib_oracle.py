#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/mafft_pairlib_oracle.tsv.gz (B09 audit A2-19b / A2-19c, stage 3: the progressive methods with a pair library, --globalpair (G-INS-1), --localpair (L-INS-1), --genafpair (E-INS-1)) with the installed MAFFT 7.505 (`mafft`, Debian 7.505-1).

Every row = one run `mafft [--quiet] <options> in.fa` on a deterministic set (python `random.Random`, fixed seed) of 2 - 28 sequences of 1 - 250 residues (related by substitutions and indels, unrelated, with IUPAC / unusual symbols, '-' in the
input, lower case): the rows of the alignment (unwrapped, in the order of the output), the return code and the SHA-256 of the program's standard error of the same run (the progress and strategy text; the full text, `\\n` `\\r` `\\b` escaped, is kept
for every 25th row; with --thread the text depends on the order the threads run: no hash).  Strategies per pair library: the default, --op / --ep, the library's own gap options (--gop --gep --gexp / --lop --lep --lexp / --LOP --LEXP), --weighti, --bl, --jtt, --tm,
--kappa, --fmodel, --nwildcard, --leavegappyregion, --termgappenalty, --unweight, --reorder, --thread 4, --anysymbol, --nuc / --amino, two sequences.
Columns: id, strategy, options, seqs ('|' separated), rc, rows ('|' separated), stderr_sha256, stderr (escaped; every 25th row).  Run: python3 scripts/alignment/generate_mafft_pairlib_oracle.py
"""
import gzip, hashlib, os, random, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'mafft_pairlib_oracle.tsv.gz')
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
    # a pair library aligns all pairs: the sizes are kept smaller than in the oracles of the progressive methods (the committed test data must run in minutes)
    return family(P, rng.choice([8, 20, 40, 80, 150]), rng.randint(3, 10))


def dna_set():
    return family(N, rng.choice([12, 30, 60, 120, 250]), rng.randint(3, 10))


def big_set():
    # a pair library aligns all n (n - 1) / 2 pairs: 21 - 28 short sequences are the largest groups kept in the committed test data
    return family(P if rng.random() < 0.5 else N, rng.choice([30, 60]), rng.randint(21, 28))


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

MODES = [('G', '--globalpair', ['--gop', '-3.0', '--gep', '0.2', '--gexp', '-0.2']),
         ('L', '--localpair', ['--lop', '-3.0', '--lep', '0.2', '--lexp', '-0.2']),
         ('E', '--genafpair', ['--lop', '-2.5', '--LOP', '-4.0', '--LEXP', '-0.1'])]
SCALE = float(os.environ.get('PAIRLIB_SCALE', '1'))
ONLY = os.environ.get('PAIRLIB_ONLY', '')
OUT = os.environ.get('PAIRLIB_OUT')


def build():
    out = []
    for tag, flag, own in MODES:
        if ONLY and tag not in ONLY:
            continue
        f = [flag]
        out += [
            (tag + '-INS-1', f, protein_set, 150), (tag + '-INS-1', f, dna_set, 150), (tag + '-INS-1', f, big_set, 40), (tag + '-INS-1', f, tiny_set, 50),
            (tag + '-INS-1', f, unrelated_set, 40), (tag + '-INS-1', ['--anysymbol'] + f, odd_set, 60), (tag + '-INS-1', f, odd_set, 12),
            (tag + '-own', f + own, protein_set, 110), (tag + '-own', f + own, dna_set, 110),
            (tag + '-op-ep', f + ['--op', '2.5', '--ep', '0.3'], protein_set, 110), (tag + '-op-ep', f + ['--op', '0.8', '--ep', '0.1'], dna_set, 110),
            (tag + '-exp', f + ['--exp', '0.2'], protein_set, 110), (tag + '-exp', f + ['--exp', '0.5', '--op', '3'], dna_set, 110),
            (tag + '-weighti', f + ['--weighti', '1.5'], protein_set, 70), (tag + '-weighti', f + ['--weighti', '4'], protein_set, 70), (tag + '-weighti', f + ['--weighti', '0.5'], dna_set, 70),
            (tag + '-bl', f + ['--bl', '80'], protein_set, 70), (tag + '-bl', f + ['--bl', '45'], protein_set, 70), (tag + '-bl', f + ['--bl', '30'], protein_set, 70),
            (tag + '-jtt', f + ['--jtt', '100'], protein_set, 110), (tag + '-jtt', f + ['--jtt', '250'], dna_set, 110), (tag + '-tm', f + ['--tm', '150'], protein_set, 220),
            (tag + '-kappa', f + ['--kappa', '4'], dna_set, 220), (tag + '-fmodel', f + ['--fmodel'], protein_set, 110), (tag + '-fmodel', f + ['--fmodel'], dna_set, 110),
            (tag + '-nwildcard', f + ['--nwildcard'], dna_set, 220),
            (tag + '-leavegappyregion', f + ['--leavegappyregion'], protein_set, 110), (tag + '-leavegappyregion', f + ['--leavegappyregion'], dna_set, 110),
            (tag + '-termgappenalty', f + ['--termgappenalty'], protein_set, 110), (tag + '-termgappenalty', f + ['--termgappenalty'], dna_set, 110),
            (tag + '-unweight', f + ['--unweight'], protein_set, 110), (tag + '-unweight', f + ['--unweight'], dna_set, 110),
            (tag + '-reorder', f + ['--reorder'], protein_set, 110), (tag + '-reorder', f + ['--reorder'], dna_set, 110),
            (tag + '-thread', f + ['--thread', '4'], protein_set, 110), (tag + '-thread', f + ['--thread', '4'], big_set, 110),
            (tag + '-two-sequences', f, pair_set, 220),
            (tag + '-nuc', f + ['--nuc'], dna_set, 220), (tag + '-amino', f + ['--amino'], protein_set, 220),
        ]
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
    for tag, flag, _own in MODES:
        if ONLY and tag not in ONLY:
            continue
        for _ in range(max(1, int(20 * SCALE))):
            k += 1
            cases.append((k, tag + '-INS-1', [flag], odd_set()))


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
