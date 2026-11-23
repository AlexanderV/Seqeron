#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/mafft_treeout_oracle.tsv.gz (B09 audit A2-19c, stage 4: `mafft --retree 0`, `--treeout`, `--nodeout`, `--distout`) with the installed MAFFT 7.505 (`mafft`, Debian 7.505-1).

Every row = one run `mafft [--quiet] <options> in.fa` on a deterministic set (python `random.Random`, fixed seed; sequences named s0, s1, ...): the rows (the sequences of --retree 0 are not aligned), the names of the output, the return code, the SHA-256 of the
standard error (the full text for every 25th row), the guide tree file `in.fa.tree` and the distance file `in.fa.hat2`. Every run is repeated; a set whose two runs differ is reported (none did).
Columns: id, strategy, options, seqs ('|' separated), rc, rows, names, stderr_sha256, stderr (escaped; every 25th row), tree (escaped), hat2 (escaped).  Run: python3 scripts/alignment/generate_mafft_treeout_oracle.py   (env TO_SCALE, TO_OUT)
"""
import gzip, hashlib, os, random, subprocess, tempfile
from concurrent.futures import ThreadPoolExecutor

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'mafft_treeout_oracle.tsv.gz')
rng = random.Random(75062)
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


def alpha():
    return P if rng.random() < 0.5 else N


def small_set():
    return family(alpha(), rng.choice([10, 30, 60, 120]), rng.randint(3, 25))


def pair_set():
    return family(alpha(), rng.choice([20, 60, 150]), 2)


def big_set():
    return family(alpha(), rng.choice([15, 40]), rng.randint(60, 120))


def tiny_set():
    a = alpha()
    return [''.join(rng.choice(a) for _ in range(rng.randint(1, 14))) for _ in range(rng.randint(3, 9))]


def odd_set():
    a = alpha()
    seqs = family(a, rng.choice([20, 40, 80]), rng.randint(3, 12))
    extra = 'XBZJOU*#.@x' if a == P else 'NRYKMSWBDHVXnryx'
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


def gapped_set():
    seqs = small_set()
    return [s[:3] + '-' * rng.randint(0, 3) + s[3:] for s in seqs]


def groups():
    """(name, base options, variants (one is chosen at random per set), [(generator, count)]): a group has at least 200 sets."""
    g = []
    sets = [(small_set, 130), (big_set, 8), (tiny_set, 20), (pair_set, 15), (gapped_set, 30)]
    g.append(('retree0', [], [['--retree', '0'], ['--retree', '0', '--treeout'], ['--retree', '0', '--nodeout'], ['--retree', '0', '--distout'], ['--retree', '0', '--treeout', '--distout'], ['--6merpair', '--retree', '0', '--treeout']], sets))
    g.append(('treeout', [], [['--treeout'], ['--nodeout'], ['--distout'], ['--treeout', '--distout'], ['--retree', '1', '--treeout'], ['--retree', '3', '--treeout', '--distout']], sets))
    g.append(('model', [], [['--nofft', '--nodeout'], ['--averagelinkage', '--treeout'], ['--minimumlinkage', '--distout'], ['--10merpair', '--nuc', '--treeout'], ['--bl', '80', '--op', '2.5', '--treeout'], ['--leavegappyregion', '--treeout']], sets))
    g.append(('anysymbol', ['--anysymbol'], [['--treeout', '--distout'], ['--retree', '0', '--distout'], ['--retree', '0', '--treeout'], ['--nodeout']], [(odd_set, 190), (small_set, 20)]))
    g.append(('plain-odd', [], [['--retree', '0']], [(odd_set, 210)]))
    return g


SCALE = float(os.environ.get('TO_SCALE', '1'))
if os.environ.get('TO_OUT'):
    out_path = os.environ['TO_OUT']


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\r', '\\r').replace('\b', '\\b').replace('\t', '\\t')


def once(opts, seqs):
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa')
        open(fa, 'w').write(''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs)))
        p = subprocess.run(['mafft'] + opts + [fa], capture_output=True, cwd=td)
        rd = lambda n: open(os.path.join(td, n), 'rb').read().decode('latin-1') if os.path.exists(os.path.join(td, n)) else ''
        return p, rd('in.fa.tree'), rd('in.fa.hat2')


def run(case):
    k, name, opts, seqs = case
    p, tree, hat = once(opts, seqs)
    p2, tree2, hat2 = once(opts, seqs)
    same = p.stdout == p2.stdout and p.stderr == p2.stderr and p.returncode == p2.returncode and tree == tree2 and hat == hat2
    rows, names = [], []
    for l in p.stdout.decode('latin-1').splitlines():
        if l.startswith('>'):
            names.append(l[1:])
            rows.append('')
        elif l.strip():
            rows[-1] += l.strip()
    full = esc(p.stderr.decode('latin-1')) if k % 25 == 0 else ''
    return same, '\t'.join([str(k), name, ' '.join(opts), '|'.join(seqs), str(p.returncode), '|'.join(rows), '|'.join(names), hashlib.sha256(p.stderr).hexdigest(), full, esc(tree), esc(hat)])


cases = []
k = 0
for name, base, variants, gens in groups():
    for gen, count in gens:
        for _ in range(max(1, int(count * SCALE))):
            k += 1
            cases.append((k, name, (['--quiet'] if k % 10 == 9 else []) + base + rng.choice(variants), gen()))
with ThreadPoolExecutor(max_workers=3) as ex:
    res = list(ex.map(run, cases))
unstable = [c[0] for c, r in zip(cases, res) if not r[0]]
print('sets', len(cases), 'runs repeated, differing between the two runs:', len(unstable), unstable[:10])
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with gzip.open(out_path, 'wt', encoding='utf-8') as f:
    f.write('\n'.join(r[1] for r in res) + '\n')
print('wrote', len(res), 'rows to', out_path)
