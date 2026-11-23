#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/mafft_parttree_oracle.tsv.gz (B09 audit A2-19c, stage 4: `mafft --parttree`, `--dpparttree`, with --partsize / --groupsize / --retree / --reorder / --treeout) with the installed MAFFT 7.505 (`mafft`, Debian 7.505-1).

Every row = one run `mafft [--quiet] <options> in.fa` on a deterministic set (python `random.Random`, fixed seed); each run is repeated and a set whose two runs differ is reported (the program is deterministic: none did). Sets: 3 - 60 sequences of 6 - 150
letters (protein and DNA, related by substitutions and indels, with duplicates, unusual symbols with --anysymbol, '-' in the input), and for every group some of 100 - 450 sequences (the representatives are chosen at random: the program's own `rand()` with seed 1).
Columns: id, strategy, options, seqs ('|' separated), rc, rows ('|' separated, in the order of the output), names (output names), stderr_sha256, stderr (escaped; every 25th row), tree (the file input.tree, escaped).
Run: python3 scripts/alignment/generate_mafft_parttree_oracle.py   (env PT_SCALE, PT_OUT for a smaller run)
"""
import gzip, hashlib, os, random, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'mafft_parttree_oracle.tsv.gz')
rng = random.Random(75061)
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


def family(alpha, n, k, dup=0.0):
    base = ''.join(rng.choice(alpha) for _ in range(n))
    r = rng.choice([0.05, 0.15, 0.3])
    ind = rng.choice([0.02, 0.05, 0.1])
    out = [mut(base, alpha, r, ind) for _ in range(k)]
    for i in range(len(out)):
        if rng.random() < dup:
            out[i] = out[rng.randrange(len(out))]
    return out


def clans(alpha, n, k):
    # several related families (clans): the classification has something to find
    out = []
    for _ in range(rng.randint(2, 5)):
        out += family(alpha, n, max(2, k // 3))
    rng.shuffle(out)
    return out


def alpha():
    return P if rng.random() < 0.5 else N


def small_set():
    return family(alpha(), rng.choice([10, 20, 40, 60]), rng.randint(3, 30), dup=rng.choice([0, 0, 0.2]))


def mid_set():
    return clans(alpha(), rng.choice([12, 24]), rng.randint(30, 60))


def big_set():
    return clans(alpha(), rng.choice([12, 18]), rng.randint(100, 200))


def tiny_set():
    a = alpha()
    return [''.join(rng.choice(a) for _ in range(rng.randint(6, 14))) for _ in range(rng.randint(2, 9))]


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


def unrelated_set():
    a = alpha()
    return [''.join(rng.choice(a) for _ in range(rng.randint(6, 80))) for _ in range(rng.randint(3, 30))]


def pair_set():
    return family(alpha(), rng.choice([20, 60, 150]), 2)


SCALE = float(os.environ.get('PT_SCALE', '1'))
OUT = os.environ.get('PT_OUT')
if OUT:
    out_path = OUT

PT, DP = ['--parttree'], ['--dpparttree']


def groups():
    """(name, base options, variants (one is chosen at random per set), [(generator, count)]): a group has at least 200 sets."""
    g = []
    for k, base in (('', PT), ('dp-', DP)):
        dp = base is DP
        big = 4 if dp else 6
        g.append((k + 'default', base, [[]], [(small_set, 130), (mid_set, 25), (tiny_set, 20), (unrelated_set, 15), (pair_set, 10), (big_set, big)]))
        g.append((k + 'partsize', base, [['--partsize', '3'], ['--partsize', '5'], ['--partsize', '8'], ['--partsize', '12'], ['--partsize', '20']], [(small_set, 120), (mid_set, 30), (tiny_set, 20), (unrelated_set, 15), (big_set, big)]))
        g.append((k + 'groupsize', base, [['--groupsize', '2'], ['--groupsize', '5'], ['--groupsize', '8', '--partsize', '10'], ['--groupsize', '20'], ['--groupsize', '40', '--partsize', '10'], ['--groupsize', '1'], ['--groupsize', '1', '--treeout']],
                  [(small_set, 120), (mid_set, 40), (tiny_set, 20), (big_set, big)]))
        g.append((k + 'reorder-treeout', base, [['--reorder'], ['--treeout'], ['--reorder', '--treeout'], ['--treeout', '--partsize', '8'], ['--reorder', '--partsize', '5']], [(small_set, 120), (mid_set, 40), (tiny_set, 20), (pair_set, 10), (big_set, big)]))
        g.append((k + 'retree', base, [['--retree', '0'], ['--retree', '0', '--treeout'], ['--retree', '1'], ['--retree', '3'], ['--retree', '3', '--reorder', '--partsize', '8'], ['--retree', '3', '--treeout']],
                  [(small_set, 130), (mid_set, 40), (tiny_set, 15), (big_set, big)]))
        g.append((k + 'anysymbol', base, [['--anysymbol'], ['--anysymbol', '--reorder'], ['--anysymbol', '--partsize', '5']], [(odd_set, 190), (small_set, 20)]))
        g.append((k + 'model', base, [['--nuc'], ['--amino'], ['--bl', '80', '--op', '2.5', '--ep', '0.3'], ['--jtt', '120', '--averagelinkage'], ['--leavegappyregion', '--nofft'], ['--minimumlinkage', '--bl', '45']],
                  [(small_set, 160), (mid_set, 40), (big_set, big)]))
    # every group has at least 205 sets: the first generator of a group is topped up
    out = []
    for name, base, variants, gens in g:
        total = sum(c for _, c in gens)
        if total < 205:
            gens = [(gens[0][0], gens[0][1] + 205 - total)] + list(gens[1:])
        out.append((name, base, variants, gens))
    return out


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\r', '\\r').replace('\b', '\\b').replace('\t', '\\t')


def parse(txt):
    rows, names = [], []
    for l in txt.splitlines():
        if l.startswith('>'):
            names.append(l[1:])
            rows.append('')
        elif l.strip():
            rows[-1] += l.strip()
    return rows, names


def once(opts, seqs):
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, 'in.fa')
        open(fa, 'w').write(''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs)))
        p = subprocess.run(['mafft'] + opts + [fa], capture_output=True, cwd=td)
        tf = os.path.join(td, 'in.fa.tree')
        tree = open(tf, 'rb').read().decode('latin-1') if os.path.exists(tf) else ''
    return p, tree


def run(case):
    k, name, opts, seqs = case
    p, tree = once(opts, seqs)
    p2, tree2 = once(opts, seqs)
    same = p.stdout == p2.stdout and p.stderr == p2.stderr and p.returncode == p2.returncode and tree == tree2
    rows, names = parse(p.stdout.decode('latin-1'))
    digest = hashlib.sha256(p.stderr).hexdigest()
    full = esc(p.stderr.decode('latin-1')) if (k % 25 == 0 and len(seqs) <= 40) else ''
    return same, '\t'.join([str(k), name, ' '.join(opts), '|'.join(seqs), str(p.returncode), '|'.join(rows), '|'.join(names), digest, full, esc(tree)])


cases = []
k = 0
for name, base, variants, gens in groups():
    for gen, count in gens:
        for _ in range(max(1, int(count * SCALE))):
            k += 1
            quiet = (k % 10 == 9)
            cases.append((k, name, (['--quiet'] if quiet else []) + base + rng.choice(variants), gen()))

with ThreadPoolExecutor(max_workers=4) as ex:
    res = list(ex.map(run, cases))
unstable = [c[0] for c, r in zip(cases, res) if not r[0]]
print('sets', len(cases), 'runs repeated, differing between the two runs:', len(unstable), unstable[:10])
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with gzip.open(out_path, 'wt', encoding='utf-8') as f:
    f.write('\n'.join(r[1] for r in res) + '\n')
print('wrote', len(res), 'rows to', out_path)
