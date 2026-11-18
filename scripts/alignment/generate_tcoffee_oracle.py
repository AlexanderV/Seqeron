#!/usr/bin/env python3
"""Oracle generator for SequenceAligner.MultipleAlignTCoffee (audit A1-13, B09 F58).

Runs the installed T-Coffee program (Debian package t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `t_coffee`)
on small random DNA sets with an EXPLICIT primary library (TC_LIB_FORMAT_01) and an EXPLICIT binary guide tree
(`-usetree`), once with `-relax_lib 1` (default) and once with `-relax_lib 0`, and writes the program's rows.

Library kinds:
  rand    : pairs of an own global (Needleman-Wunsch) + local (Smith-Waterman) alignment of every sequence pair, weight 10 x
            percent identity (signal-added)
  binary  : the library the program itself builds with `-method slow_pair,lalign_id_pair` (`-lib_only -out_lib`)

Output TSV columns: kind, relax, seqs (|), tree (nested "(a,b)" of 0-based indices), library ("a:p:b:q:w;" 1-based positions as in
the TC_LIB_FORMAT_01 file), rows (|).

usage: generate_tcoffee_oracle.py OUT.tsv [SETS=300] [SEED=1]   (needs `t_coffee` on PATH; takes ~10 minutes for 300 sets)
"""
import os, random, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

M, X, G = 1, -1, -2


def nw(a, b):
    m, n = len(a), len(b)
    F = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1): F[i][0] = i * G
    for j in range(n + 1): F[0][j] = j * G
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            F[i][j] = max(F[i-1][j-1] + (M if a[i-1] == b[j-1] else X), F[i-1][j] + G, F[i][j-1] + G)
    i, j, A, B = m, n, [], []
    while i > 0 or j > 0:
        if i > 0 and j > 0 and F[i][j] == F[i-1][j-1] + (M if a[i-1] == b[j-1] else X):
            A.append(a[i-1]); B.append(b[j-1]); i -= 1; j -= 1
        elif i > 0 and F[i][j] == F[i-1][j] + G:
            A.append(a[i-1]); B.append('-'); i -= 1
        else:
            A.append('-'); B.append(b[j-1]); j -= 1
    return ''.join(reversed(A)), ''.join(reversed(B)), 0, 0


def sw(a, b):
    m, n = len(a), len(b)
    F = [[0] * (n + 1) for _ in range(m + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            F[i][j] = max(0, F[i-1][j-1] + (M if a[i-1] == b[j-1] else X), F[i-1][j] + G, F[i][j-1] + G)
            if F[i][j] > best: best, bi, bj = F[i][j], i, j
    i, j, A, B, si, sj = bi, bj, [], [], None, None
    while i > 0 and j > 0 and F[i][j] > 0:
        si, sj = i, j
        if F[i][j] == F[i-1][j-1] + (M if a[i-1] == b[j-1] else X):
            A.append(a[i-1]); B.append(b[j-1]); i -= 1; j -= 1
        elif F[i][j] == F[i-1][j] + G:
            A.append(a[i-1]); B.append('-'); i -= 1
        else:
            A.append('-'); B.append(b[j-1]); j -= 1
    return ''.join(reversed(A)), ''.join(reversed(B)), (si - 1 if si else 0), (sj - 1 if sj else 0)


def pid(A, B):
    al = sum(1 for x, y in zip(A, B) if x != '-' and y != '-')
    idn = sum(1 for x, y in zip(A, B) if x != '-' and x == y)
    return 0 if al == 0 else (100 * idn) // al


def rand_library(seqs):
    k, lib = len(seqs), {}
    for i in range(k):
        for j in range(i + 1, k):
            for f in (nw, sw):
                A, B, pi, pj = f(seqs[i], seqs[j])
                w = 10 * pid(A, B)
                if w == 0: continue
                for x, y in zip(A, B):
                    if x != '-' and y != '-':
                        lib[(i, pi + 1, j, pj + 1)] = lib.get((i, pi + 1, j, pj + 1), 0) + w; pi += 1; pj += 1
                    elif x != '-': pi += 1
                    elif y != '-': pj += 1
    return lib


def write_fasta(path, seqs):
    with open(path, 'w') as f:
        for n, s in enumerate(seqs): f.write(f'>s{n}\n{s}\n')


def write_lib(path, seqs, lib):
    k = len(seqs)
    with open(path, 'w') as f:
        f.write('! TC_LIB_FORMAT_01\n%d\n' % k)
        for n, s in enumerate(seqs): f.write(f's{n} {len(s)} {s}\n')
        for i in range(k):
            for j in range(i + 1, k):
                f.write(f'# {i+1} {j+1}\n')
                for (a, p, b, q), w in sorted(lib.items()):
                    if a == i and b == j: f.write(f'{p} {q} {w} 0 0\n')
        f.write('! SEQ_1_TO_N\n')


def read_lib(path):
    lib, cur = [], None
    for l in open(path):
        l = l.strip()
        if l.startswith('#'):
            t = l[1:].split(); cur = (int(t[0]) - 1, int(t[1]) - 1)
        elif l and not l.startswith('!') and cur is not None:
            t = l.split()
            if len(t) >= 3 and t[0].isdigit(): lib.append((cur[0], int(t[0]), cur[1], int(t[1]), int(t[2])))
    return lib


def newick(t):
    return f's{t}:0.1' if isinstance(t, int) else f'({newick(t[0])},{newick(t[1])}):0.1'


def plain(t):
    return str(t) if isinstance(t, int) else f'({plain(t[0])},{plain(t[1])})'


def rand_tree(rng, idx):
    idx = list(idx)
    if len(idx) == 1: return idx[0]
    rng.shuffle(idx); m = rng.randint(1, len(idx) - 1)
    return (rand_tree(rng, idx[:m]), rand_tree(rng, idx[m:]))


def run_set(args):
    sid, seed = args
    rng = random.Random(seed * 100003 + sid)
    k = rng.choice([3, 3, 4, 4, 5, 6])
    kind = 'rand' if sid % 2 == 0 else 'binary'
    base = ''.join(rng.choice('ACGT') for _ in range(rng.randint(5, 11)))
    seqs = []
    for _ in range(k):                       # half related (mutated copies of one base), half independent
        if rng.random() < 0.5:
            s = list(base)
            for _ in range(rng.randint(0, 3)):
                p = rng.randrange(len(s)); op = rng.choice('sid')
                if op == 's': s[p] = rng.choice('ACGT')
                elif op == 'i': s.insert(p, rng.choice('ACGT'))
                elif len(s) > 3: del s[p]
            seqs.append(''.join(s))
        else:
            seqs.append(''.join(rng.choice('ACGT') for _ in range(rng.randint(3, 11))))
    tree = rand_tree(rng, range(k))
    with tempfile.TemporaryDirectory() as d:
        fa, lb, dn = d + '/s.fa', d + '/l.tc_lib', d + '/t.dnd'
        write_fasta(fa, seqs)
        env = dict(os.environ, HOME=d)
        if kind == 'rand':
            write_lib(lb, seqs, rand_library(seqs))
        else:
            subprocess.run(['t_coffee', '-seq', fa, '-method', 'slow_pair,lalign_id_pair', '-lib_only', '-out_lib', lb, '-n_core', '1', '-quiet'],
                           capture_output=True, cwd=d, env=env)
            if not os.path.exists(lb): return []
        lib = read_lib(lb)
        open(dn, 'w').write(newick(tree) + ';\n')
        out = []
        for relax in (1, 0):
            of = d + f'/o{relax}.fa'
            subprocess.run(['t_coffee', '-seq', fa, '-lib', lb, '-usetree', dn, '-relax_lib', str(relax), '-outfile', of,
                            '-output', 'fasta_aln', '-quiet', '-n_core', '1'], capture_output=True, cwd=d, env=env)
            if not os.path.exists(of): continue
            rows, cur = {}, None
            for l in open(of):
                l = l.strip()
                if l.startswith('>'): cur = l[1:]; rows[cur] = ''
                elif cur: rows[cur] += l
            out.append('\t'.join([kind, str(relax), '|'.join(seqs), plain(tree),
                                  ';'.join(f'{a}:{p}:{b}:{q}:{w}' for a, p, b, q, w in lib),
                                  '|'.join(rows[f's{i}'] for i in range(k))]))
        return out


def main():
    path = sys.argv[1]
    sets = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 1
    with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as ex:
        results = list(ex.map(run_set, [(i, seed) for i in range(sets)]))
    with open(path, 'w') as f:
        f.write('kind\trelax\tseqs\ttree\tlibrary\trows\n')
        for r in results:
            for line in r: f.write(line + '\n')


if __name__ == '__main__':
    main()
