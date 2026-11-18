#!/usr/bin/env python3
"""Oracle generator for SequenceAligner.MultipleAlignTCoffeeProgram (audit A1-13b, B09 F66).

Runs the installed T-Coffee program (Debian package t-coffee 13.45.0.4846264+really13.41.0.28bdc39+dfsg-1, `t_coffee`) on random DNA and
PROTEIN sets and records, for the program's OWN primary-library methods (`proba_pair` = the default of `t_coffee seqs.fa`, `slow_pair`,
`lalign_id_pair` and their combinations):

  tcoffee_program_library_oracle.tsv : small sets (2-4 sequences), the primary library `t_coffee -seq s.fa -method M -lib_only -out_lib` writes
                                       (every printed pair "a:p:b:q:weight:count", 1-based residues, s1 < s2 side as in the file).
  tcoffee_program_msa_oracle.tsv     : sets of 2-8 sequences, the guide tree the program used (parsed from its "Group n: [Group r (k seq)] with
                                       [Group l (k seq)]" log lines, nested "(left,right)" of 0-based input indices) and the alignment
                                       `t_coffee -seq s.fa [-method M] -output fasta_aln` writes (rows in input order).

  tcoffee_program_big_oracle.tsv     : (only with BIG=1) two sets of 100 and 101 sequences: the program switches to its fast neighbour joining
                                       (`fast_nj_tree`) from 100 sequences on.

Sequence names are q000, q001, ... : the program mis-handles names of which one is a prefix of another (s1 / s10: its name lookup compares prefixes), which
changes the guide tree for 11 or more sequences named s0, s1, ...; the port assumes distinct names without that relation.
Everything runs with `-n_core 1 -no_warning`; the library does not depend on -n_core (checked with 2, 3, 4 and 8 cores).

usage: [BIG=1] generate_tcoffee_program_oracle.py OUT_DIR [SCALE=1.0] [SEED=7]   (needs `t_coffee` on PATH; ~25 minutes for SCALE 1)
"""
import os, random, re, subprocess, sys, tempfile
from concurrent.futures import ThreadPoolExecutor

DNA = 'ACGT'
PROTEIN = 'ACDEFGHIKLMNPQRSTVWY'

# (kind, alphabet, method list (None = plain t_coffee), sets, k range, length range)
LIBRARY_PLAN = [
    ('dna', DNA, 'slow_pair', 90, (2, 4), (3, 14)),
    ('dna', DNA, 'lalign_id_pair', 90, (2, 4), (3, 14)),
    ('dna', DNA, 'proba_pair', 90, (2, 4), (3, 14)),
    ('dna', DNA, 'slow_pair,lalign_id_pair', 120, (2, 4), (3, 12)),
    ('dna', DNA, 'lalign_id_pair,slow_pair', 60, (2, 4), (3, 12)),
    ('dna', DNA, 'proba_pair,slow_pair,lalign_id_pair', 60, (2, 3), (3, 12)),
    ('protein', PROTEIN, 'slow_pair', 70, (2, 4), (3, 14)),
    ('protein', PROTEIN, 'lalign_id_pair', 70, (2, 4), (3, 14)),
    ('protein', PROTEIN, 'proba_pair', 70, (2, 4), (3, 14)),
    ('protein', PROTEIN, 'slow_pair,lalign_id_pair', 100, (2, 4), (3, 12)),
]
MSA_PLAN = [
    ('dna', DNA, 'plain', 260, (2, 8), (3, 36)),
    ('protein', PROTEIN, 'plain', 260, (2, 8), (3, 36)),
    ('dna', DNA, 'slow_pair,lalign_id_pair', 260, (2, 8), (3, 36)),
    ('protein', PROTEIN, 'slow_pair,lalign_id_pair', 260, (2, 8), (3, 36)),
    ('dna', DNA, 'slow_pair', 80, (2, 8), (3, 30)),
    ('dna', DNA, 'lalign_id_pair', 80, (2, 8), (3, 30)),
    ('dna', DNA, 'lalign_id_pair,slow_pair', 80, (2, 7), (3, 30)),
    ('dna', DNA, 'proba_pair,slow_pair,lalign_id_pair', 80, (2, 7), (3, 30)),
    ('protein', PROTEIN, 'slow_pair', 60, (2, 7), (3, 30)),
    ('protein', PROTEIN, 'lalign_id_pair', 60, (2, 7), (3, 30)),
]


def make_set(rng, alphabet, krange, lrange):
    k = rng.randint(*krange)
    base = ''.join(rng.choice(alphabet) for _ in range(rng.randint(*lrange)))
    seqs = []
    for _ in range(k):                       # 60 % related (mutated copies of one base), the others independent
        if rng.random() < 0.6:
            s = list(base)
            for _ in range(rng.randint(0, max(1, len(base) // 6))):
                p = rng.randrange(len(s)); op = rng.choice('sid')
                if op == 's': s[p] = rng.choice(alphabet)
                elif op == 'i': s.insert(p, rng.choice(alphabet))
                elif len(s) > 2: del s[p]
            seqs.append(''.join(s))
        else:
            seqs.append(''.join(rng.choice(alphabet) for _ in range(rng.randint(max(2, lrange[0] - 1), lrange[1]))))
    return seqs


def run_tc(seqs, method, want_lib):
    with tempfile.TemporaryDirectory() as d:
        fa = d + '/s.fa'
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>q{i:03d}\n{s}\n')
        env = dict(os.environ, HOME=d)
        cmd = ['t_coffee', '-seq', fa, '-n_core', '1', '-no_warning']
        if method != 'plain': cmd += ['-method', method]
        if want_lib:
            cmd += ['-lib_only', '-out_lib', d + '/o.lib']
            subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env)
            lib, cur = [], None
            try:
                for l in open(d + '/o.lib'):
                    l = l.strip()
                    if l.startswith('#'):
                        t = l[1:].split(); cur = (int(t[0]), int(t[1]))
                    elif l and not l.startswith('!') and cur is not None:
                        t = l.split()
                        if len(t) >= 4 and t[0].isdigit(): lib.append(f'{cur[0]-1}:{t[0]}:{cur[1]-1}:{t[1]}:{t[2]}:{t[3]}')
            except FileNotFoundError:
                return None
            return ';'.join(lib)
        cmd += ['-outfile', d + '/o.fa', '-output', 'fasta_aln']
        r = subprocess.run(cmd, capture_output=True, text=True, cwd=d, env=env)
        try:
            rows, cur = {}, None
            for l in open(d + '/o.fa'):
                l = l.strip()
                if l.startswith('>'): cur = l[1:]; rows[cur] = ''
                elif cur: rows[cur] += l
            rows = [rows[f'q{i:03d}'] for i in range(len(seqs))]
        except FileNotFoundError:
            return None
        groups = re.findall(r'Group\s+(\d+): \[Group\s+(\d+) \(\s*\d+ seq\)\] with \[Group\s+(\d+) \(\s*\d+ seq\)\]', (r.stderr + r.stdout).replace('\r', '\n'))
        node = {i: str(i - 1) for i in range(1, len(seqs) + 1)}
        for p, rg, lf in groups: node[int(p)] = '(' + node[int(lf)] + ',' + node[int(rg)] + ')'
        tree = node[int(groups[-1][0])] if groups else ('0' if len(seqs) == 1 else '')
        return tree, '|'.join(rows)


def job(args):
    kind, alphabet, method, idx, krange, lrange, seed, want_lib = args
    rng = random.Random(seed * 1000003 + idx * 7919 + sum(map(ord, kind + method)) * 104729 + (0 if want_lib else 17))
    seqs = make_set(rng, alphabet, krange, lrange)
    res = run_tc(seqs, method, want_lib)
    if res is None: return None
    if want_lib: return '\t'.join([kind, method, '|'.join(seqs), res])
    tree, rows = res
    return '\t'.join([kind, method, '|'.join(seqs), tree, rows])


BIG_PLAN = [('dna', DNA, 'slow_pair', 1, (100, 100), (6, 10)), ('dna', DNA, 'proba_pair', 1, (101, 101), (6, 10))]


def main():
    out = sys.argv[1]
    if os.environ.get('BIG'):
        seed = int(sys.argv[3]) if len(sys.argv) > 3 else 7
        os.makedirs(out, exist_ok=True)
        jobs = [(kind, alphabet, method, i, kr, lr, seed, False) for kind, alphabet, method, n, kr, lr in BIG_PLAN for i in range(n)]
        with open(os.path.join(out, 'tcoffee_program_big_oracle.tsv'), 'w') as f:
            f.write('type\tmethods\tseqs\ttree\trows\n')
            for j in jobs:
                r = job(j)
                if r: f.write(r + '\n')
        return
    scale = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 7
    os.makedirs(out, exist_ok=True)
    for want_lib, plan, name, header in ((True, LIBRARY_PLAN, 'tcoffee_program_library_oracle.tsv', 'type\tmethods\tseqs\tlibrary'),
                                         (False, MSA_PLAN, 'tcoffee_program_msa_oracle.tsv', 'type\tmethods\tseqs\ttree\trows')):
        jobs = [(kind, alphabet, method, i, kr, lr, seed, want_lib) for kind, alphabet, method, n, kr, lr in plan for i in range(max(1, int(n * scale)))]
        with ThreadPoolExecutor(max_workers=os.cpu_count() or 4) as ex:
            results = list(ex.map(job, jobs))
        with open(os.path.join(out, name), 'w') as f:
            f.write(header + '\n')
            for r in results:
                if r: f.write(r + '\n')


if __name__ == '__main__':
    main()
