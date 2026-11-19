#!/usr/bin/env python3
"""B09 audit A2-1: a FRESH random sample of `muscle3` runs (other seeds than the committed oracles) for re-verifying SequenceAligner.Muscle* against the installed Debian binary.

    python3 scripts/alignment/muscle3_fresh_sample.py --seed 7 --count 300 --out /tmp/fresh7.tsv

Needs `muscle3` (apt-get install muscle3, MUSCLE v3.8.1551). Draws `count` random `align` sets (protein / DNA / RNA, 1-9 sequences of 12-150 residues, random option mixtures:
iterations, distance1, cluster1/2, root1/2, weights, objective scores, term gaps, anchors, -brenner, -le/-sp/-sv/-spn, gap parameters, -matrix blosum62 / pam250 / blosum62nostar on every
alphabet), runs them through the package harness (handoff/A1-10c-muscle-r3/run_sets.py) and writes the oracle TSV in the layout of muscle3_r3_oracle.tsv (same 16 columns) to --out.
The same file is what the (uncommitted) comparison test reads; the committed rounds are muscle3_oracle / muscle3_r2_oracle / muscle3_r3_oracle.
"""
import argparse, importlib.util, os, random, re, shutil, subprocess, sys, tempfile

ap = argparse.ArgumentParser()
ap.add_argument('--seed', type=int, required=True)
ap.add_argument('--count', type=int, default=200)
ap.add_argument('--out', required=True)
ap.add_argument('--kinds', default='align', help='comma list: align (default); others reuse the sets of the committed rounds only')
a = ap.parse_args()

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
pkg = os.path.join(root, 'docs', 'Validation', 'review-2026-09', 'handoff', 'A1-10c-muscle-r3')
binary = shutil.which('muscle3') or sys.exit('muscle3 not installed')
ver = subprocess.run([binary, '-version'], capture_output=True, text=True).stdout
assert 'MUSCLE v3.8.1551' in ver, ver
spec = importlib.util.spec_from_file_location('run_sets', os.path.join(pkg, 'run_sets.py'))
rs = importlib.util.module_from_spec(spec); spec.loader.exec_module(rs)
rs.BIN = binary
rs.HERE = pkg    # matrices live in the package

rnd = random.Random(a.seed)
FREQ = dict(A=8.25, R=5.53, N=4.06, D=5.45, C=1.37, Q=3.93, E=6.75, G=7.07, H=2.27, I=5.96, L=9.66, K=5.84, M=2.42, F=3.86, P=4.70, S=6.56, T=5.34, W=1.08, Y=2.92, V=6.87)
AA, WT = list(FREQ), list(FREQ.values())
MIN_LEN = 12


def res(typ):
    return rnd.choices(AA, WT)[0] if typ == 'protein' else rnd.choice('ACGU' if typ == 'rna' else 'ACGT')


def mutate(typ, s, rate):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.55: out.append(res(typ))
        elif r < rate * 0.75: continue
        elif r < rate:
            out.append(ch); out += [res(typ) for _ in range(rnd.choice([1, 1, 2, 4]))]
        else: out.append(ch)
    while len(out) < MIN_LEN: out.append(res(typ))
    return ''.join(out)


def seqs_for(typ):
    k = rnd.choice([1, 2, 2, 3, 3, 4, 5, 6, 7, 9])
    L = rnd.randint(MIN_LEN + 4, 150 if typ != 'protein' else 110)
    base = ''.join(res(typ) for _ in range(L))
    div = rnd.choice([0.03, 0.1, 0.2, 0.4, 0.6])
    out = [mutate(typ, base, rnd.choice([div, div / 2, div * 1.3])) if rnd.random() < 0.92 else ''.join(res(typ) for _ in range(rnd.randint(MIN_LEN, 60))) for _ in range(k)]
    r = rnd.random()
    if typ == 'protein':
        if r < 0.25: out = [('M' + s if rnd.random() < 0.7 else s) for s in out]
        elif r < 0.35: out = [''.join(c if rnd.random() > 0.03 else rnd.choice('BXZ') for c in s) for s in out]
    elif r < 0.3:
        out = [''.join(c if rnd.random() > 0.04 else rnd.choice('NRYXKMSWBDHV') for c in s) for s in out]
    if rnd.random() < 0.1 and typ != 'protein':   # starts with the wildcard M (the amino-acid methionine hack of -matrix on nucleotides)
        out = [('M' + s if rnd.random() < 0.6 else s) for s in out]
    return out


def options(typ):
    d = {}
    if rnd.random() < 0.5: d['maxiters'] = rnd.choice([1, 2, 3, 4, 8])
    if rnd.random() < 0.25: d['distance1'] = rnd.choice(['kmer6_6', 'kmer20_3', 'kmer20_4', 'kbit20_3', 'pctidkimura', 'pctidlog'] + (['kmer4_6'] if typ != 'protein' else []))
    if rnd.random() < 0.2: d['cluster1'] = rnd.choice(['upgmb', 'upgma', 'upgmamin', 'upgmamax', 'neighborjoining'])
    if rnd.random() < 0.15: d['cluster2'] = rnd.choice(['upgmb', 'upgma', 'neighborjoining'])
    if rnd.random() < 0.15: d['root1'] = rnd.choice(['pseudo', 'midlongestspan', 'minavgleafdist'])
    if rnd.random() < 0.15: d['weight1'] = rnd.choice(['none', 'henikoff', 'henikoffpb', 'gsc', 'clustalw', 'threeway'])
    if rnd.random() < 0.15: d['objscore'] = rnd.choice(['sp', 'dp', 'xp', 'spf', 'spm'])
    if rnd.random() < 0.1: d['termgaps'] = rnd.choice(['half', 'full', 'ext'])
    if rnd.random() < 0.1: d['noanchors'] = ''
    if rnd.random() < 0.08: d['brenner'] = ''
    if rnd.random() < 0.2: d[rnd.choice(['le', 'sp', 'sv', 'spn'])] = ''
    if rnd.random() < 0.2: d['gapopen'] = rnd.choice(['-1.5', '-2.9', '-4', '-6', '-200', '-400', '-1000', '-1439'])
    if rnd.random() < 0.1: d['gapextend'] = rnd.choice(['-0.1', '-0.5', '-1', '-10', '-30'])
    if rnd.random() < 0.1: d['center'] = rnd.choice(['-0.3', '-0.52', '-1', '0', '-100'])
    if rnd.random() < 0.35: d['matrix'] = rnd.choice(['@blosum62', '@pam250', '@blosum62nostar'])
    return ','.join(k if v == '' else f'{k}={v}' for k, v in d.items())


matrices = {'@' + n[:-4]: open(os.path.join(pkg, 'inputs', 'matrices', n)).read() for n in os.listdir(os.path.join(pkg, 'inputs', 'matrices'))}
esc = lambda s: s.replace('\\', '\\\\').replace('\t', '\\t').replace('\n', '\\n')
n = 0
with open(a.out, 'w', newline='\n') as fh:
    fh.write('id\tkind\ttype\toptions\tdata\trc\tnames\trows\ttree1\ttree2\tdist1\tdist2\tweights\tusedtree\tstdout\tmatrixtext\n')
    for i in range(a.count):
        typ = rnd.choice(['protein', 'protein', 'dna', 'rna'])
        data = ';'.join(seqs_for(typ))
        opts = options(typ)
        sid = f'F{a.seed}_{i}'
        row = rs.run_one(sid, 'align', typ, data, opts)
        mtext = ''
        for key, text in matrices.items():
            if re.search('=' + re.escape(key) + '(,|$)', opts):
                mtext = text
        fh.write('\t'.join([row[0], row[1], row[2], row[3], data] + row[4:] + [esc(mtext)]) + '\n')
        n += 1
print(n, 'sets ->', a.out)
