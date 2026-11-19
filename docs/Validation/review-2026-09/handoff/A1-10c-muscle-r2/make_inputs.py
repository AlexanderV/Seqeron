#!/usr/bin/env python3
"""A1-10c-muscle ROUND 2 hand-off inputs: deterministic sets for the options and command lines of MUSCLE 3.8.1551 that round 1 (420 sets, equal on all) did not cover.

Writes inputs/sets.tsv with the columns   id  kind  type  data  options
  kind  = align | profile | refine | spscore | maketree | ppscore | sw | refinew | profdb
  type  = protein | dna | rna          (passed as -seqtype)
  data  = align: the unaligned sequences ';'-joined
          profile / ppscore / sw: two alignments  rows1 '|' rows2  (rows ';'-joined, '-' = gap)
          refine / spscore / maketree / refinew: one alignment, rows ';'-joined
          profdb: an alignment '|' unaligned sequences
  options = comma separated  key=value  (or  key  for a flag): every key is passed to muscle as  -key value ; special values:
          @blosum62 / @blosum62nostar / @pam250  -> the path of the matrix file written by run_sets.py (inputs/matrices)
          @own  (usetree)  -> the Newick tree that `muscle -clusteronly -tree1` writes for the same sequences (run first by run_sets.py)
          tree1=1, tree2=1, distmx1=1, distmx2=1, computeweights=1 -> -tree1 / -tree2 / -distmx1 / -distmx2 / -computeweights <file>, the file is returned
Usage: python3 make_inputs.py     (deterministic: the committed inputs/sets.tsv is reproduced exactly)
"""
import os, random

rnd = random.Random(20260901)
HERE = os.path.dirname(os.path.abspath(__file__))
FREQ = dict(A=8.25, R=5.53, N=4.06, D=5.45, C=1.37, Q=3.93, E=6.75, G=7.07, H=2.27, I=5.96, L=9.66, K=5.84, M=2.42, F=3.86, P=4.70, S=6.56, T=5.34, W=1.08, Y=2.92, V=6.87)
AA = list(FREQ)
WT = [FREQ[a] for a in AA]
MIN_LEN = 12     # kmer20_3 / kmer20_4 / kbit20_3 are undefined (they read outside the buffer) for sequences shorter than 2 / 4 / 2 residues


def rprot(): return rnd.choices(AA, WT)[0]


def alpha(typ): return {'protein': AA, 'dna': list('ACGT'), 'rna': list('ACGU')}[typ]


def rres(typ): return rprot() if typ == 'protein' else rnd.choice(alpha(typ))


def mutate(typ, s, rate):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.55: out.append(rres(typ))
        elif r < rate * 0.75: continue
        elif r < rate:
            out.append(ch)
            for _ in range(rnd.choice([1, 1, 1, 2, 3, 5])): out.append(rres(typ))
        else: out.append(ch)
    while len(out) < MIN_LEN: out.append(rres(typ))
    return ''.join(out)


def family(typ, maxlen, k):
    div = rnd.choice([0.02, 0.08, 0.15, 0.3, 0.5, 0.7])
    L = rnd.randint(MIN_LEN + 8, maxlen if typ == 'protein' else maxlen * 2)
    if typ == 'protein':
        kind = rnd.choice([0, 0, 0, 0, 1, 2])
        if kind == 0: base = ''.join(rprot() for _ in range(L))
        elif kind == 1: base = ''.join(rnd.choice('AAAL') for _ in range(L))
        else: base = (''.join(rprot() for _ in range(rnd.randint(2, 6))) * L)[:L]
    else:
        base = ''.join(rres(typ) for _ in range(L))
    seqs = []
    for _ in range(k):
        if rnd.random() < 0.9: seqs.append(mutate(typ, base, rnd.choice([div, div * 0.5, div * 1.3])))
        else: seqs.append(''.join(rres(typ) for _ in range(rnd.randint(MIN_LEN, maxlen))))
    return seqs


def decorate(typ, seqs):
    r = rnd.random()
    if typ == 'protein':
        if r < 0.2: seqs = [('M' + s if rnd.random() < 0.7 else s) for s in seqs]
        elif r < 0.3: seqs = [''.join(c if rnd.random() > 0.03 else rnd.choice('BXZJO') for c in s) for s in seqs]
        elif r < 0.34: seqs = [s.lower() if rnd.random() < 0.5 else s for s in seqs]
    else:
        if r < 0.1: seqs = [''.join(c if rnd.random() > 0.03 else rnd.choice('NRYXKMSW') for c in s) for s in seqs]
    if rnd.random() < 0.05 and len(seqs) >= 2:
        seqs = seqs + [seqs[rnd.randrange(len(seqs))] for _ in range(rnd.randint(1, 2))]
    return seqs


def msa_family(typ, ncols, k, div=None):
    """A true alignment by construction: columns of a base sequence, per-sequence substitutions / deletions and inserted columns that only some sequences fill."""
    div = div if div is not None else rnd.choice([0.05, 0.12, 0.25, 0.4])
    cols = [[rres(typ)] * k for _ in range(ncols)]
    rows = [[] for _ in range(k)]
    for col in cols:
        if rnd.random() < div / 3:     # an inserted column: a few sequences have a residue
            who = set(rnd.sample(range(k), rnd.randint(1, max(1, k // 2))))
            for i in range(k): rows[i].append(rres(typ) if i in who else '-')
        for i in range(k):
            r = rnd.random()
            rows[i].append('-' if r < div / 2 else rres(typ) if r < div else col[i])
    rows = [''.join(r) for r in rows]
    keep = [c for c in range(len(rows[0])) if any(r[c] != '-' for r in rows)]
    rows = [''.join(r[c] for c in keep) for r in rows]
    for i, r in enumerate(rows):                       # every sequence keeps at least MIN_LEN residues
        if len(r.replace('-', '')) < MIN_LEN:
            rows[i] = ''.join(rres(typ) if ch == '-' and rnd.random() < 0.9 else ch for ch in r)
    return rows


def windows_ok(rows, w):
    n = len(rows[0])
    return all(any(r[c] != '-' for c in range(s, min(s + w, n))) for r in rows for s in range(0, n, w))


def opt(d): return ','.join(k if v is None else f'{k}={v}' for k, v in d.items())


def gaps(typ, d):
    if rnd.random() < 0.25: d['gapopen'] = rnd.choice(['-1.5', '-2.9', '-4', '-6'] if typ == 'protein' else ['-200', '-400', '-700'])
    if rnd.random() < 0.12: d['gapextend'] = rnd.choice(['-0.1', '-0.5', '-1'] if typ == 'protein' else ['-10', '-30', '-60'])
    if rnd.random() < 0.1 and typ == 'protein': d['center'] = rnd.choice(['-0.3', '-0.52', '-1'])


WEIGHTS = ['none', 'henikoff', 'henikoffpb', 'gsc', 'clustalw', 'threeway']
OBJ = ['sp', 'dp', 'xp', 'spf', 'spm']
rows_out = []


def add(kind, typ, data, d):
    rows_out.append((f'R{len(rows_out)}', kind, typ, data, opt(d)))


def unaligned(typ, kmin=3, kmax=9, maxlen=None):
    maxlen = maxlen or rnd.choice([40, 70, 70, 130])
    return decorate(typ, family(typ, maxlen, rnd.randint(kmin, kmax)))


def types(n_protein, n_dna, n_rna):
    return ['protein'] * n_protein + ['dna'] * n_dna + ['rna'] * n_rna


def write_matrices():
    """inputs/matrices/*.txt in the NCBI matrix format that readmx.cpp reads (generated once with Biopython's tables; committed, so the user needs no Python packages)."""
    out = os.path.join(HERE, 'inputs', 'matrices')
    os.makedirs(out, exist_ok=True)
    try:
        from Bio.Align import substitution_matrices as sm
    except ImportError:
        print('Biopython not installed: the committed matrix files are kept')
        return
    for name, fn in (('BLOSUM62', 'blosum62'), ('PAM250', 'pam250')):
        m = sm.load(name)
        letters = list(m.alphabet)
        lines = [f'# {name} (Biopython substitution_matrices, NCBI format)', '   ' + '  '.join(f'{c:>2}' for c in letters)]
        for a in letters:
            lines.append(f'{a}  ' + '  '.join(f'{int(m[a, b]):>2}' for b in letters))
        open(os.path.join(out, fn + '.txt'), 'w').write('\n'.join(lines) + '\n')
    # without a trailing '*' column the program reads HeadingCount - 1 columns (readmx.cpp): a 20-residue heading loses the last column
    m = sm.load('BLOSUM62')
    res = list('ARNDCQEGHILKMFPSTWYV')
    lines = ['# BLOSUM62, 20 residues, no * column', '   ' + '  '.join(f'{c:>2}' for c in res)]
    for a in res:
        lines.append(f'{a}  ' + '  '.join(f'{int(m[a, b]):>2}' for b in res))
    open(os.path.join(out, 'blosum62nostar.txt'), 'w').write('\n'.join(lines) + '\n')


def main():
    # A. iterations 3 and above (default 8 when -maxiters is absent), anchors, objective scores, weights, terminal gaps
    for typ in types(75, 30, 10):
        d = {}
        d['maxiters'] = rnd.choice([3, 3, 4, 5, 8, 16]) if rnd.random() < 0.85 else None
        if d['maxiters'] is None: del d['maxiters']
        if rnd.random() < 0.3: d['noanchors'] = None
        elif rnd.random() < 0.25:
            d['anchorspacing'] = rnd.choice([8, 16, 48])
        if rnd.random() < 0.15: d['smoothwindow'] = rnd.choice([3, 5, 9])
        if rnd.random() < 0.1: d['smoothscoreceil'] = rnd.choice(['2', '5', '100', '999'])
        if rnd.random() < 0.1: d['minbestcolscore'] = rnd.choice(['1', '3', '80'])
        if rnd.random() < 0.1: d['minsmoothscore'] = rnd.choice(['0.5', '2', '60'])
        if rnd.random() < 0.4: d['objscore'] = rnd.choice(OBJ + (['ps'] if typ == 'protein' else []))
        if rnd.random() < 0.3: d['weight2'] = rnd.choice(WEIGHTS)
        if rnd.random() < 0.2: d['termgaps'] = rnd.choice(['full', 'half', 'ext'])
        if typ == 'protein' and rnd.random() < 0.15: d[rnd.choice(['sp', 'sv', 'le'])] = None
        gaps(typ, d)
        add('align', typ, ';'.join(unaligned(typ)), d)
    # B. distances of the two trees, clustering and rooting
    for typ in types(50, 20, 10):
        d = {'maxiters': rnd.choice([1, 2, 3, 4])}
        if rnd.random() < 0.6: d['distance1'] = rnd.choice(['kmer20_3', 'kmer20_4', 'kbit20_3', 'pwkimura', 'pwscoredist', 'kmer6_6'] if typ == 'protein' else ['kmer20_3', 'kmer20_4', 'kbit20_3', 'pwkimura', 'kmer4_6'])
        if rnd.random() < 0.5: d['distance2'] = rnd.choice(['pctidkimura', 'pctidlog', 'edit'] + (['scoredist'] if typ == 'protein' else []))
        if rnd.random() < 0.5: d['cluster1'] = rnd.choice(['neighborjoining', 'upgma', 'upgmb', 'upgmamin', 'upgmamax'])
        if rnd.random() < 0.5: d['cluster2'] = rnd.choice(['neighborjoining', 'upgma', 'upgmb', 'upgmamin', 'upgmamax'])
        if rnd.random() < 0.5: d['root1'] = rnd.choice(['pseudo', 'midlongestspan', 'minavgleafdist'])
        if rnd.random() < 0.5: d['root2'] = rnd.choice(['pseudo', 'midlongestspan', 'minavgleafdist'])
        if rnd.random() < 0.15: d['sueff'] = rnd.choice(['0.0', '0.3', '1.0'])
        add('align', typ, ';'.join(unaligned(typ)), d)
    # C. sequence weights (-weight sets both, -weight1 / -weight2 one each) and -maxtrees (more than one round only without the ClustalW weights)
    for typ in types(55, 18, 7):
        d = {'maxiters': rnd.choice([2, 3, 5])}
        r = rnd.random()
        if r < 0.3: d['weight'] = rnd.choice(WEIGHTS)
        else:
            d['weight1'] = rnd.choice(WEIGHTS)
            if r < 0.7: d['weight2'] = rnd.choice(WEIGHTS)
        if rnd.random() < 0.6: d['maxtrees'] = rnd.choice([0, 1, 2, 3, 5])
        add('align', typ, ';'.join(unaligned(typ)), d)
    # D. options that the binary should ignore (compiled out): diagonals, hydrophobicity, double affine gaps, resources, output switches
    ignored = [{'diags': None}, {'diags1': None}, {'diags2': None, 'diaglength': 12}, {'diagbreak': 3, 'diagmargin': 9}, {'hydro': 9}, {'hydro': 3, 'hydrofactor': 2.5}, {'tomhydro': None},
               {'gapopen2': '-5', 'gapextend2': '-1'}, {'gapambig': '0.5'}, {'maxmb': 5}, {'maxhours': 0.5}, {'core': None}, {'nocore': None}, {'group': None}, {'dimer': None}, {'verbose': None}]
    for i, ign in enumerate(ignored * 2):
        typ = ['protein', 'protein', 'dna', 'protein'][i % 4]
        d = {'maxiters': rnd.choice([2, 4, 8])}
        d.update(ign)
        add('align', typ, ';'.join(unaligned(typ)), d)
    # E. -brenner (root alignment), -matrix, -usetree
    for typ in types(14, 4, 2):
        add('align', typ, ';'.join(unaligned(typ)), {'maxiters': rnd.choice([1, 2, 4]), 'brenner': None})
    for i in range(30):
        d = {'maxiters': rnd.choice([1, 2, 3, 6]), 'matrix': rnd.choice(['@blosum62', '@blosum62', '@blosum62nostar', '@pam250'])}
        if rnd.random() < 0.3: d['gapopen'] = rnd.choice(['-5', '-10', '-14'])
        if rnd.random() < 0.2: d[rnd.choice(['le', 'sv', 'sp'])] = None
        add('align', 'protein', ';'.join(unaligned('protein')), d)
    for typ in types(24, 4, 2):
        d = {'usetree': '@own', 'maxiters': rnd.choice([1, 2, 3, 5])}
        if rnd.random() < 0.3: d['weight1'] = rnd.choice(WEIGHTS)
        add('align', typ, ';'.join(unaligned(typ)), d)
    # F. the files of the program: -tree1, -tree2, -distmx1, -distmx2, -computeweights, -clusteronly
    for typ in types(30, 8, 4):
        d = {'maxiters': rnd.choice([1, 2, 3]), 'tree1': 1, 'distmx1': 1}
        if d['maxiters'] >= 2: d['tree2'] = 1; d['distmx2'] = 1
        if rnd.random() < 0.3: d['cluster1'] = 'neighborjoining'
        if rnd.random() < 0.3: d['distance1'] = rnd.choice(['kmer20_4', 'kbit20_3', 'kmer20_3'])
        add('align', typ, ';'.join(unaligned(typ)), d)
    for typ in types(8, 2, 1):
        add('align', typ, ';'.join(unaligned(typ)), {'computeweights': 1, 'maxiters': rnd.choice([1, 2]), **({'weight1': rnd.choice(WEIGHTS)} if rnd.random() < 0.5 else {})})
    for typ in types(5, 2, 1):
        add('align', typ, ';'.join(unaligned(typ)), {'clusteronly': None, 'tree1': 1})
    add('align', 'protein', ';'.join(unaligned('protein')), {'stable': None})          # "-stable not supported in this version of muscle"
    # G. the other command lines
    for typ in types(24, 8, 4):
        k1, k2 = rnd.randint(1, 5), rnd.randint(1, 4)
        n = rnd.randint(25, 70)
        both = msa_family(typ, n, k1 + k2)
        a, b = both[:k1], both[k1:]
        a = [''.join(c for c in r) for r in a]
        # profiles of the same length are not required; drop the all-gap columns of each half
        def trim(rows):
            keep = [c for c in range(len(rows[0])) if any(r[c] != '-' for r in rows)]
            return [''.join(r[c] for c in keep) for r in rows]
        d = {}
        if rnd.random() < 0.3: d['weight1'] = rnd.choice(WEIGHTS)
        if rnd.random() < 0.2: d['weight2'] = rnd.choice(WEIGHTS)
        if rnd.random() < 0.2: d['cluster2'] = rnd.choice(['upgma', 'upgmb', 'neighborjoining'])
        if rnd.random() < 0.15: d['root1'] = rnd.choice(['midlongestspan', 'minavgleafdist'])
        gaps(typ, d)
        add('profile', typ, ';'.join(trim(a)) + '|' + ';'.join(trim(b)), d)
    for typ in types(24, 8, 3):
        rows = msa_family(typ, rnd.randint(30, 90), rnd.randint(3, 9), rnd.choice([0.15, 0.3, 0.45]))
        d = {'maxiters': rnd.choice([1, 2, 4, 8])}
        if rnd.random() < 0.3: d['noanchors'] = None
        if rnd.random() < 0.3: d['objscore'] = rnd.choice(OBJ)
        if rnd.random() < 0.2: d['weight1'] = rnd.choice(WEIGHTS)
        if rnd.random() < 0.2: d['cluster2'] = rnd.choice(['upgma', 'upgmb', 'neighborjoining'])
        add('refine', typ, ';'.join(rows), d)
    for typ in types(14, 5, 3):
        rows = msa_family(typ, rnd.randint(25, 70), rnd.randint(2, 8))
        d = {}
        if rnd.random() < 0.3: d['weight1'] = rnd.choice(WEIGHTS)
        if rnd.random() < 0.25: d['cluster2'] = rnd.choice(['upgma', 'upgmb', 'neighborjoining'])
        if rnd.random() < 0.25: d['distance2'] = rnd.choice(['pctidlog', 'edit'])
        gaps(typ, d)
        add('spscore', typ, ';'.join(rows), d)
    for typ in types(14, 6, 3):
        rows = msa_family(typ, rnd.randint(25, 70), rnd.randint(3, 9))
        d = {}
        if rnd.random() < 0.4: d['cluster2'] = rnd.choice(['upgma', 'upgmb', 'upgmamin', 'upgmamax', 'neighborjoining'])
        if rnd.random() < 0.3: d['root2'] = rnd.choice(['midlongestspan', 'minavgleafdist'])
        if rnd.random() < 0.3: d['distance2'] = rnd.choice(['pctidlog', 'edit'])
        add('maketree', typ, ';'.join(rows), d)
    for typ in types(5, 16, 6):                         # protein: expected to crash (no score matrix), recorded
        both = msa_family(typ, rnd.randint(25, 60), rnd.randint(2, 7))
        k = len(both) // 2
        d = {}
        if rnd.random() < 0.3: d['weight1'] = rnd.choice(WEIGHTS)
        add('ppscore', typ, ';'.join(both[:k]) + '|' + ';'.join(both[k:2 * k]), d)
    for typ in types(3, 10, 5):
        for _ in range(100):
            rows = msa_family(typ, rnd.randint(40, 90), rnd.randint(3, 8), 0.12)
            w = rnd.choice([15, 20, 30])
            if windows_ok(rows, w): break
        d = {'refinewindow': w}
        if rnd.random() < 0.4: d['maxiters'] = rnd.choice([1, 2, 3])
        add('refinew', typ, ';'.join(rows), d)
    for typ in types(3, 14, 6):
        a = family(typ, 60, 1)[0]
        b = mutate(typ, a[rnd.randint(0, 10):], 0.15)
        multi = rnd.random() < 0.35
        d = {}
        if multi:
            d['distance1'] = 'pctidkimura'
            rows = msa_family(typ, rnd.randint(25, 45), 3)
            rows2 = msa_family(typ, rnd.randint(25, 45), 2)
            add('sw', typ, ';'.join(rows) + '|' + ';'.join(rows2), d)
        else:
            add('sw', typ, a + '|' + b, d)
    for typ in ['protein', 'dna', 'protein', 'dna', 'dna']:
        rows = msa_family(typ, rnd.randint(25, 45), 3)
        add('profdb', typ, ';'.join(rows) + '|' + ';'.join(unaligned(typ, 2, 3, 40)), {})
    os.makedirs(os.path.join(HERE, 'inputs'), exist_ok=True)
    write_matrices()
    with open(os.path.join(HERE, 'inputs', 'sets.tsv'), 'w') as f:
        f.write('id\tkind\ttype\tdata\toptions\n')
        for r in rows_out:
            f.write('\t'.join(r) + '\n')
    from collections import Counter
    print(len(rows_out), 'sets:', dict(Counter(r[1] for r in rows_out)))


if __name__ == '__main__':
    main()
