#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-10c part b) oracle: ClustalW 2.1 multiple-alignment stage for PROTEIN (`-TYPE=PROTEIN -ALIGN`).

Seqeron API: SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWProteinOptions {...}).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7) `-TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` with the option variants
-MATRIX={GONNET,BLOSUM,PAM,ID} -PWMATRIX={GONNET,BLOSUM,PAM,ID} -GAPOPEN -GAPEXT -PWGAPOPEN -PWGAPEXT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE
(-KTUPLE -WINDOW -TOPDIAGS -PAIRGAP) -NOPGAP -NOHGAP -NOVGAP -GAPDIST -HGAPRESIDUES -ENDGAPS (a no-op on the command line: it only writes
useEndGaps = false, so the rows must not change). The binary computes its own pairwise distances, neighbour-joining
guide tree, sequence weights and profile alignments (matrix series chosen by the percent identity of the two profiles, residue-specific / hydrophilic /
variable-region gap penalties, gap-separation distance); the C# port must reproduce the FASTA rows exactly (input order).

Output TSV (tests/.../TestData/Alignment/clustalw_protein_msa_oracle.tsv):  id  seqs(;)  rows(;)  options(k=v,k=v)
Second output (python3 ... --series OUTDIR MATRICES_H EMBOSS_DATA_DIR): clustalw_protein_series_reference.tsv  matrix a b value - the canonical references of the
matrix SERIES tables added to SubstitutionMatrix (ClustalWBlosum40/45/80, ClustalWPam20/60/120/350): every one of the 23 x 23 cells of matrices.h blosum40mt,
blosum45mt, blosum80mt, pam20mt, pam60mt, pam120mt, pam350mt equals the EMBOSS data file EBLOSUM40/45/80, EPAM20/60/120/350 (asserted here); blosum62mt2 is
asserted to be exactly 2 x EBLOSUM62. (GONNET40/80/120/160/350 have no independent local reference: matrices.h only; GONNET250 is checked by F62.)

Usage: python3 scripts/alignment/generate_clustalw_protein_msa_oracle.py OUTDIR [nSets] [seed]      (needs `clustalw` on PATH)
Committed file: OUTDIR/clustalw_protein_msa_oracle.tsv from `... OUTDIR 1000` (default seed).
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261126)
# Swiss-Prot background frequencies (percent, rounded) of the 20 amino acids
FREQ = dict(A=8.25, R=5.53, N=4.06, D=5.45, C=1.37, Q=3.93, E=6.75, G=7.07, H=2.27, I=5.96, L=9.66, K=5.84, M=2.42, F=3.86, P=4.70, S=6.56, T=5.34, W=1.08, Y=2.92, V=6.87)
AA = list(FREQ)
W = [FREQ[a] for a in AA]


def rres(): return rnd.choices(AA, W)[0]


def gen(kind, L):
    if kind == 0: return ''.join(rres() for _ in range(L))
    if kind == 1: return ''.join(rnd.choice("GPSNDQEKR") for _ in range(L))             # hydrophilic run material
    if kind == 2: return ''.join(rnd.choice("AAAL") for _ in range(L))                   # low complexity
    if kind == 3:
        u = ''.join(rres() for _ in range(rnd.randint(1, 5)))
        return (u * L)[:L]
    # blocks of hydrophilic stretches in random sequence (exercises -HGAPRESIDUES runs of > 3)
    out = []
    while len(out) < L:
        if rnd.random() < 0.3: out += [rnd.choice("GPSNDQEKR") for _ in range(rnd.randint(2, 9))]
        else: out += [rres() for _ in range(rnd.randint(3, 25))]
    return ''.join(out[:L])


def mutate(s, rate):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.55: out.append(rres())
        elif r < rate * 0.75: continue
        elif r < rate:
            out.append(ch)
            for _ in range(rnd.choice([1, 1, 1, 2, 3, 5])): out.append(rres())
        else: out.append(ch)
    return ''.join(out) or 'A'


def make_set(maxlen, kmin, kmax):
    k = rnd.randint(kmin, kmax)
    base = gen(rnd.choice([0, 0, 0, 0, 1, 2, 3, 4, 4]), rnd.randint(10, maxlen))
    seqs = []
    div = rnd.choice([0.02, 0.08, 0.15, 0.3, 0.5, 0.7, 0.85])
    for _ in range(k):
        if rnd.random() < 0.85: seqs.append(mutate(base, rnd.choice([div, div * 0.5, div * 1.3])))
        else: seqs.append(gen(rnd.choice([0, 0, 1, 2, 3, 4]), rnd.randint(5, maxlen)))
    return seqs


def options():
    o = {}
    if rnd.random() < 0.3: return o
    if rnd.random() < 0.3: o['matrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
    if rnd.random() < 0.15: o['pwmatrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
    if rnd.random() < 0.3: o['gapopen'] = rnd.choice([2, 5, 8, 10, 12.5, 15, 20, 30])
    if rnd.random() < 0.3: o['gapext'] = rnd.choice([0.05, 0.1, 0.2, 0.5, 1, 2, 3.3])
    if rnd.random() < 0.2: o['maxdiv'] = rnd.choice([0, 10, 30, 50, 70, 100])
    if rnd.random() < 0.1: o['noweights'] = None
    if rnd.random() < 0.15: o['negative'] = None
    if rnd.random() < 0.2: o['nopgap'] = None
    if rnd.random() < 0.2: o['nohgap'] = None
    if rnd.random() < 0.25: o['novgap'] = None
    if rnd.random() < 0.25: o['gapdist'] = rnd.choice([0, 1, 2, 4, 6, 8, 12, 20])
    if rnd.random() < 0.08: o['endgaps'] = None
    if rnd.random() < 0.15: o['hgapresidues'] = ''.join(rnd.sample("ACDEFGHIKLMNPQRSTVWY", rnd.randint(1, 12)))
    if rnd.random() < 0.2:
        o['quick'] = None
        if rnd.random() < 0.5:
            k = rnd.choice([1, 2]); o['ktuple'] = k
            if rnd.random() < 0.5: o['window'] = rnd.randint(k + 1, k + 8)
            if rnd.random() < 0.5: o['topdiags'] = rnd.randint(k + 1, k + 8)
            if rnd.random() < 0.5: o['pairgap'] = rnd.randint(k + 1, k + 10)
    if rnd.random() < 0.12: o['pwgapopen'] = rnd.choice([3, 5, 10, 15, 20])
    if rnd.random() < 0.12: o['pwgapext'] = rnd.choice([0.05, 0.1, 0.5, 1, 2])
    return o


FLAGS = {'gapopen': 'GAPOPEN', 'gapext': 'GAPEXT', 'maxdiv': 'MAXDIV', 'noweights': 'NOWEIGHTS', 'negative': 'NEGATIVE', 'quick': 'QUICKTREE',
         'ktuple': 'KTUPLE', 'window': 'WINDOW', 'topdiags': 'TOPDIAGS', 'pairgap': 'PAIRGAP', 'pwgapopen': 'PWGAPOPEN', 'pwgapext': 'PWGAPEXT',
         'matrix': 'MATRIX', 'pwmatrix': 'PWMATRIX', 'nopgap': 'NOPGAP', 'nohgap': 'NOHGAP', 'novgap': 'NOVGAP', 'gapdist': 'GAPDIST',
         'hgapresidues': 'HGAPRESIDUES', 'endgaps': 'ENDGAPS'}


def run(seqs, o):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
        out = os.path.join(d, 'o.fa')
        args = ['clustalw', f'-INFILE={fa}', '-TYPE=PROTEIN', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}']
        for k, v in o.items():
            args.append(f'-{FLAGS[k]}' + ('' if v is None else f'={v}'))
        subprocess.run(args, capture_output=True, text=True, cwd=d)
        rows = {}
        name = None
        for line in open(out):
            line = line.rstrip('\n')
            if line.startswith('>'): name = line[1:].strip(); rows[name] = ''
            elif name: rows[name] += line.strip()
    return [rows[f's{i}'] for i in range(len(seqs))]


def ambiguous(seqs):
    return [''.join(c if rnd.random() > 0.04 else rnd.choice("BZXOU") for c in s) for s in seqs]


def main(outdir, nsets):
    os.makedirs(outdir, exist_ok=True)
    kinds = [(0.5, 70), (0.3, 130), (0.15, 220), (0.05, 320)]
    with open(os.path.join(outdir, 'clustalw_protein_msa_oracle.tsv'), 'w') as f:
        f.write('id\tseqs\trows\toptions\n')
        for i in range(nsets):
            r = rnd.random(); acc = 0
            for p, ml in kinds:
                acc += p
                if r < acc: break
            seqs = make_set(ml, 2, rnd.choice([5, 8, 8, 12, 25]))
            if rnd.random() < 0.1: seqs = ambiguous(seqs)
            if rnd.random() < 0.05: seqs = seqs + [seqs[rnd.randrange(len(seqs))] for _ in range(rnd.randint(1, 3))]
            if rnd.random() < 0.04: seqs = [gen(rnd.randint(0, 3), rnd.randint(1, 6)) for _ in range(rnd.randint(2, 6))]
            o = options()
            rows = run(seqs, o)
            assert all(len(x) == len(rows[0]) for x in rows), (i, o)
            ostr = ','.join(k if v is None else f'{k}={v}' for k, v in o.items())
            f.write(f"{i}\t{';'.join(seqs)}\t{';'.join(rows)}\t{ostr}\n")
    print(nsets, 'sets')


SERIES = [('BLOSUM40', 'blosum40mt', 'EBLOSUM40', 1), ('BLOSUM45', 'blosum45mt', 'EBLOSUM45', 1), ('BLOSUM80', 'blosum80mt', 'EBLOSUM80', 1),
          ('PAM20', 'pam20mt', 'EPAM20', 1), ('PAM60', 'pam60mt', 'EPAM60', 1), ('PAM120', 'pam120mt', 'EPAM120', 1), ('PAM350', 'pam350mt', 'EPAM350', 1),
          ('BLOSUM62', 'blosum62mt2', 'EBLOSUM62', 2)]
ORDER = "ABCDEFGHIKLMNPQRSTVWXYZ"


def series_reference(outdir, matrices_h, emboss_dir):
    import re
    text = open(matrices_h).read()
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'clustalw_protein_series_reference.tsv'), 'w') as f:
        f.write('matrix\ta\tb\tvalue\n')
        for name, table, emb, mult in SERIES:
            body = re.sub(r'/\*.*?\*/', '', re.search(r'short\s+' + table + r'\[\]\s*=\s*\{([^}]*)\}', text, re.S).group(1), flags=re.S)
            tri = [int(x) for x in re.findall(r'-?\d+', body)]
            assert len(tri) == 276
            rows = [l.split() for l in open(os.path.join(emboss_dir, emb)) if not l.startswith('#') and l.strip()]
            cols = rows[0]
            E = {(r[0], c): int(r[1 + k]) for r in rows[1:] for k, c in enumerate(cols)}
            ix = 0
            for i in range(23):
                for j in range(i + 1):
                    a, b = ORDER[i], ORDER[j]
                    assert tri[ix] == E[(a, b)] * mult == E[(b, a)] * mult, (name, a, b)
                    ix += 1
            for a in ORDER:
                for b in ORDER:
                    f.write(f'{name}\t{a}\t{b}\t{E[(a, b)]}\n')
    print('series reference written')


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == '--series':
        series_reference(sys.argv[2], sys.argv[3], sys.argv[4]); sys.exit(0)
    if len(sys.argv) > 3: rnd.seed(int(sys.argv[3]))
    main(sys.argv[1] if len(sys.argv) > 1 else '.', int(sys.argv[2]) if len(sys.argv) > 2 else 700)
