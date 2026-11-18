#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-15b part b) oracle: ClustalW 2.1 PROTEIN pairwise stages (`-TYPE=PROTEIN`).

Seqeron API: SequenceAligner.ClustalWPairwiseAlign(a, b, ClustalWProteinOptions) (full stage), SequenceAligner.ClustalWQuickPairwiseAlign(a, b,
ClustalWSequenceType.Protein, ClustalWQuickOptions, maxSequenceLength) (fast stage) and MsaOptions { Distance = ClustalWAffine | ClustalWQuick } on
protein inputs.

References (all run here):
  * `clustalw` 2.1 (Debian 2.1+lgpl-7) `-TYPE=PROTEIN -ALIGN [-PWMATRIX=GONNET|BLOSUM|PAM|ID] [-PWGAPOPEN -PWGAPEXT]` and `-QUICKTREE [-KTUPLE -WINDOW
    -TOPDIAGS -PAIRGAP -SCORE]`: "Sequences (i:j) Aligned. Score: N" (full: (int) of the single-precision percent identity of the Smith-Waterman-
    Gotoh / Myers-Miller alignment; quick: %lg of the double percent fragment score, "Not Aligned" for <= 0.1).
  * The full stage is also reproduced by `CW` of generate_clustalw_pairwise_oracle.py (line-by-line port of FullPairwiseAlign.cpp, matrix and gap
    penalties replaced by the protein ones computed here from SubMatrix.cpp::getMatrix / getPairwiseMatrix: matrices.h blosum30mt / pam350mt / gon250mt /
    idmat times intScale 100 (GONNET: 10), matrixAvgScore = -int(sum of the off-diagonal cells over the 25 residue codes / ((22^2 - 22) / 2.0f)), gap open
    = (int)(2 * avg * (PWGAPOPEN + ln MIN(n, m)) * 1.0) (2 * (int)((PWGAPOPEN + ln MIN) * intScale) when avg <= 0), gap extend = (int)(PWGAPEXT * intScale)
    in single precision). Every row is asserted equal to the binary before it is written.
  * matrices.h: sources.debian.org clustalw 2.1+lgpl-7 src/substitutionMatrix/matrices.h (path via env CLUSTALW_MATRICES_H; the four triangles are
    committed in the C# class, this script re-derives them).

Output TSV files (tests/.../TestData/Alignment/):
  clustalw_protein_pairwise_oracle.tsv: matrix(GONNET|BLOSUM|PAM|ID) gapOpen gapExtend seq1 seq2 clustalwScore identities dpScore gapOpenInt gapExtendInt matAvg
  clustalw_protein_quick_oracle.tsv   : kind ktuple window topdiags pairgap percent maxLen seq1 seq2 printed
  msa_matrix_clustalw_oracle.tsv      : model matrix open ext pwMatrix pwOpen pwExt seqs(;) progressive(;) iterative(;) consistency(;)
        protein MSA rows (independent Python of generate_msa_matrix_oracle.py) whose guide tree is built from the ClustalW PROTEIN pairwise distances
        (1 - identities / MIN(len) of the port above, exact rationals, UPGMA)  == MsaOptions { Distance = ClustalWAffine, ClustalWProtein = ... }
  clustalw_protein_matrices_reference.tsv: matrix a b value - the canonical references of the matrix tables: BLOSUM30 = EMBOSS EBLOSUM30 (all 23 x 23
        cells), GONNET250 = round(10 x Biopython 1.88 GONNET1992) on the 20 amino acids (the two tables are asserted equal to matrices.h here)
Usage: python3 scripts/alignment/generate_clustalw_protein_pairwise_oracle.py OUTDIR MATRICES_H [EBLOSUM30]     (needs `clustalw` on PATH, Biopython)
"""
import importlib.util, math, os, random, re, subprocess, sys, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('cwp', os.path.join(HERE, 'generate_clustalw_pairwise_oracle.py'))
cwp = importlib.util.module_from_spec(spec); spec.loader.exec_module(cwp)

CODES = "ABCDEFGHIKLMNOPQRSTUVWXYZ-"     # UserParameters::aminoAcidCodes
AA_ORDER = "ABCDEFGHIKLMNPQRSTVWXYZ"     # SubMatrix::setUpCrossReferences aminoAcidOrder (the packing order of matrices.h)
MAX_AA = len(CODES) - 2                  # 24
NAMES = {'GONNET': 'gon250mt', 'BLOSUM': 'blosum30mt', 'PAM': 'pam350mt', 'ID': 'idmat'}


def parse_triangles(path):
    text = open(path).read()
    out = {}
    for name in NAMES.values():
        m = re.search(r'short\s+' + name + r'\[\]\s*=\s*\{([^}]*)\}', text)
        out[name] = [int(x) for x in re.findall(r'-?\d+', m.group(1))]
        assert len(out[name]) == 276, name
    return out


def pairwise_matrix(tri, matrix):
    """SubMatrix::getMatrix for the pairwise matrices (negFlag = true) -> (matrix over residue codes, intScale, matrixAvgScore)."""
    scale = 100 // 10 if matrix == 'GONNET' else 100
    n = len(CODES)
    xref = [-1] * (MAX_AA + 2)
    for i, c in enumerate(AA_ORDER): xref[i] = CODES.index(c)
    M = [[0] * n for _ in range(n)]
    ix = 0; max_res = 0
    for i in range(MAX_AA + 1):
        ti = xref[i]
        for j in range(i + 1):
            tj = xref[j]
            if ti != -1 and tj != -1:
                k = tri[NAMES[matrix]][ix]
                if ti == tj:
                    M[ti][ti] = k * scale; max_res += 1
                else:
                    M[ti][tj] = k * scale; M[tj][ti] = k * scale
                ix += 1
    max_res -= 1
    av3 = 0
    for i in range(MAX_AA + 1):
        for j in range(i + 1):
            if i != j: av3 += M[i][j]
    av3 = int(np.float32(av3) / (np.float32(max_res * max_res - max_res) / np.float32(2)))
    return M, scale, -av3


def gap_penalties(n, m, mat_avg, scale, pw_open, pw_ext):
    go32, ge32 = np.float32(pw_open), np.float32(pw_ext)
    ml = math.log(float(min(n, m)))
    if mat_avg <= 0:
        go = 2 * int((float(go32) + ml) * scale)
    else:
        go = int(2 * mat_avg * (float(go32) + ml) * 1.0)
    return go, int(np.float32(ge32 * np.float32(scale)))


rnd = random.Random(20261125)
AA20 = "ARNDCQEGHILKMFPSTWYV"
AA23 = "ABCDEFGHIKLMNPQRSTVWXYZ"


def gen(kind, L, alpha):
    if kind == 0: return ''.join(rnd.choice(alpha) for _ in range(L))
    if kind == 1: return ''.join(rnd.choice(alpha[:2]) for _ in range(L))
    if kind == 2: return ''.join(rnd.choice([alpha[0]] * 3 + [alpha[1]]) for _ in range(L))
    u = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 4)))
    return (u * L)[:L]


def mutate(s, rate, alpha):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: continue
        elif r < rate: out.append(ch); out.append(rnd.choice(alpha))
        else: out.append(ch)
    return ''.join(out) or alpha[0]


def make_set(maxlen, alpha, kmin=3, kmax=7):
    k = rnd.randint(kmin, kmax)
    base = gen(rnd.choice([0, 0, 0, 1, 2, 3]), rnd.randint(5, maxlen), alpha)
    seqs = []
    for _ in range(k):
        if rnd.random() < 0.6: seqs.append(mutate(base, rnd.choice([0.05, 0.15, 0.3, 0.5]), alpha))
        else: seqs.append(gen(rnd.choice([0, 0, 1, 2, 3]), rnd.randint(3, maxlen), alpha))
    return seqs


def run_full(seqs, matrix, pw_open, pw_ext):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
        args = ['clustalw', f'-INFILE={fa}', '-TYPE=PROTEIN', '-ALIGN', f'-OUTFILE={d}/o.aln']
        if matrix: args.append(f'-PWMATRIX={matrix}')
        if pw_open is not None: args.append(f'-PWGAPOPEN={pw_open}')
        if pw_ext is not None: args.append(f'-PWGAPEXT={pw_ext}')
        out = subprocess.run(args, capture_output=True, text=True, cwd=d).stdout
    res = {(int(a) - 1, int(b) - 1): int(c) for a, b, c in re.findall(r'Sequences \((\d+):(\d+)\) Aligned\. Score:\s+(-?\d+)', out)}
    assert len(res) == len(seqs) * (len(seqs) - 1) // 2, out[-400:]
    return res


def run_quick(seqs, p):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
        args = ['clustalw', f'-INFILE={fa}', '-TYPE=PROTEIN', '-QUICKTREE', '-ALIGN', f'-OUTFILE={d}/o.aln']
        if p['ktuple']: args.append(f"-KTUPLE={p['ktuple']}")
        if p['window']: args.append(f"-WINDOW={p['window']}")
        if p['topdiags']: args.append(f"-TOPDIAGS={p['topdiags']}")
        if p['pairgap']: args.append(f"-PAIRGAP={p['pairgap']}")
        if not p['percent']: args.append('-SCORE=ABSOLUTE')
        out = subprocess.run(args, capture_output=True, text=True, cwd=d).stdout
    res = {}
    for a, b, sc in re.findall(r'Sequences \((\d+):(\d+)\) Aligned\. Score:\s+(\S+)', out): res[(int(a) - 1, int(b) - 1)] = sc
    for a, b in re.findall(r'Sequences \((\d+):(\d+)\) Not Aligned', out): res[(int(a) - 1, int(b) - 1)] = 'NA'
    assert len(res) == len(seqs) * (len(seqs) - 1) // 2, out[-500:]
    return res


def quick_params():
    r = rnd.random()
    if r < 0.4: return dict(ktuple=None, window=None, topdiags=None, pairgap=None, percent=True)
    kt = rnd.choice([None, 1, 2])
    k = kt or 1
    return dict(ktuple=kt,
                window=rnd.choice([None, rnd.randint(k + 1, k + 9)]),
                topdiags=rnd.choice([None, rnd.randint(k + 1, k + 10)]),
                pairgap=rnd.choice([None, rnd.randint(k + 1, k + 9)]),
                percent=rnd.random() < 0.8)


def effective(p):
    k = p['ktuple'] or 1
    return dict(ktuple=k,
                window=p['window'] if p['window'] and p['window'] > k else 5,
                topdiags=p['topdiags'] if p['topdiags'] and p['topdiags'] > k else 5,
                pairgap=p['pairgap'] if p['pairgap'] and p['pairgap'] > k else (k + 3 if p['ktuple'] else 3),
                percent=p['percent'])


def cw_distance(M, scale, avg, a, b, po, pe):
    """Exact-rational ClustalW distance 1 - identities / MIN(len) of the port (1 for an empty sequence)."""
    from fractions import Fraction
    if not a or not b: return Fraction(1)
    go, ge = gap_penalties(len(a), len(b), avg, scale, po, pe)
    cwp.MAT = M
    c = cwp.CW(a, b, go, ge)
    c.forward(); c.reverse(); c.printPtr = 1; c.lastPrint = 0
    c.diff(c.sb1 - 1, c.sb2 - 1, c.se1 - c.sb1 + 1, c.se2 - c.sb2 + 1, 0, 0)
    return 1 - Fraction(c.identities(), min(len(a), len(b)))


def msa_oracle(outdir, mats):
    spec2 = importlib.util.spec_from_file_location('mo', os.path.join(HERE, 'generate_msa_matrix_oracle.py'))
    mo = importlib.util.module_from_spec(spec2); spec2.loader.exec_module(mo)
    r2 = random.Random(20261126)
    rows = []
    pw_variants = [('GONNET', 10, 0.1), ('GONNET', 10, 0.1), ('BLOSUM', 10, 0.1), ('PAM', 5, 0.5), ('ID', 10, 0.1), ('GONNET', 20, 0.2)]
    for t in range(60):
        model = 'linear' if t % 2 == 0 else 'affine'
        scheme = mo.SCHEMES[t % len(mo.SCHEMES)]
        pm, po, pe = pw_variants[t % len(pw_variants)]
        M, scale, avg = mats[pm]
        mo.dist = lambda a, b, sc, model_, M=M, scale=scale, avg=avg, po=po, pe=pe: cw_distance(M, scale, avg, a, b, po, pe)
        mo.rnd = r2
        seqs = mo.make_set(r2, 36, mo.AA)
        sc = mo.Sc(*scheme)
        pr, it = mo.progressive_iterative(seqs, sc, model)
        co = mo.consistency(seqs, sc, model)
        rows.append('\t'.join([model, scheme[0], str(scheme[1]), str(scheme[2]), pm, str(po), str(pe), ';'.join(seqs), ';'.join(pr), ';'.join(it), ';'.join(co)]))
    with open(os.path.join(outdir, 'msa_matrix_clustalw_oracle.tsv'), 'w') as f:
        f.write('model\tmatrix\topen\text\tpwMatrix\tpwOpen\tpwExt\tseqs\tprogressive\titerative\tconsistency\n' + '\n'.join(rows) + '\n')
    print('msa:', len(rows), 'sets')


def matrix_reference(outdir, tri, eblosum30):
    from Bio.Align import substitution_matrices as sm
    order = AA_ORDER
    lines = []
    for name, key in (('BLOSUM30', 'blosum30mt'), ('GONNET250', 'gon250mt')):
        t = tri[key]; ix = 0; full = {}
        for i, a in enumerate(order):
            for j in range(i + 1):
                full[(a, order[j])] = t[ix]; full[(order[j], a)] = t[ix]; ix += 1
        if name == 'BLOSUM30':
            ref = {}
            rl = [l for l in open(eblosum30) if not l.startswith('#') and l.strip()]
            cols = rl[0].split()
            for l in rl[1:]:
                p = l.split()
                for c, v in zip(cols, p[1:]): ref[(p[0], c)] = int(v)
            keys = [(a, b) for a in order for b in order]
        else:
            g = sm.load('GONNET1992')
            ref = {(a, b): int(round(10 * g[a][b])) for a in AA20 for b in AA20}
            keys = [(a, b) for a in AA20 for b in AA20]
        for a, b in keys:
            assert full[(a, b)] == ref[(a, b)], (name, a, b, full[(a, b)], ref[(a, b)])
            lines.append(f'{name}\t{a}\t{b}\t{ref[(a, b)]}')
    with open(os.path.join(outdir, 'clustalw_protein_matrices_reference.tsv'), 'w') as f:
        f.write('matrix\ta\tb\tvalue\n' + '\n'.join(lines) + '\n')
    print('matrix reference cells:', len(lines))


def main(outdir, matrices_h, eblosum30=None):
    tri = parse_triangles(matrices_h)
    mats = {m: pairwise_matrix(tri, m) for m in NAMES}
    for m, (M, scale, avg) in mats.items(): print(m, 'intScale', scale, 'matrixAvgScore', avg)

    rows = []
    sets = 0
    gaps = [(None, None), (None, None), (10, 0.1), (5, 0.5), (20, 0.05), (0, 0.1), (10, 1.0), (7.5, 0.25), (15, 0.2), (2, 2)]
    plan = [('rand', 150, 70, AA20), ('mut', 130, 150, AA20), ('x', 40, 90, AA23), ('rep', 40, 80, "AL"), ('long', 6, 400, AA20)]
    for kind, n, maxlen, alpha in plan:
        for _ in range(n):
            matrix = rnd.choice(list(NAMES))
            po, pe = rnd.choice(gaps)
            seqs = make_set(maxlen, alpha)
            res = run_full(seqs, matrix, po, pe)
            M, scale, avg = mats[matrix]
            cwp.MAT = M
            sets += 1
            for (i, j), sc in sorted(res.items()):
                a, b = seqs[i], seqs[j]
                go, ge = gap_penalties(len(a), len(b), avg, scale, 10.0 if po is None else po, 0.1 if pe is None else pe)
                c = cwp.CW(a, b, go, ge)
                c.forward(); c.reverse(); c.printPtr = 1; c.lastPrint = 0
                c.diff(c.sb1 - 1, c.sb2 - 1, c.se1 - c.sb1 + 1, c.se2 - c.sb2 + 1, 0, 0)
                cnt = c.identities()
                mm = np.float32(np.float32(100.0) * np.float32(cnt)) / np.float32(min(len(a), len(b)))
                assert int(np.float32(mm)) == sc, (matrix, po, pe, a, b, sc, int(mm))
                rows.append((matrix, 10.0 if po is None else po, 0.1 if pe is None else pe, a, b, sc, cnt, c.maxScore, go, ge, avg))
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'clustalw_protein_pairwise_oracle.tsv'), 'w') as f:
        f.write('matrix\tgapOpen\tgapExtend\tseq1\tseq2\tclustalwScore\tidentities\tdpScore\tgapOpenInt\tgapExtendInt\tmatAvg\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
    print(f'full: {sets} sets, {len(rows)} pairs, all equal to the port')

    qrows = []
    qsets = 0
    plan = [('rand', 140, 70, AA20), ('mut', 100, 140, AA20), ('x', 30, 80, AA23), ('rep', 60, 70, "AL"), ('long', 5, 500, AA20)]
    for kind, n, maxlen, alpha in plan:
        for _ in range(n):
            p = quick_params()
            seqs = make_set(maxlen, alpha)
            if kind == 'rep':
                seqs = [gen(rnd.choice([2, 3, 1]), rnd.randint(20, maxlen), alpha) for _ in range(rnd.randint(2, 5))] + seqs[:2]
            maxlen_set = max(len(s) for s in seqs)
            res = run_quick(seqs, p)
            e = effective(p)
            qsets += 1
            for (i, j), sc in sorted(res.items()):
                qrows.append((kind, e['ktuple'], e['window'], e['topdiags'], e['pairgap'], int(e['percent']), maxlen_set, seqs[i], seqs[j], sc))
    with open(os.path.join(outdir, 'clustalw_protein_quick_oracle.tsv'), 'w') as f:
        f.write('kind\tktuple\twindow\ttopdiags\tpairgap\tpercent\tmaxLen\tseq1\tseq2\tprinted\n')
        for r in qrows: f.write('\t'.join(str(x) for x in r) + '\n')
    print(f'quick: {qsets} sets, {len(qrows)} pairs; not aligned: {sum(1 for r in qrows if r[-1] == "NA")}')
    msa_oracle(outdir, mats)
    if eblosum30: matrix_reference(outdir, tri, eblosum30)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
