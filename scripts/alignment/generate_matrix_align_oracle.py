#!/usr/bin/env python3
"""ALIGN-GLOBAL/LOCAL/SEMI (audit A1-2) oracle: pairwise alignment scores with a substitution matrix and
real-valued gap scores - EMBOSS needle/water 6.6.0, Biopython 1.88 PairwiseAligner(substitution_matrix=...),
parasail 1.3.4 (integer gaps only).

Seqeron convention (MatrixScoring): GapOpen/GapExtend are scores (<= 0); affine gap of length k = GapOpen + k*GapExtend;
linear = GapExtend per position. Biopython: open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend.
EMBOSS: needle -endweight Y -endopen G -endextend E (needle's DEFAULT leaves end gaps free) with -gapopen G -gapextend E, G = -(GapOpen + GapExtend), E = -GapExtend (G = cost of a gap of length 1);
linear model = needle/water with gapopen == gapextend.
parasail: open = -(GapOpen + GapExtend), extend = -GapExtend (nw / sw / sg_dx).
semi-global = fitting (seq1 inside seq2, free seq2 flanks) = Biopython end_insertion_score = 0 = parasail sg_dx
(not available in EMBOSS).

Output TSV (embedded in tests/.../TestData/Alignment/matrix_align_oracle.tsv); columns
 mode(global|local|semi) model(linear|affine) matrix gapOpen gapExtend seq1 seq2 biopython emboss parasail
(NA when the tool cannot express the row). Rows are asserted equal across the available tools.

Usage: python3 scripts/alignment/generate_matrix_align_oracle.py [n_pairs=420] > matrix_align_oracle.tsv
"""
import random, subprocess, re, sys, tempfile
from Bio.Align import PairwiseAligner, substitution_matrices as sm
import parasail

rnd = random.Random(20261118)
PROT = 'ACDEFGHIKLMNPQRSTVWY'
IUPAC = 'ACGTRYSWKMBDHVN'
EMBOSS_MATRIX = {'BLOSUM62': 'EBLOSUM62', 'PAM250': 'EPAM250', 'EDNAFULL': 'EDNAFULL'}
BIO_MATRIX = {'BLOSUM62': 'BLOSUM62', 'PAM250': 'PAM250', 'EDNAFULL': 'NUC.4.4'}
PARASAIL = {'BLOSUM62': parasail.blosum62, 'PAM250': parasail.pam250, 'EDNAFULL': parasail.nuc44}
HBA = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR'
HBB = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'

# (GapOpen, GapExtend) in Seqeron convention: integers, EMBOSS defaults, fractional
GAPS = [(-2, -1), (-9, -1), (-10, -1), (-5, -2), (0, -4), (-9.5, -0.5), (-9.9, -0.1), (-11.5, -0.5), (-4.5, -0.5),
        (-3.25, -0.25), (-7.5, -1.5), (-0.5, -0.5), (-2, -0.1), (-12, -0.3), (-3, -2.5), (-6.75, -0.75)]


def mutate(s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def gen():
    kind = rnd.choice(['prot', 'prot', 'protx', 'dna', 'iupac'])
    if kind == 'dna': alpha, mat = 'ACGT', 'EDNAFULL'
    elif kind == 'iupac': alpha, mat = IUPAC, 'EDNAFULL'
    else:
        alpha = PROT + ('BZX' if kind == 'protx' else '')
        mat = rnd.choice(['BLOSUM62', 'PAM250'])
    a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 70)))
    b = mutate(a, alpha, rnd.choice([0.1, 0.3, 0.6]))
    if rnd.random() < 0.3: b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 50)))
    if rnd.random() < 0.35: b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 25))) + b + ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 25)))
    return a, b, mat, kind


def bio_score(mode, model, mat, go, ge, a, b):
    al = PairwiseAligner()
    al.substitution_matrix = sm.load(BIO_MATRIX[mat])
    al.mode = 'local' if mode == 'local' else 'global'
    ext = ge
    opn = ge if model == 'linear' else go + ge
    al.open_gap_score = opn
    al.extend_gap_score = ext
    if mode == 'semi':
        al.end_insertion_score = 0
    return float(al.score(a, b))


def emboss_score(prog, mat, go, ge, model, a, b, kind, d):
    # EMBOSS numbers are positive: G = cost of a gap of length 1
    if model == 'linear': G, E = -ge, -ge
    else: G, E = -(go + ge), -ge
    for v in (G, E):
        if abs(round(v, 1) - v) > 1e-9: return None  # needle/water print one decimal
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind in ('dna', 'iupac') else ['-sprotein1', '-sprotein2']
    cmd = [prog, '-asequence', d + '/a.fa', '-bsequence', d + '/b.fa', '-gapopen', repr(round(G, 1)), '-gapextend', repr(round(E, 1)),
           '-datafile', EMBOSS_MATRIX[mat], '-awidth', '1000', '-outfile', d + '/o.txt', '-auto'] + typ
    if prog == 'needle':
        # needle's default (-endweight N) makes end gaps FREE on both sequences; the unit's global alignment charges them
        cmd += ['-endweight', 'Y', '-endopen', repr(round(G, 1)), '-endextend', repr(round(E, 1))]
    if subprocess.run(cmd, capture_output=True, text=True).returncode != 0: return None
    m = re.search(r'# Score: ([-0-9.]+)', open(d + '/o.txt').read())
    return float(m.group(1)) if m else None


def parasail_score(mode, model, mat, go, ge, a, b):
    if go != int(go) or ge != int(ge): return None
    opn = -(0 if model == 'linear' else go) - ge
    ext = -ge
    # linear: open = extend = -ge  (first gap position already costs ge)
    if model == 'linear': opn = ext
    f = {'global': parasail.nw, 'local': parasail.sw, 'semi': parasail.sg_dx}[mode]
    return float(f(a, b, int(opn), int(ext), PARASAIL[mat]).score)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 420
    pairs = [(HBA, HBB, 'BLOSUM62', 'prot')] + [gen() for _ in range(n)]
    print('mode\tmodel\tmatrix\tgapOpen\tgapExtend\tseq1\tseq2\tbiopython\temboss\tparasail')
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        for idx, (a, b, mat, kind) in enumerate(pairs):
            for mode in ('global', 'local', 'semi'):
                for model in ('linear', 'affine'):
                    go, ge = rnd.choice(GAPS)
                    if idx == 0: go, ge = -9.5, -0.5  # needle defaults gapopen 10 gapextend 0.5
                    if a == '' or b == '': continue
                    bs = bio_score(mode, model, mat, go, ge, a, b)
                    es = None
                    if mode != 'semi' and len(a) > 1 and len(b) > 1:  # water finds no 1-residue local alignment
                        es = emboss_score('needle' if mode == 'global' else 'water', mat, go, ge, model, a, b, kind, d)
                    ps = parasail_score(mode, model, mat, go, ge, a, b)
                    for name, v in (('emboss', es), ('parasail', ps)):
                        if v is not None and abs(v - bs) > 0.051:
                            bad += 1
                            print(f'# MISMATCH {name} {mode} {model} {mat} {go} {ge} {a} {b} bio={bs} {name}={v}', file=sys.stderr)
                    print('\t'.join([mode, model, mat, repr(float(go)), repr(float(ge)), a, b, repr(bs),
                                     'NA' if es is None else repr(es), 'NA' if ps is None else repr(ps)]))
    print(f'# oracle disagreements: {bad}', file=sys.stderr)


main()
