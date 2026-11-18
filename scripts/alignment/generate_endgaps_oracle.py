#!/usr/bin/env python3
"""ALIGN-SEMI-001 (audit A1-3) oracle: ends-free (generalised semi-global) pairwise alignment, all 16 end-gap policies.

Seqeron API: SequenceAligner.SemiGlobalAlign / SemiGlobalAlignAffine(seq1, seq2, scoring, EndGaps) where EndGaps is a
flag set {FreeStart1=1, FreeEnd1=2, FreeStart2=4, FreeEnd2=8}; a free flag means the corresponding leading/trailing
overhang costs nothing (everything else is an ordinary gap: linear d per position, affine GapOpen + k*GapExtend).
Sequence 1 = Biopython target = parasail s1/query, sequence 2 = Biopython query = parasail s2/db:
  FreeStart1 = left_deletion_score = 0   = parasail _qb      FreeEnd1 = right_deletion_score = 0   = _qe
  FreeStart2 = left_insertion_score = 0  = parasail _db      FreeEnd2 = right_insertion_score = 0  = _de
Gap convention (MatrixScoring / ScoringMatrix): affine gap of length k = GapOpen + k*GapExtend; linear = GapExtend per
position. Biopython: open_gap_score = GapOpen + GapExtend, extend_gap_score = GapExtend (linear: both = GapExtend).
EMBOSS needle (default, i.e. -endweight N: end gaps free on BOTH sequences = EndGaps.All): -gapopen G -gapextend E with
G = -(GapOpen + GapExtend), E = -GapExtend.

Output TSV (embedded: tests/.../TestData/Alignment/endgaps_oracle.tsv), columns
  kind(scalar|matrix) model(linear|affine) flags(0..15) matrix match mismatch gapOpen gapExtend seq1 seq2
  biopython biopythonScoreApi parasail emboss alignments
biopython = PairwiseAligner.align(...).score (the optimum); biopythonScoreApi = PairwiseAligner.score(...) (quirk on degenerate
all-gap paths, see bio()); alignments = Biopython co-optimal gapped row pairs "row1|row2;..." (NA when more than 40 co-optimal alignments).
parasail: NA when the policy has no sg_* variant (3 flags) or the row has non-integer scores; the single-end variants
of parasail (sg, _qe, _de, _qb_de, _qe_db, _qe_de) never score ABOVE the optimum but are lower on rows whose optimum ends in the
border cell (m,0)/(0,n) (parasail does not consider those end cells) - the test documents this, Biopython and the exhaustive
reference agree with Seqeron on every row.

Usage: python3 scripts/alignment/generate_endgaps_oracle.py [n_pairs=70] > endgaps_oracle.tsv
"""
import random, re, subprocess, sys, tempfile
from Bio.Align import PairwiseAligner, substitution_matrices as sm
import parasail

rnd = random.Random(20261118)
PROT = 'ACDEFGHIKLMNPQRSTVWY'
SCHEMES = [(2, -3, 0, -2), (2, -3, -3, -1), (1, -1, -2, -1), (5, -4, -10, -1), (3, -2, 0, -4), (1, -2, -5, -2),
           (2, -3, -2, 0), (1, 0, -1, -1), (4, -5, 0, -6), (2, -1, -4, -3)]
MATRIX_GAPS = [(-9.5, -0.5), (-4.5, -0.5), (-2, -1), (0, -4), (-3.25, -0.25), (-7.5, -1.5), (-0.5, -0.5), (-11.5, -0.5)]
PARASAIL_NAME = {0: 'nw', 1: 'sg_qb', 2: 'sg_qe', 3: 'sg_qx', 4: 'sg_db', 5: 'sg_qb_db', 6: 'sg_qe_db', 8: 'sg_de',
                 9: 'sg_qb_de', 10: 'sg_qe_de', 12: 'sg_dx', 15: 'sg'}


def mutate(s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def gen_pair(alpha, maxlen):
    a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))
    b = mutate(a, alpha, rnd.choice([0.1, 0.3, 0.6]))
    r = rnd.random()
    if r < 0.2: b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))
    elif r < 0.45:  # dovetail: suffix of a == prefix of b
        k = rnd.randint(1, len(a)); b = a[len(a) - k:] + ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, maxlen // 2)))
    elif r < 0.6:  # containment
        i = rnd.randint(0, len(a)); j = rnd.randint(i, len(a)); b = a[i:j] or a[:1]
    return a, b


def configure(al, flags, model, go, ge):
    al.open_gap_score = ge if model == 'linear' else go + ge
    al.extend_gap_score = ge
    if flags & 1: al.left_deletion_score = 0
    if flags & 2: al.right_deletion_score = 0
    if flags & 4: al.left_insertion_score = 0
    if flags & 8: al.right_insertion_score = 0


def rows_of(al, target, query, alignment):
    """gapped (target, query) rows from the alignment coordinates"""
    c = alignment.coordinates
    r1, r2 = [], []
    for k in range(c.shape[1] - 1):
        dt = c[0][k + 1] - c[0][k]; dq = c[1][k + 1] - c[1][k]
        if dt and dq:
            assert dt == dq
            r1.append(target[c[0][k]:c[0][k + 1]]); r2.append(query[c[1][k]:c[1][k + 1]])
        elif dt:
            r1.append(target[c[0][k]:c[0][k + 1]]); r2.append('-' * dt)
        elif dq:
            r1.append('-' * dq); r2.append(query[c[1][k]:c[1][k + 1]])
    return ''.join(r1), ''.join(r2)


def bio(a, b, flags, model, go, ge, matrix, ma, mi, with_rows):
    al = PairwiseAligner()
    al.mode = 'global'
    if matrix: al.substitution_matrix = sm.load(matrix)
    else: al.match_score = ma; al.mismatch_score = mi
    configure(al, flags, model, go, ge)
    # The alignment score of align() is the optimum (= exhaustive search on every row); PairwiseAligner.score() has a quirk
    # on the degenerate "one sequence entirely against gaps" paths with free ends (Biopython 1.88: score('AC', 'T') with
    # left_deletion_score = right_insertion_score = 0 is -2.0 while the alignment 'AC-'/'--T' it returns scores 0.0), so
    # the oracle is align().score; score() is kept in its own column to document the quirk.
    alns = al.align(a, b)
    score = float(alns.score)
    score_api = float(al.score(a, b))
    rows = 'NA'
    if with_rows:
        try:
            if len(alns) <= 40:
                rows = ';'.join(sorted({'|'.join(rows_of(al, a, b, x)) for x in alns}))
        except OverflowError:
            pass
    return score, score_api, rows


def parasail_score(a, b, flags, model, go, ge, ma, mi):
    if flags not in PARASAIL_NAME or go != int(go) or ge != int(ge): return None
    mat = parasail.matrix_create('ACGT', ma, mi)
    opn = -ge if model == 'linear' else -(go + ge)
    return float(getattr(parasail, PARASAIL_NAME[flags])(a, b, int(opn), int(-ge), mat).score)


def emboss_needle(a, b, go, ge, model, matrix, kind, d):
    G, E = (-ge, -ge) if model == 'linear' else (-(go + ge), -ge)
    for v in (G, E):
        if abs(round(v, 1) - v) > 1e-9: return None
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind == 'dna' else ['-sprotein1', '-sprotein2']
    cmd = ['needle', '-asequence', d + '/a.fa', '-bsequence', d + '/b.fa', '-gapopen', repr(round(G, 1)), '-gapextend',
           repr(round(E, 1)), '-datafile', {'BLOSUM62': 'EBLOSUM62', 'EDNAFULL': 'EDNAFULL'}[matrix], '-awidth', '1000',
           '-outfile', d + '/o.txt', '-auto'] + typ  # default -endweight N: free end gaps on both sequences
    if subprocess.run(cmd, capture_output=True, text=True).returncode != 0: return None
    m = re.search(r'# Score: ([-0-9.]+)', open(d + '/o.txt').read())
    return float(m.group(1)) if m else None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 70
    print('kind\tmodel\tflags\tmatrix\tmatch\tmismatch\tgapOpen\tgapExtend\tseq1\tseq2\tbiopython\tbiopythonScoreApi\tparasail\temboss\talignments')
    bad = 0
    with tempfile.TemporaryDirectory() as d:
        for idx in range(n):
            a, b = gen_pair('ACGT' if rnd.random() < 0.7 else 'AC', 24)
            for model in ('linear', 'affine'):
                for flags in range(16):
                    ma, mi, go, ge = rnd.choice(SCHEMES)
                    bs, bsapi, rows = bio(a, b, flags, model, go, ge, None, ma, mi, with_rows=len(a) <= 22 and len(b) <= 22)
                    ps = parasail_score(a, b, flags, model, go, ge, ma, mi)
                    if ps is not None and ps > bs + 1e-9:
                        bad += 1; print(f'# parasail ABOVE biopython {flags} {model} {a} {b}', file=sys.stderr)
                    print('\t'.join(['scalar', model, str(flags), '-', str(ma), str(mi), str(go), str(ge), a, b, repr(bs), repr(bsapi),
                                     'NA' if ps is None else repr(ps), 'NA', rows]))
        # protein / IUPAC matrix rows (fractional gap scores) incl. EMBOSS needle default (end gaps free on both sequences)
        for idx in range(n // 2):
            kind = rnd.choice(['prot', 'dna'])
            alpha, mat = (PROT, 'BLOSUM62') if kind == 'prot' else ('ACGT', 'EDNAFULL')
            a, b = gen_pair(alpha, 45)
            for model in ('linear', 'affine'):
                for flags in range(16):
                    go, ge = rnd.choice(MATRIX_GAPS)
                    bs, bsapi, _ = bio(a, b, flags, model, go, ge, 'BLOSUM62' if mat == 'BLOSUM62' else 'NUC.4.4', 0, 0, False)
                    es = emboss_needle(a, b, go, ge, model, mat, kind, d) if flags == 15 and len(a) > 1 and len(b) > 1 else None
                    print('\t'.join(['matrix', model, str(flags), mat, '0', '0', repr(float(go)), repr(float(ge)), a, b, repr(bs), repr(bsapi), 'NA',
                                     'NA' if es is None else repr(es), 'NA']))
        # more EMBOSS needle (default -endweight N) rows: protein / nucleotide, flags = All
        for idx in range(n * 2):
            kind = rnd.choice(['prot', 'prot', 'dna'])
            alpha, mat = (PROT, 'BLOSUM62') if kind == 'prot' else ('ACGT', 'EDNAFULL')
            a, b = gen_pair(alpha, 60)
            if len(a) < 2 or len(b) < 2: continue
            for model in ('linear', 'affine'):
                go, ge = rnd.choice(MATRIX_GAPS)
                es = emboss_needle(a, b, go, ge, model, mat, kind, d)
                if es is None: continue
                bs, bsapi, _ = bio(a, b, 15, model, go, ge, 'BLOSUM62' if mat == 'BLOSUM62' else 'NUC.4.4', 0, 0, False)
                print('\t'.join(['matrix', model, '15', mat, '0', '0', repr(float(go)), repr(float(ge)), a, b, repr(bs), repr(bsapi), 'NA',
                                 repr(es), 'NA']))
    print(f'# parasail-above violations: {bad}', file=sys.stderr)


main()
