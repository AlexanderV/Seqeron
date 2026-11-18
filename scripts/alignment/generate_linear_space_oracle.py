#!/usr/bin/env python3
"""ALIGN-GLOBAL/LOCAL/SEMI-001 (audit A1-6) oracle: exact linear-space pairwise alignment
(SequenceAligner.{Global,Local,SemiGlobal}Align[Affine]LinearSpace and the automatic switch above 25 M cells).

Three independent references (the TSV is embedded: tests/.../TestData/Alignment/linear_space_oracle.tsv):

  bio rows       Biopython 1.88 PairwiseAligner (mode global / local, left/right_*_score = 0 for the ends-free policies) and
                 parasail 1.3.4 (nw / sw / sg_dx / sg_qb_de ... where a variant exists). `biopython` = align().score (the optimum);
                 `alignments` = every co-optimal gapped row pair (NA when more than 40). Pairs: random, mutated, planted long
                 indels (5-40 nt), tie-rich (alphabets A and AC), lengths up to 120, scalar schemes incl. GapOpen 0 / GapExtend 0.
  stretcher rows EMBOSS 6.6.0 stretcher = Myers & Miller (1988) linear-space global alignment, end gaps penalised. EMBOSS gap
                 convention: a gap of length k costs gapopen + (k-1)*gapextend, i.e. Seqeron GapOpen = -(G - E), GapExtend = -E
                 with MatrixScoring (EDNAFULL / EBLOSUM62). Columns: score + the two aligned rows.
  matcher rows   EMBOSS 6.6.0 matcher = Huang & Miller (1991) SIM, linear-space local alignment (first alternative: the optimal local
                 alignment). Same gap convention. Columns: score, aligned rows, 1-based start/end of both sequences.

Usage: python3 scripts/alignment/generate_linear_space_oracle.py [n_pairs=160] > linear_space_oracle.tsv
Output columns: kind family model flags matrix match mismatch gapOpen gapExtend seq1 seq2 score parasail extra
  kind = bio | stretcher | matcher; family = global | local | semi; flags = EndGaps value (bio semi rows) or 0;
  extra = bio: co-optimal rows "r1|r2;..." (NA > 40); stretcher: "row1|row2|u"; matcher: "row1|row2|s1|e1|s2|e2|u" where u = 1 when
          Biopython 1.88 (same matrix, open = -G, extend = -E) finds exactly one optimal alignment, else 2 (the EMBOSS row pair must then
          equal Seqeron's; co-optimal alignments may be tie-broken differently).
"""
import random, re, subprocess, sys, tempfile
from Bio.Align import PairwiseAligner, substitution_matrices as sm
import parasail

rnd = random.Random(20261118 + 6)
PROT = 'ACDEFGHIKLMNPQRSTVWY'
SCHEMES = [(2, -3, 0, -2), (2, -3, -3, -1), (1, -1, -2, -1), (5, -4, -10, -1), (3, -2, 0, -4), (1, -2, -5, -2),
           (2, -3, -2, 0), (1, 0, -1, -1), (4, -5, 0, -6), (2, -1, -4, -3), (1, -1, 0, 0)]
EMBOSS_GAPS = [(16, 4), (10, 0.5), (12, 2), (8, 8), (20, 1)]
PARASAIL_FLAGS = {0: 'nw', 12: 'sg_dx', 4: 'sg_db', 1: 'sg_qb', 15: 'sg'}  # parasail variants that agree with the optimum


def mutate(s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def plant_long_gaps(s, alpha):
    """insert / delete one or two blocks of 5-40 nt"""
    s = list(s)
    for _ in range(rnd.randint(1, 2)):
        pos = rnd.randint(0, len(s))
        if rnd.random() < 0.5:
            s[pos:pos] = [rnd.choice(alpha) for _ in range(rnd.randint(5, 40))]
        else:
            k = rnd.randint(5, 40)
            del s[pos:pos + k]
    return ''.join(s) or rnd.choice(alpha)


def gen_pair(alpha, maxlen):
    a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))
    r = rnd.random()
    if r < 0.15: b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, maxlen)))
    elif r < 0.55: b = plant_long_gaps(mutate(a, alpha, rnd.choice([0.05, 0.15])), alpha)
    elif r < 0.75:  # containment / overlap shapes
        i = rnd.randint(0, len(a)); j = rnd.randint(i, len(a)); b = (a[i:j] or a[:1])
        if rnd.random() < 0.5: b = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 12))) + b
    else: b = mutate(a, alpha, rnd.choice([0.1, 0.3, 0.6]))
    return a, b


def configure(al, flags, model, go, ge):
    al.open_gap_score = ge if model == 'linear' else go + ge
    al.extend_gap_score = ge
    if flags & 1: al.left_deletion_score = 0
    if flags & 2: al.right_deletion_score = 0
    if flags & 4: al.left_insertion_score = 0
    if flags & 8: al.right_insertion_score = 0


def rows_of(target, query, alignment):
    c = alignment.coordinates
    r1, r2 = [], []
    for k in range(c.shape[1] - 1):
        dt = c[0][k + 1] - c[0][k]; dq = c[1][k + 1] - c[1][k]
        if dt and dq:
            r1.append(target[c[0][k]:c[0][k + 1]]); r2.append(query[c[1][k]:c[1][k + 1]])
        elif dt:
            r1.append(target[c[0][k]:c[0][k + 1]]); r2.append('-' * dt)
        elif dq:
            r1.append('-' * dq); r2.append(query[c[1][k]:c[1][k + 1]])
    return ''.join(r1), ''.join(r2)


def bio(a, b, mode, flags, model, go, ge, ma, mi, with_rows):
    al = PairwiseAligner()
    al.mode = mode
    al.match_score = ma; al.mismatch_score = mi
    configure(al, flags, model, go, ge)
    alns = al.align(a, b)
    score = float(alns.score)
    rows = 'NA'
    if with_rows:
        try:
            if len(alns) <= 40:
                rows = ';'.join(sorted({'|'.join(rows_of(a, b, x)) for x in alns}))
        except OverflowError:
            pass
    return score, rows


def parasail_score(a, b, mode, flags, model, go, ge, ma, mi):
    name = 'sw' if mode == 'local' else PARASAIL_FLAGS.get(flags)
    if name is None: return None
    mat = parasail.matrix_create('ACGT', ma, mi)
    opn = -ge if model == 'linear' else -(go + ge)
    s = float(getattr(parasail, name)(a, b, int(opn), int(-ge), mat).score)
    return s


def bio_unique(a, b, mode, matrix, G, E):
    """1 when Biopython finds exactly ONE optimal alignment (same matrix, open = -G, extend = -E), else 2"""
    al = PairwiseAligner()
    al.mode = mode
    al.substitution_matrix = sm.load('BLOSUM62' if matrix == 'EBLOSUM62' else 'NUC.4.4')
    al.open_gap_score = -float(G)
    al.extend_gap_score = -float(E)
    alns = al.align(a, b)
    try:
        return 1 if len(alns) == 1 else 2
    except OverflowError:
        return 2


def emboss_rows(text):
    rows = {}
    for line in text.splitlines():
        m = re.match(r'^([ab])\s+(\d+) ([A-Za-z-]+)\s+(\d+)\s*$', line)
        if m: rows.setdefault(m.group(1), []).append((int(m.group(2)), m.group(3), int(m.group(4))))
    return rows


def emboss(prog, a, b, G, E, matrix, kind, d):
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind == 'dna' else ['-sprotein1', '-sprotein2']
    cmd = [prog, d + '/a.fa', d + '/b.fa', '-gapopen', repr(float(G)) if G != int(G) else str(int(G)),
           '-gapextend', repr(float(E)) if E != int(E) else str(int(E)),
           '-datafile', matrix, '-aformat3', 'srspair', '-awidth', '100000', '-outfile', d + '/o.txt', '-auto'] + typ
    if prog == 'matcher': cmd += ['-alternatives', '1']
    if subprocess.run(cmd, capture_output=True, text=True).returncode != 0: return None
    text = open(d + '/o.txt').read()
    ms = re.search(r'# Score: ([-0-9.]+)', text)
    if not ms: return None
    rows = emboss_rows(text)
    if 'a' not in rows or 'b' not in rows or len(rows['a']) != 1: return None
    (s1, r1, e1), (s2, r2, e2) = rows['a'][0], rows['b'][0]
    return float(ms.group(1)), r1.upper(), r2.upper(), s1, e1, s2, e2


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    print('kind\tfamily\tmodel\tflags\tmatrix\tmatch\tmismatch\tgapOpen\tgapExtend\tseq1\tseq2\tscore\tparasail\textra')
    with tempfile.TemporaryDirectory() as d:
        # ---- Biopython / parasail ----
        for idx in range(n):
            alpha = rnd.choice(['ACGT', 'ACGT', 'AC', 'A'])
            a, b = gen_pair(alpha, 120 if alpha == 'ACGT' else 40)
            for model in ('linear', 'affine'):
                for fam, mode, flags in (('global', 'global', 0), ('local', 'local', 0), ('semi', 'global', 12),
                                         ('semi', 'global', rnd.choice([1, 2, 3, 4, 5, 6, 9, 10, 15]))):
                    ma, mi, go, ge = rnd.choice(SCHEMES)
                    with_rows = len(a) <= 40 and len(b) <= 40
                    bs, rows = bio(a, b, mode, flags, model, go, ge, ma, mi, with_rows)
                    ps = parasail_score(a, b, mode, flags, model, go, ge, ma, mi)
                    if ps is not None and ps > bs + 1e-9: print(f'# parasail ABOVE biopython {fam} {flags} {a} {b}', file=sys.stderr)
                    print('\t'.join(['bio', fam, model, str(flags), '-', str(ma), str(mi), str(go), str(ge), a, b, repr(bs),
                                     'NA' if ps is None else repr(ps), rows]))
        # ---- EMBOSS stretcher (Myers-Miller) and matcher (Huang-Miller SIM) ----
        done = {'stretcher': 0, 'matcher': 0}
        while done['stretcher'] < 160 or done['matcher'] < 160:
            kind = rnd.choice(['prot', 'dna', 'dna'])
            alpha, mat = (PROT, 'EBLOSUM62') if kind == 'prot' else ('ACGT', 'EDNAFULL')
            a, b = gen_pair(alpha, rnd.choice([40, 120, 300]))
            if len(a) < 2 or len(b) < 2: continue
            G, E = rnd.choice(EMBOSS_GAPS)
            for prog in ('stretcher', 'matcher'):
                if done[prog] >= 160: continue
                r = emboss(prog, a, b, G, E, mat, kind, d)
                if r is None: continue
                score, r1, r2, s1, e1, s2, e2 = r
                done[prog] += 1
                uniq = bio_unique(a, b, 'global' if prog == 'stretcher' else 'local', mat, G, E)
                extra = f'{r1}|{r2}|{uniq}' if prog == 'stretcher' else f'{r1}|{r2}|{s1}|{e1}|{s2}|{e2}|{uniq}'
                print('\t'.join([prog, 'global' if prog == 'stretcher' else 'local', 'affine', '0', ('BLOSUM62' if mat == 'EBLOSUM62' else 'EDNAFULL'),
                                 '0', '0', repr(-(float(G) - float(E))), repr(-float(E)), a, b, repr(score), 'NA', extra]))


main()
