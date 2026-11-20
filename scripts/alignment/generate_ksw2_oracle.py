#!/usr/bin/env python3
"""ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit A2-29b, B09 F125): ksw2 extension alignment (ksw_extz, ksw_extz2_sse, ksw_extd, ksw_extd2_sse of lh3/ksw2 and minimap2).

Two outputs from one deterministic case generator (python3 generate_ksw2_oracle.py MODE OUTFILE [seed]):
  cases   : every case (exact ones and the heuristic ones: Z-drop, approximate max / drop) as the input of the C harness of docs/Validation/review-2026-09/handoff/A2-29b-ksw2/ (the real ksw2 functions are not
            runnable in the review sandbox: they have to be compiled from the published source, which the sandbox refuses). Line:  id variant flag w zdrop end_bonus q e q2 e2 m mat... query target
  oracle  : the EXACT cases (no Z-drop, no approximate flags) with the expected values from an independent Gotoh dynamic program over the band (Python, written from the definition: H over the cells |i - j| <= w,
            gap of length k costs q + k e or min(q + k e, q2 + k e2)): max, mqe, mte and score. The program does not look at the Seqeron or the ksw2 code.
            TSV (tests/.../TestData/Alignment/ksw2_exact_oracle.tsv): id variant flag w end_bonus q e q2 e2 m mat query target max mqe mte score
Sequences are digit strings (A = 0, C = 1, G = 2, T = 3, N = 4), "-" = empty.
"""
import random, sys

NEG = -(1 << 30)
SCORE_ONLY, RIGHT, GENERIC, APPROX_MAX, APPROX_DROP, EXTZ_ONLY, REV_CIGAR = 1, 2, 4, 8, 16, 0x40, 0x80
VARIANTS = ['extz', 'extz2', 'extd', 'extd2']


def simple_mat(m, a, b, ambi):
    mat = [0] * (m * m)
    for i in range(m - 1):
        for j in range(m - 1):
            mat[i * m + j] = a if i == j else -b
        mat[i * m + m - 1] = -ambi
    for j in range(m):
        mat[(m - 1) * m + j] = -ambi
    return mat


def mutate(rnd, s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate: out.append(rnd.randrange(4))
        elif r < rate * 1.3: continue
        elif r < rate * 1.6:
            out.append(c); out.extend(rnd.randrange(4) for _ in range(rnd.randint(1, 8)))
        else: out.append(c)
    return out


def make_case(rnd, variant, exact):
    simd = variant.endswith('2')
    dual = variant.startswith('extd')
    a = rnd.choice([1, 2, 2, 3, 5])
    b = rnd.randint(1, 6)
    ambi = rnd.choice([0, 1, 2, 3])
    q = rnd.randint(0, 8)
    e = rnd.randint(2, 4) if dual else rnd.randint(1, 4)
    if simd:
        b = min(b, 2 * (q + e))
        ambi = min(ambi, 2 * (q + e))       # ksw2 returns at once when a score is below -2 (q + e): the dynamic range of the 8-bit differences
    if dual:
        e2 = rnd.randint(1, e - 1)
        q2 = q + (e - e2) + rnd.randint(0, 30)      # q + e <= q2 + e2: the short-gap piece is the cheaper one for k = 1 (ksw_extd2_sse swaps the pieces otherwise)
         # the second piece must have the smaller extension cost (minimap2 -E 2,1): with e2 == e the long-gap boundary values of ksw_extd2_sse are not defined
    else:
        q2 = e2 = 0
    flag = 0
    if rnd.random() < 0.15: flag |= SCORE_ONLY
    if rnd.random() < 0.4: flag |= RIGHT
    if rnd.random() < 0.3: flag |= EXTZ_ONLY
    if rnd.random() < 0.3: flag |= REV_CIGAR
    m = 5
    mat = simple_mat(m, a, b, ambi)
    if rnd.random() < 0.3:
        flag |= GENERIC
        for i in range(4):
            for j in range(4):
                if (i, j) != (0, 0):
                    mat[i * m + j] = rnd.randint(max(-b, -2 * (q + e)), a - 1) if i != j else rnd.randint(max(1, a - 2), a)
    n = rnd.choice([1, 2, 3, 5, 8, 12, 20, 33, 50, 70, 120])
    base = [rnd.randrange(4) if rnd.random() > 0.04 else 4 for _ in range(n)]
    t = mutate(rnd, base, rnd.choice([0.0, 0.05, 0.15, 0.3]))
    if rnd.random() < 0.15: t = [rnd.randrange(4) for _ in range(rnd.choice([1, 4, 17, 60]))]
    if rnd.random() < 0.3: t = t[rnd.randint(0, len(t) // 3):]
    if rnd.random() < 0.3: t = t[:max(1, len(t) - rnd.randint(0, len(t) // 3))]
    qs = base
    if rnd.random() < 0.3: qs = base[:max(1, len(base) - rnd.randint(0, len(base) // 2))]
    if not t: t = [0]
    # the SIMD functions compute whole 16-lane blocks, so their band is blurred to multiples of 16 at the edges: a banded SIMD case is only decided by the real program (heuristic set)
    banded_ok = (not simd) or (not exact)
    w = -1 if (rnd.random() < 0.45 or not banded_ok) else rnd.choice([1, 2, 3, 5, 8, 15, 30] if simd else [0, 1, 2, 3, 5, 8, 15, 30])
    if w >= 0 and not simd and not (flag & (EXTZ_ONLY | SCORE_ONLY)):
        w = max(w, abs(len(qs) - len(t)))       # the scalar functions read the backtrack matrix unchecked: the corner must be inside the band
    end_bonus = rnd.choice([0, 0, 3, 10, 40]) if simd else 0
    zdrop = -1
    if not exact:                                    # heuristic cases: Z-drop and the approximate flags of the SIMD functions
        r = rnd.random()
        if r < 0.55:
            zdrop = rnd.choice([0, 2, 5, 10, 20, 50, 100])
        if simd and r > 0.35:
            flag |= APPROX_MAX | (APPROX_DROP if rnd.random() < 0.6 else 0)
        if zdrop < 0 and not (flag & (APPROX_MAX | APPROX_DROP)):
            zdrop = rnd.choice([5, 30, 100])
    # the SIMD functions need a + 2 (q + e) to fit a signed byte and the matrix to fit the 8-bit range
    qe = q + e
    if simd and (mat[0] + 2 * qe > 120): q = 3; e = 1
    return dict(variant=variant, flag=flag, w=w, zdrop=zdrop, end_bonus=end_bonus, q=q, e=e, q2=q2, e2=e2, m=m, mat=mat, query=qs, target=t)


def gotoh(c):
    """Best scores over the band: (max, mqe, mte, score) by the definition."""
    q, e, q2, e2 = c['q'], c['e'], c['q2'], c['e2']
    dual = c['variant'].startswith('extd')
    simd = c['variant'].endswith('2')
    qs, ts, m, mat, flag = c['query'], c['target'], c['m'], c['mat'], c['flag']
    qlen, tlen = len(qs), len(ts)
    w = c['w'] if c['w'] >= 0 else max(qlen, tlen)

    def g(k):
        return min(q + k * e, q2 + k * e2) if dual else q + k * e

    def sc(i, j):
        s, t = ts[i], qs[j]
        if simd and not (flag & GENERIC):                       # match / mismatch / wildcard of the SIMD functions
            if s == m - 1 or t == m - 1:
                nn = mat[m * m - 1]
                return -(e2 if dual else e) if nn == 0 else nn
            return mat[0] if s == t else mat[1]
        return mat[s * m + t]

    H = [[NEG] * (qlen + 1) for _ in range(tlen + 1)]
    H[0][0] = 0
    mx = 0
    for i in range(tlen + 1):
        for j in range(qlen + 1):
            if (i == 0 and j == 0) or abs(i - j) > w:
                continue
            best = NEG
            if i > 0 and j > 0 and H[i - 1][j - 1] > NEG:
                best = H[i - 1][j - 1] + sc(i - 1, j - 1)
            for k in range(1, j + 1):
                if H[i][j - k] > NEG: best = max(best, H[i][j - k] - g(k))
            for k in range(1, i + 1):
                if H[i - k][j] > NEG: best = max(best, H[i - k][j] - g(k))
            H[i][j] = best
            if i > 0 and j > 0 and best > mx: mx = best
    mqe = max((H[i][qlen] for i in range(1, tlen + 1) if H[i][qlen] > NEG), default=NEG)
    mte = max((H[tlen][j] for j in range(1, qlen + 1) if H[tlen][j] > NEG), default=NEG)
    return mx, mqe, mte, H[tlen][qlen]


def seq(s): return ''.join(map(str, s)) if s else '-'


def main(mode, outfile, seed):
    rnd = random.Random(seed)
    rows = []
    n = 0
    for variant in VARIANTS:
        for k in range(330 if mode == 'oracle' else 330 + 170):
            exact = k < 330
            c = make_case(rnd, variant, exact)
            if mode == 'oracle' and not exact: continue
            mx, mqe, mte, score = gotoh(c) if exact else (0, 0, 0, 0)
            head = [str(n), variant]
            if mode == 'cases':
                rows.append(' '.join(map(str, [n, variant, c['flag'], c['w'], c['zdrop'], c['end_bonus'], c['q'], c['e'], c['q2'], c['e2'], c['m']] + c['mat'] + [seq(c['query']), seq(c['target'])])))
            else:
                rows.append('\t'.join(map(str, [n, variant, c['flag'], c['w'], c['end_bonus'], c['q'], c['e'], c['q2'], c['e2'], c['m'], ','.join(map(str, c['mat'])), seq(c['query']), seq(c['target']),
                                                 mx, mqe if mqe > NEG else 'NEG', mte if mte > NEG else 'NEG', score if score > NEG else 'NEG'])))
            n += 1
    with open(outfile, 'w') as f:
        if mode == 'oracle': f.write('id\tvariant\tflag\tw\tend_bonus\tq\te\tq2\te2\tm\tmat\tquery\ttarget\tmax\tmqe\tmte\tscore\n')
        f.write('\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 20261126)
