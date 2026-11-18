#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-11) oracle: EXACT sum-of-pairs multiple alignment.

Seqeron API: SequenceAligner.MultipleAlignExact(sequences, scoring, MsaOptions{GapModel, ExactMaxStateCells, ExactCarrilloLipman}) and
MultipleAlignExactWithStatistics. Objective = SumOfPairsScore: linear = column sum of pairs (residue pair Match/Mismatch, residue vs gap
GapExtend, gap vs gap 0); affine = NATURAL affine sum of pairs (Altschul 1989: project every pair, drop gap-gap columns, every maximal
run of the same row's gap pays GapOpen + len*GapExtend).

Independent references (written from the definitions; nothing of the C# code is reused):
  * exact_linear(): DP over the k-dimensional lattice, 2^k - 1 neighbours (Altschul & Lipman 1989), memoised recursion.
  * exact_optimum() of generate_msa_affine_oracle.py: DP over (position vector, run state of EVERY pair) - the natural affine SP
    (3^pairs states, no weak-order reduction as in the C# program).
  * brute_force(): enumeration of EVERY column sequence of tiny sets, scored by sp_linear / sp_affine (asserted equal to both DPs).
  * k = 2: Biopython 1.88 PairwiseAligner(mode='global') (linear: open = extend = GapExtend; affine: open = GapOpen + GapExtend).
Output TSV (embedded: tests/.../TestData/Alignment/msa_exact_oracle.tsv), columns
  model(linear|affine) match mismatch open ext seqs(;-separated) optimum optimalRows(;-separated, one optimal alignment)
Usage: python3 scripts/alignment/generate_msa_exact_oracle.py > msa_exact_oracle.tsv
"""
import importlib.util, itertools, os, random, sys
from functools import lru_cache
from multiprocessing import Pool
from Bio.Align import PairwiseAligner

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('aff', os.path.join(HERE, 'generate_msa_affine_oracle.py'))
aff = importlib.util.module_from_spec(spec)
spec.loader.exec_module(aff)

sys.setrecursionlimit(1000000)
SCHEMES = [(1, -1, -2, -1), (2, -3, -5, -2), (5, -4, -10, -1), (1, -2, -3, -1), (3, -2, 0, -4), (2, -1, -6, -2), (1, -1, -4, 0), (1, -1, -1, -1),
           (4, -5, -3, -1), (0, 0, -2, -1), (1, 0, 0, -1)]


def exact_linear(seqs, sc):
    M, X, O, E = sc; k = len(seqs); L = tuple(len(s) for s in seqs)
    moves = [mv for mv in itertools.product((0, 1), repeat=k) if any(mv)]
    pairs = list(itertools.combinations(range(k), 2))

    @lru_cache(maxsize=None)
    def f(pos):
        if pos == L: return 0, ()
        best = None
        for mv in moves:
            npos = tuple(p + d for p, d in zip(pos, mv))
            if any(a > l for a, l in zip(npos, L)): continue
            col = tuple(seqs[r][pos[r]] if mv[r] else '-' for r in range(k))
            cost = 0
            for a, b in pairs:
                x, y = col[a], col[b]
                if x == '-' and y == '-': continue
                cost += E if (x == '-' or y == '-') else (M if x == y else X)
            sub, tail = f(npos)
            if best is None or cost + sub > best[0]: best = (cost + sub, (col,) + tail)
        return best
    v, cols = f(tuple([0] * k))
    rows = [''.join(c[r] for c in cols) for r in range(k)]
    assert aff.sp_linear(rows, sc) == v
    return v, rows


def brute_force(seqs, sc, model):
    k = len(seqs); L = tuple(len(s) for s in seqs)
    moves = [mv for mv in itertools.product((0, 1), repeat=k) if any(mv)]
    score = aff.sp_affine if model == 'affine' else aff.sp_linear
    best = [None]

    def rec(pos, cols):
        if pos == L:
            rows = [''.join(c[r] for c in cols) for r in range(k)]
            v = score(rows, sc)
            if best[0] is None or v > best[0]: best[0] = v
            return
        for mv in moves:
            npos = tuple(p + d for p, d in zip(pos, mv))
            if any(a > l for a, l in zip(npos, L)): continue
            rec(npos, cols + [tuple(seqs[r][pos[r]] if mv[r] else '-' for r in range(k))])
    rec(tuple([0] * k), [])
    return best[0]


def bio_pair(a, b, sc, model):
    M, X, O, E = sc
    if not a or not b: return (O * (1 if model == 'affine' else 0) + E * max(len(a), len(b))) if (a or b) else 0
    og = O + E if model == 'affine' else E
    return int(round(PairwiseAligner(mode='global', match_score=M, mismatch_score=X, open_gap_score=og, extend_gap_score=E).score(a, b)))


def solve(job):
    model, sc, seqs = job
    k = len(seqs)
    v, rows = exact_linear(seqs, sc) if model == 'linear' else aff.exact_optimum(seqs, sc)
    if k == 2: assert v == bio_pair(seqs[0], seqs[1], sc, model), (model, sc, seqs)
    if sum(map(len, seqs)) <= (10 if k >= 4 else 9): assert v == brute_force(seqs, sc, model), (model, sc, seqs, v)
    return '\t'.join(map(str, [model, *sc, ';'.join(seqs), v, ';'.join(rows)]))


def rseq(rnd, n, alpha): return ''.join(rnd.choice(alpha) for _ in range(n))


def mutate(rnd, s, rate):
    out = []
    for ch in s:
        r = rnd.random()
        if r < rate * 0.4: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.7: continue
        elif r < rate: out.append(ch); out.extend(rnd.choice('ACGT') for _ in range(rnd.randint(1, 3)))
        else: out.append(ch)
    return ''.join(out)


def make_jobs():
    rnd = random.Random(20261122)
    jobs = []
    # k = 3, lengths <= 12: related (mutated copies) and unrelated, both models
    for t in range(520):
        model = 'linear' if t % 2 == 0 else 'affine'
        sc = rnd.choice(SCHEMES)
        if rnd.random() < 0.6:
            base = rseq(rnd, rnd.randint(3, 10), rnd.choice(['ACGT', 'AC', 'AAAC']))
            seqs = [mutate(rnd, base, rnd.choice([0.1, 0.25, 0.4])) for _ in range(3)]
            seqs = [s[:12] for s in seqs]
        else:
            seqs = [rseq(rnd, rnd.randint(0 if t % 17 == 0 else 1, 12), rnd.choice(['ACGT', 'AC'])) for _ in range(3)]
        jobs.append((model, sc, seqs))
    # k = 2 (Biopython), k = 4 tiny (brute force), k = 5 linear tiny, k = 4 slightly longer
    for t in range(60):
        jobs.append(('linear' if t % 2 else 'affine', rnd.choice(SCHEMES), [rseq(rnd, rnd.randint(0, 14), 'ACGT') for _ in range(2)]))
    for t in range(80):
        model = 'linear' if t % 2 == 0 else 'affine'
        jobs.append((model, rnd.choice(SCHEMES), [rseq(rnd, rnd.randint(1, 3), rnd.choice(['ACGT', 'AC'])) for _ in range(4)]))
    for t in range(60):
        model = 'linear' if t % 2 == 0 else 'affine'
        base = rseq(rnd, rnd.randint(3, 6), 'ACGT')
        jobs.append((model, rnd.choice(SCHEMES), [mutate(rnd, base, 0.25)[:7] or 'A' for _ in range(4)]))
    for t in range(30):
        jobs.append(('linear', rnd.choice(SCHEMES), [rseq(rnd, rnd.randint(1, 4), 'ACGT') for _ in range(5)]))
    for t in range(14):
        jobs.append(('affine', rnd.choice(SCHEMES), [rseq(rnd, rnd.randint(1, 3), 'ACGT') for _ in range(5)]))
    # k = 3, longer related sets (length 18..26): pruning at a size where the lattice has ~10^4 cells
    for t in range(40):
        model = 'linear' if t % 2 == 0 else 'affine'
        base = rseq(rnd, rnd.randint(18, 24), 'ACGT')
        jobs.append((model, rnd.choice(SCHEMES), [mutate(rnd, base, rnd.choice([0.1, 0.25]))[:26] or 'A' for _ in range(3)]))
    return jobs


if __name__ == '__main__':
    jobs = make_jobs()
    with Pool(os.cpu_count()) as pool:
        rows = pool.map(solve, jobs, chunksize=4)
    print('model\tmatch\tmismatch\topen\text\tseqs\toptimum\toptimalRows')
    print('\n'.join(rows))
    print(f'{len(rows)} rows', file=sys.stderr)
