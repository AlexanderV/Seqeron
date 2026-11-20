#!/usr/bin/env python3
"""GENOMIC-SIMILARITY-001 oracle (A3-11): sourmash 4.9.4 ANI estimates with confidence intervals.

For random related sequence pairs (400-4000 nt, point-mutated copies at 0.5-25 %, plus unrelated pairs, optional N), k in {5..31}, FracMinHash scaled S in {1..20}, confidence in
{0.90, 0.95, 0.99}: sourmash `MinHash(n=0, ksize=k, scaled=S, seed=42)` of both sequences and, per pair, the four ANI results
  c_ab : A.containment_ani(B, estimate_ci=True, confidence=c)
  c_ba : B.containment_ani(A, estimate_ci=True, confidence=c)
  max  : A.max_containment_ani(B, estimate_ci=True, confidence=c)
  jac  : A.jaccard_ani(B)
(sourmash.distance_utils; n_unique_kmers = len(sketch) * scaled exactly as MinHash.containment_ani / jaccard_ani do), each as a row.
TSV (tests/.../TestData/Analysis/ani_ci_oracle.tsv):
  seqA seqB k scaled confidence kind  dist  dist_low  dist_high  p_nothing_in_common  p_exceeds_threshold(0|1)  size_is_inaccurate(0|1)  jaccard_error  je_exceeds_threshold(0|1)
(`.` = None: no confidence interval could be estimated / not applicable). Repr of the floats is exact.
Never aborts: a failing call is reported on stderr and the row skipped.
Usage: python3 scripts/analysis/generate_ani_ci_oracle.py [pairs] > ani_ci_oracle.tsv
"""
import random, sys, warnings
from sourmash import MinHash

warnings.filterwarnings('ignore')
rnd = random.Random(20261124)
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s))


def rand(n):
    return ''.join(rnd.choice('ACGT') for _ in range(n))


def mutate(s, rate):
    return ''.join(rnd.choice('ACGT') if rnd.random() < rate else c for c in s)


def sm(seq, k, scaled):
    mh = MinHash(n=0, ksize=k, scaled=scaled, seed=42)
    mh.add_sequence(seq, force=True)
    return mh


def f(x):
    return '.' if x is None else repr(float(x))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 260
    rows = 0
    import io, contextlib
    for _ in range(n):
        L = rnd.choice([400, 800, 1500, 4000])
        a = rand(L)
        t = rnd.random()
        if t < 0.8:
            b = mutate(a, rnd.choice([0.005, 0.01, 0.03, 0.06, 0.1, 0.15, 0.25]))
        elif t < 0.9:
            b = rc(mutate(a, 0.02))
        else:
            b = rand(rnd.choice([L, L // 2]))
        if rnd.random() < 0.1:
            i = rnd.randrange(len(a)); a = a[:i] + 'N' + a[i + 1:]
        k = rnd.choice([5, 7, 9, 11, 13, 15, 17, 21, 25, 31])
        S = rnd.choice([1, 1, 2, 3, 5, 10, 20])
        conf = rnd.choice([0.90, 0.95, 0.95, 0.99])
        A, B = sm(a, k, S), sm(b, k, S)
        if len(A) == 0 or len(B) == 0:
            continue
        try:
            with contextlib.redirect_stderr(io.StringIO()):
                res = {
                    'c_ab': A.containment_ani(B, estimate_ci=True, confidence=conf),
                    'c_ba': B.containment_ani(A, estimate_ci=True, confidence=conf),
                    'max': A.max_containment_ani(B, estimate_ci=True, confidence=conf),
                    'jac': A.jaccard_ani(B),
                }
        except Exception as e:
            sys.stderr.write('FAIL %s\n' % e)
            continue
        for kind, r in res.items():
            if kind == 'jac':
                je = r.jaccard_error
                print('\t'.join([a, b, str(k), str(S), repr(conf), kind, f(r.dist), '.', '.', f(r.p_nothing_in_common), str(int(r.p_exceeds_threshold)),
                                 str(int(r.size_is_inaccurate)), f(je), str(int(r.je_exceeds_threshold))]))
            else:
                print('\t'.join([a, b, str(k), str(S), repr(conf), kind, f(r.dist), f(r.dist_low), f(r.dist_high), f(r.p_nothing_in_common),
                                 str(int(r.p_exceeds_threshold)), str(int(r.size_is_inaccurate)), '.', '.']))
            rows += 1
    sys.stderr.write('rows=%d\n' % rows)


if __name__ == '__main__':
    main()
