#!/usr/bin/env python3
"""GENOMIC-SIMILARITY-001 oracle (A2-31): MinHash / FracMinHash / containment / ANI, sourmash 4.9.4 vs independent Python sets.

Output TSV (tests/.../TestData/Analysis/similarity_sketch_oracle.tsv):
  seqA  seqB  k  kind  param  jaccard  contAinB  contBinA  maxContainment  aniAinB  aniBinA  jaccardAni
kind:
  exact   (param 0)   exact canonical k-mer sets (Python sets of min(w, revcomp w), windows with a non-ACGT symbol skipped);
                      containment C(A,B) = |A&B|/|A|; ANI = sourmash containment_to_distance(...).ani (= C^(1/k), 0 for C = 0, 1 for C = 1);
                      jaccardAni = 1 - sourmash jaccard_to_distance(J, k, 1).dist = (2J/(1+J))^(1/k)   (asserted equal to the formula)
  scaled  (param S)   sourmash MinHash(n=0, ksize=k, scaled=S, seed=42): jaccard, contained_by both ways, max_containment (sourmash's own values)
  num     (param s)   sourmash MinHash(n=s, ksize=k, seed=42) bottom-s sketch: jaccard (== the `mash dist` estimator for k >= 17; Mash truncates the
                      hash to 32 bits for k <= 16, sourmash does not, so num rows use k >= 17); the other columns are '.'
For `exact`, sourmash MinHash(scaled=1) is asserted to give the same Jaccard / containments (it holds every k-mer).
Usage: python3 scripts/analysis/generate_similarity_sketch_oracle.py > similarity_sketch_oracle.tsv   (never aborts: ERR rows)
"""
import random, sys
from sourmash import MinHash
from sourmash.distance_utils import containment_to_distance, jaccard_to_distance

rnd = random.Random(20261122)
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s))


def kmers(s, k):
    out = set()
    for i in range(len(s) - k + 1):
        w = s[i:i + k]
        if all(c in COMP for c in w):
            out.add(min(w, rc(w)))
    return out


def rand(n):
    return ''.join(rnd.choice('ACGT') for _ in range(n))


def mutate(s, rate):
    return ''.join(rnd.choice('ACGT') if rnd.random() < rate else c for c in s)


def pair():
    n = rnd.randint(30, 160)
    a = rand(n)
    kind = rnd.random()
    if kind < 0.45:
        b = mutate(a, rnd.choice([0.01, 0.03, 0.08, 0.2]))
    elif kind < 0.65:
        b = rc(mutate(a, rnd.choice([0.0, 0.05])))
    elif kind < 0.85:
        cut = rnd.randint(1, n - 1)
        b = a[:cut] + rand(rnd.randint(5, 80))
    else:
        b = rand(rnd.randint(30, 150))
    if rnd.random() < 0.15:
        i = rnd.randrange(len(a)); a = a[:i] + 'N' + a[i + 1:]
    return a, b


def sm(seq, k, n=0, scaled=0):
    mh = MinHash(n=n, ksize=k, scaled=scaled, seed=42)
    mh.add_sequence(seq, force=True)
    return mh


def f(x):
    return repr(float(x))


def main():
    rows = 0
    for _ in range(330):
        a, b = pair()
        for k in rnd.sample([3, 5, 7, 9, 11, 15, 21, 31], 3):
            A, B = kmers(a, k), kmers(b, k)
            if not A and not B:
                continue
            inter, union = len(A & B), len(A | B)
            J = inter / union if union else 0.0
            cab = inter / len(A) if A else 0.0
            cba = inter / len(B) if B else 0.0
            mx = max(cab, cba)
            try:
                ma, mb = sm(a, k, scaled=1), sm(b, k, scaled=1)
                assert abs(ma.jaccard(mb) - J) < 1e-12 and abs(ma.contained_by(mb) - cab) < 1e-12 and abs(mb.contained_by(ma) - cba) < 1e-12, (a, b, k)
                def ani(c, size):
                    return containment_to_distance(c, k, 1, n_unique_kmers=max(size, 1)).ani
                ja = 1.0 - jaccard_to_distance(J, k, 1, n_unique_kmers=max(len(A), 1)).dist
                assert abs(ja - (0.0 if J == 0 else 1.0 if J == 1 else (2 * J / (1 + J)) ** (1.0 / k))) < 1e-12
                print('\t'.join([a, b, str(k), 'exact', '0', f(J), f(cab), f(cba), f(mx), f(ani(cab, len(A))), f(ani(cba, len(B))), f(ja)]))
                rows += 1
                for S in rnd.sample([2, 3, 5, 10], 2):
                    sa, sb = sm(a, k, scaled=S), sm(b, k, scaled=S)
                    c1, c2 = sa.contained_by(sb), sb.contained_by(sa)
                    j = sa.jaccard(sb)
                    mc = sa.max_containment(sb)
                    ja = 1.0 - jaccard_to_distance(j, k, S, n_unique_kmers=max(len(A), 1)).dist
                    print('\t'.join([a, b, str(k), 'scaled', str(S), f(j), f(c1), f(c2), f(mc), f(ani(c1, len(A))), f(ani(c2, len(B))), f(ja)]))
                    rows += 1
            except Exception as e:
                sys.stderr.write('ERR %s %s\n' % (k, e))
                print('\t'.join([a, b, str(k), 'ERR', '0', '.', '.', '.', '.', '.', '.', '.']))
        for k in rnd.sample([17, 19, 21, 25, 31], 2):
            for s in rnd.sample([5, 10, 30, 100, 500], 2):
                try:
                    j = sm(a, k, n=s).jaccard(sm(b, k, n=s))
                except Exception as e:
                    sys.stderr.write('ERR %s\n' % e); continue
                print('\t'.join([a, b, str(k), 'num', str(s), f(j), '.', '.', '.', '.', '.', '.']))
                rows += 1
    sys.stderr.write('rows=%d\n' % rows)


if __name__ == '__main__':
    main()
