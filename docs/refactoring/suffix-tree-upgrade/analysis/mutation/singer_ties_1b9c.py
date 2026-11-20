#!/usr/bin/env python3
"""Work package 1b.9c, the argument for the mutant D10 (DifferenceCover.NearestSingerPeriod keeps the first modulus on a tie of the ratios): "<" and "<=" differ only for a target that has two Singer moduli (q^2 + q + 1, q prime, q < 200: 46 values from 7 to 39,801)
   at exactly the same ratio max(x / t, t / x) in double arithmetic. This runs both rules (the same IEEE operations as the C# code) for every target t = 1 .. 300,000 and prints the targets where they differ (none). A tie of a modulus below t and one above t needs t^2 = x1 * x2 (up to rounding), so t <= 39,801; above that
   the nearest modulus is the largest one under either rule: the range covers every target that can tie. usage: singer_ties_1b9c.py"""
def prime(v):
    if v < 2:
        return False
    d = 2
    while d * d <= v:
        if v % d == 0:
            return False
        d += 1
    return True


xs = [q * q + q + 1 for q in range(2, 200) if prime(q)]


def nearest(t, strict):
    best, best_ratio = 7, float('inf')
    for x in xs:
        r = max(x / t, t / x)
        if (r < best_ratio) if strict else (r <= best_ratio):
            best_ratio, best = r, x
    return best


diff = [t for t in range(1, 300001) if nearest(t, True) != nearest(t, False)]
print(len(xs), 'Singer moduli, targets 1 .. 300000 where "<" and "<=" differ:', diff)
