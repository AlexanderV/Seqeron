#!/usr/bin/env python3
"""Exact emulation of R's mean() (summary.c real_mean: long-double sum / n, then a second pass adding sum(x - s) / n,
rounded once to double) with rational arithmetic and a 64-bit mantissa (x87 extended, round-to-nearest-even, exponent range
ignored). Shows that R's algorithm is NOT always within one ulp of the correctly rounded mean, which is why
StatisticsHelperProperties.ExtendedPrecisionMean_IntegerInput_WithinOneUlpOfExact uses the bound
  1 ulp + 2^-62 * mean(|x - exact|).

Usage:  python3 -I scripts/verification/r_mean_longdouble_emulation.py [samples]
Prints the 22-integer example (CI counterexample, FsCheck 3.6.0, 2026-11-04) and a Monte-Carlo estimate over random
integer vectors (n 1..40, |x| <= 1e6, the same distribution as the property).
"""
import math, random, sys
from fractions import Fraction


def rnd(q, p=64):
    if q == 0:
        return Fraction(0)
    sgn = 1 if q > 0 else -1
    a = abs(q)
    e = a.numerator.bit_length() - a.denominator.bit_length() - p
    while True:
        m = a / (Fraction(2) ** e)
        if m >= 2 ** p:
            e += 1
        elif m < 2 ** (p - 1):
            e -= 1
        else:
            break
    fl = m.numerator // m.denominator
    r = m - fl
    if r > Fraction(1, 2) or (r == Fraction(1, 2) and fl % 2 == 1):
        fl += 1
    return sgn * Fraction(fl) * Fraction(2) ** e


def r_mean(xs):
    n = len(xs)
    s = Fraction(0)
    for x in xs:
        s = rnd(s + x)
    s = rnd(s / n)
    t = Fraction(0)
    for x in xs:
        t = rnd(t + rnd(Fraction(x) - s))
    return float(rnd(s + rnd(t / n)))


example = [939853, 812265, -544392, 558617, -82511, -513122, 30411, 527327, -409764, -715409, -11039, 518869, 830861,
           -695734, 221617, -818309, -108322, 987640, -43404, 471213, -970525, -985298]
exact = sum(example) / len(example)
m = r_mean(example)
u = math.ulp(exact)
print(f"example: sum={sum(example)} exact={exact!r} R-algorithm={m!r} deviation={abs(m - exact) / u:.1f} ulp")

samples = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
random.seed(11)
over = 0
worst = 0.0
worst_ratio = 0.0
for _ in range(samples):
    xs = [random.randint(-1_000_000, 1_000_000) for _ in range(random.randint(1, 40))]
    ex = sum(xs) / len(xs)
    mm = r_mean(xs)
    uu = math.ulp(ex) if ex != 0 else 5e-324
    d = abs(mm - ex) / uu
    over += d > 1
    worst = max(worst, d)
    mad = sum(abs(x - ex) for x in xs) / len(xs)
    worst_ratio = max(worst_ratio, abs(mm - ex) / (uu + mad * 2.0 ** -62))
print(f"{samples} random vectors: more than 1 ulp in {over} ({100 * over / samples:.2f} %), worst {worst:.1f} ulp, "
      f"max |err| / (1 ulp + 2^-62 * meanAbsDev) = {worst_ratio:.3f}")
