#!/usr/bin/env python3
"""B09 WP67 (F308, R14-13): an INDEPENDENT implementation of the conservation scores of Capra and Singh (2007) that produced the expected values of ConservationScores_Tests.cs / AlignInfoToolsTests.cs. It is written from the
definitions of the paper (frequency with the pseudo count 1e-7 over the 20 amino acids and the gap, h = sum f ln f / ln(min(21, rows)), JSD = (sum f log2(f / m) + sum b log2(b / m)) / 2 with m = (f + b) / 2 and the BLOSUM62 background
of the paper, the gap penalty 1 - weighted gap fraction), not from the program (score_conservation.py of compbio.cs.princeton.edu/conservation is a Python 2 CGI script and is not run here).
Usage: python3 scripts/alignment/conservation_scores_reference.py        (prints column, shannon, shannon without gap penalty, jsd, jsd without gap penalty)"""
import math

AA = 'ARNDCQEGHILKMFPSTWYV-'
BG = [0.078, 0.051, 0.041, 0.052, 0.024, 0.034, 0.059, 0.083, 0.025, 0.062, 0.092, 0.056, 0.024, 0.044, 0.043, 0.059, 0.055, 0.014, 0.034, 0.072]
PC = 1e-7


def freq(col, w=None):
    n = len(col)
    w = w or [1.0] * n
    c = [PC] * 21
    for k, a in enumerate(AA):
        for j in range(n):
            if col[j] == a:
                c[k] += w[j]
    tot = 0.0
    for x in w:
        tot += x
    d = tot + 21 * PC
    return [x / d for x in c]


def gap_penalty(col, w=None):
    n = len(col)
    w = w or [1.0] * n
    g = 0.0
    for i in range(n):
        if col[i] == '-':
            g += w[i]
    tot = 0.0
    for x in w:
        tot += x
    return 1 - g / tot


def shannon(col, penalty=True, w=None):
    f = freq(col, w)
    h = 0.0
    for x in f:
        if x != 0:
            h += x * math.log(x)
    h /= math.log(min(21, len(col)))
    s = 1 - (-1 * h)
    return s * gap_penalty(col, w) if penalty else s


def jsd(col, penalty=True, w=None, b=BG):
    f = freq(col, w)
    if len(b) == 20:
        f = f[:-1]
        s = 0.0
        for x in f:
            s += x
        f = [x / s for x in f]
    r = [.5 * x + .5 * y for x, y in zip(f, b)]
    d = 0.0
    for i in range(len(f)):
        if r[i] != 0:
            if f[i] == 0:
                d += b[i] * math.log(b[i] / r[i]) / math.log(2)
            elif b[i] == 0:
                d += f[i] * math.log(f[i] / r[i]) / math.log(2)
            else:
                d += f[i] * math.log(f[i] / r[i]) / math.log(2) + b[i] * math.log(b[i] / r[i]) / math.log(2)
    d /= 2
    return d * gap_penalty(col, w) if penalty else d


if __name__ == '__main__':
    for col in ['AAAAA', 'AAAAC', 'ACDEF', 'AAA--', 'ARNDC', 'LLLLLLLLLL', 'AALL']:
        print(col, repr(shannon(col)), repr(shannon(col, False)), repr(jsd(col)), repr(jsd(col, False)))
