#!/usr/bin/env python3
"""GENOMIC-TANDEM-001 oracle: exact tandem repeats = maximal repetitions (runs) with a primitive unit.

Definition (Kolpakov & Kucherov 1999 "Finding maximal repetitions in a word in linear time"; Gusfield 1997 sec. 9
"tandem arrays"; TRF, Benson 1999: "two or more contiguous copies of a pattern"; microsatellite = a *primitive* unit
repeated, MISA "reject false type motifs"): for unit length p, a repeat is a maximal interval S[a..e) with period p
(S[x] = S[x+p] for all a <= x < e-p) that cannot be extended left or right, has at least minReps complete copies
(floor((e-a)/p) >= minReps) and whose unit S[a..a+p) is primitive (not a power of a shorter word). It is reported once,
at its left end a, with floor((e-a)/p) copies. Output order: unit length ascending, then start ascending.

Computed THREE independent ways per (sequence, minUnit, minReps, maxUnit) row, asserted equal here:
  A. Kolpakov-Kucherov runs by brute force: every interval [a,e), its *smallest* period by definition, maximality by
     trying to extend; the interval is a run iff e-a >= 2*smallest period; every unit length p >= smallest period that
     is a *multiple*-free report is derived from the run: only p == smallest period has a primitive unit (a unit of
     length p = m*q of a run with smallest period q is a power of the length-q unit => not primitive).
  B. Boolean-array formulation: eq_p[x] = (S[x]==S[x+p]); maximal True-runs of length L give the interval
     [a, a+L+p); complete copies = (L+p)//p; primitive test via (u+u).find(u,1) == len(u).
  C. Regex: for unit u in all words of length p occurring in S: re.finditer of (?<!...) is awkward, so use the
     lookahead formulation  (?=((.{p})\\2{k-1,}))  restricted to left-maximal starts (S[a-1] != S[a-1+p] or a==0).
Usage: python3 generate_tandem_oracle.py > tandem_oracle.tsv
Columns: sequence minUnit minReps maxUnit(0 = unbounded) count hits("unit:pos:copies;..." | "-")

Ambiguous mode (A1-21, GenomicAnalyzer.FindTandemRepeats(string ...)): python3 generate_tandem_oracle.py --ambiguous > tandem_iupac_oracle.tsv
Same three independent methods, sequences over alphabets with N and other IUPAC symbols; the definition adds the microsatellite
criterion of MISA (`[acgt]`) / pytrf (skips N) / RepeatFinder.FindMicrosatellites: only a unit spelled with A, C, G, T can repeat, so
a window containing any other symbol never forms a repeat (the N-containing candidates of every method are dropped).
"""
import random
import re
import sys

AMBIGUOUS = "--ambiguous" in sys.argv
rnd = random.Random(20261125 if AMBIGUOUS else 20261118)


def primitive(u):
    # AMBIGUOUS: a unit must also be made of A/C/G/T only (MISA [acgt])
    return (u + u).find(u, 1) == len(u) and (not AMBIGUOUS or set(u) <= set("ACGT"))


def method_a(s, min_u, min_r, max_u):
    n = len(s)
    out = set()
    for a in range(n):
        for e in range(a + 2, n + 1):
            t = s[a:e]
            m = len(t)
            q = next(q for q in range(1, m + 1) if all(t[x] == t[x + q] for x in range(m - q)))
            if 2 * q > m:
                continue
            # maximal for period q: cannot extend left / right while keeping period q
            if a > 0 and s[a - 1] == s[a - 1 + q]:
                continue
            if e < n and s[e] == s[e - q]:
                continue
            copies = m // q
            if AMBIGUOUS and not set(s[a:a + q]) <= set("ACGT"):
                continue
            if q >= min_u and (max_u == 0 or q <= max_u) and copies >= min_r:
                out.add((q, a, s[a:a + q], copies))
    return sorted((p, a, u, c) for p, a, u, c in out)


def method_b(s, min_u, min_r, max_u):
    n = len(s)
    out = []
    top = n // min_r if max_u == 0 else min(max_u, n // min_r)
    for p in range(min_u, top + 1):
        x = 0
        while x + p < n:
            if s[x] != s[x + p]:
                x += 1
                continue
            y = x
            while y + p < n and s[y] == s[y + p]:
                y += 1
            length = (y - x) + p
            copies = length // p
            u = s[x:x + p]
            if copies >= min_r and primitive(u):
                out.append((p, x, u, copies))
            x = y
    return sorted(out)


def method_c(s, min_u, min_r, max_u):
    n = len(s)
    out = []
    top = n // min_r if max_u == 0 else min(max_u, n // min_r)
    for p in range(min_u, top + 1):
        pat = re.compile(r'(?=((.{%d})\2{%d,}))' % (p, min_r - 1), re.S)
        for m in pat.finditer(s):
            a = m.start()
            if a > 0 and s[a - 1] == s[a - 1 + p]:
                continue
            e = a + p
            while e < n and s[e] == s[e - p]:
                e += 1
            u = s[a:a + p]
            if primitive(u):
                out.append((p, a, u, (e - a) // p))
    return sorted(out)


def fmt(hits):
    return ";".join(f"{u}:{a}:{c}" for p, a, u, c in hits) or "-"


def row(s, min_u, min_r, max_u):
    a = method_a(s, min_u, min_r, max_u)
    b = method_b(s, min_u, min_r, max_u)
    c = method_c(s, min_u, min_r, max_u)
    assert a == b == c, (s, min_u, min_r, max_u, a, b, c)
    return f"{s}\t{min_u}\t{min_r}\t{max_u}\t{len(a)}\t{fmt(a)}"


def main():
    if AMBIGUOUS:
        main_ambiguous()
        return
    crafted = [
        "ATATATAT", "AAAAA", "AAAAAA", "ATATATA", "TATATAT", "ATGATGATG", "ATTCGATTCGATTCG", "CCATGATGATGCC",
        "AAATTT", "AAAAAATTTTT", "CAGCAGCAGCAGCAG", "TTAGGGTTAGGGTTAGGGTTAGGG", "GATAGATAGATA", "ACGTACACGTACACGTAC",
        "ACACACACGTGTGTGT", "AAAAAAAAAAAA", "ACGACGACGACGAC", "GCATATATGC", "TTCAGCAGCAGTT", "ACGT", "AC", "A",
        "ACACACACACACAC", "AATAATAATAAT", "GAGAGAGTGAGAGAGA", "CACACACATATATATA", "AAGAAGAAGAAGAAG", "TGTGTGTGAAA",
        "ACGTACGTACGTAAAA", "CGCGCGCGCGAA", "ATATATATATATATATATAT", "AACAACAACAACAACAAC",
    ]
    rows = []
    for s in crafted:
        for min_u, min_r, max_u in [(1, 2, 0), (2, 2, 0), (1, 3, 0), (2, 3, 0), (3, 2, 0), (1, 2, 2), (2, 2, 4), (4, 2, 0)]:
            rows.append(row(s, min_u, min_r, max_u))
    for _ in range(900):
        alpha = rnd.choice(["AC", "ACG", "ACGT", "AT", "ACGT"])
        n = rnd.randint(0, 36)
        s = "".join(rnd.choice(alpha) for _ in range(n))
        # plant a repeat sometimes
        if n >= 8 and rnd.random() < 0.6:
            u = "".join(rnd.choice(alpha) for _ in range(rnd.randint(1, 4)))
            k = rnd.randint(2, 5)
            pos = rnd.randint(0, n - 1)
            s = s[:pos] + u * k + s[pos:]
        min_u = rnd.randint(1, 3)
        min_r = rnd.randint(2, 4)
        max_u = rnd.choice([0, 0, 1, 2, 3, 5])
        if max_u and max_u < min_u:
            max_u = min_u
        rows.append(row(s, min_u, min_r, max_u))
    sys.stdout.write("\n".join(rows) + "\n")


def main_ambiguous():
    crafted = [
        "ACACNACAC", "NNNNNNNN", "ACNACNACN", "ACGNACGNACGN", "NACACACN", "ANANANAN", "NANANANA", "AAAANAAAA", "ATATATRTATAT",
        "ACGTNNNNACGTACGT", "CAGCAGCAGNNNCAGCAGCAG", "RYRYRYRY", "ACACRCACAC", "GAGAGAGANGAGAGAGA", "TTTTTNTTTTT", "NNNACACACNNN",
        "ACGTACGTYACGTACGT", "AAAAAAAAAAAAN", "NAAAAAAAAAAAA", "ACNGTACNGT", "CCNCCNCCN", "ATGATGATGNATGATGATG", "KKKKKKK", "ACACACACAC",
    ]
    rows = []
    for s in crafted:
        for min_u, min_r, max_u in [(1, 2, 0), (2, 2, 0), (1, 3, 0), (2, 3, 0), (3, 2, 0), (1, 2, 2), (2, 2, 4), (4, 2, 0)]:
            rows.append(row(s, min_u, min_r, max_u))
    for _ in range(700):
        alpha = rnd.choice(["ACN", "ACGTN", "ACGTN", "ATN", "ACGTNRY", "ACGN", "AN", "ACGTACGTACGTN"])
        n = rnd.randint(0, 40)
        s = "".join(rnd.choice(alpha) for _ in range(n))
        if n >= 8 and rnd.random() < 0.7:
            u = "".join(rnd.choice(alpha) for _ in range(rnd.randint(1, 4)))
            k = rnd.randint(2, 5)
            pos = rnd.randint(0, n - 1)
            s = s[:pos] + u * k + s[pos:]
        min_u = rnd.randint(1, 3)
        min_r = rnd.randint(2, 4)
        max_u = rnd.choice([0, 0, 1, 2, 3, 5])
        if max_u and max_u < min_u:
            max_u = min_u
        rows.append(row(s, min_u, min_r, max_u))
    sys.stdout.write("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
