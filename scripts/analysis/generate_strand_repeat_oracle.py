#!/usr/bin/env python3
"""GENOMIC-REPEAT-001 oracle (A2-33 / A2-35): both-strand repeats and repeat pairs, MUMmer 3.23 `repeat-match`.

Two independent references per (sequence, L):
  1. `repeat-match -n L` (default = BOTH strands) and `repeat-match -f -n L` (forward only), parsed to 0-based pairs
     (i, j, length, strand): a reverse-complement line `Start1 Start2r Length` means S[i..i+L) == revcomp(S[k..k+L)) with
     i = Start1-1 and k = Start2 - Length (Start2 is the 1-based position of the LAST base of the second copy);
  2. a brute force of the both-strand repeat classes over T = S # revcomp(S) (non-ACGT symbols are unique
     separators: they never match and always end a repeat): for every canonical substring alpha = min(a, revcomp a) of
     length >= L, the forward starts F (S[i..] == alpha) and reverse starts R (S[i..] == revcomp(alpha)), the class is a
     repeat when |F| + |R| >= 2 (a palindrome counts both strands) and a MAXIMAL repeat when, over the occurrences in T,
     the left contexts and the right contexts both take >= 2 distinct values (every out-of-range / non-ACGT neighbour is
     its own value; an R occurrence reads the complement of the opposite neighbour).
The script reports (stderr) every ACGT sequence where the class set of the repeat-match pairs differs from the brute-force maximal classes
(expected: none). repeat-match is run on S + 'N' (B04 F25: without the sentinel S and revcomp(S) share the final terminator, the tool drops palindromic pairs such as (0,0,2) of ATAAT and aborts on exact-palindrome genomes).

TSV (tests/.../TestData/Analysis/strand_repeat_oracle.tsv), one line per (sequence, L):
  sequence  L  maximal(csv  alpha:F:R ...)  all(alpha:F:R ...)  pairs(i:j:len:s;... s=f|r, from repeat-match)  forward_pairs(from -f)
  'all' is '?' when not computed (only sequences <= 26 nt, L in {2,4}); with F / R = comma lists of 0-based forward-strand window starts ('-' when empty); pairs are '-' for sequences with N.
Usage: python3 scripts/analysis/generate_strand_repeat_oracle.py > strand_repeat_oracle.tsv   (never aborts: a failed tool call gives ERR)
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261120)
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}


def rc(s):
    return ''.join(COMP.get(c, c) for c in reversed(s))


def rand(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


def classes(S, L, maximal_only):
    n = len(S)
    cand = set()
    for l in range(L, n + 1):
        for i in range(n - l + 1):
            w = S[i:i + l]
            if all(c in COMP for c in w):
                cand.add(min(w, rc(w)))
    out = {}
    for a in cand:
        l = len(a)
        b = rc(a)
        F = [i for i in range(n - l + 1) if S[i:i + l] == a]
        R = [i for i in range(n - l + 1) if S[i:i + l] == b]
        if len(F) + len(R) < 2:
            continue
        if maximal_only:
            left, right = set(), set()
            for i in F:  # occurrence of a read on the forward strand
                left.add(S[i - 1] if i > 0 and S[i - 1] in COMP else ('u', 'L', i))
                right.add(S[i + l] if i + l < n and S[i + l] in COMP else ('u', 'R', i))
            for i in R:  # occurrence of a on the reverse strand: neighbours are complements of the opposite ones
                left.add(COMP[S[i + l]] if i + l < n and S[i + l] in COMP else ('v', 'L', i))
                right.add(COMP[S[i - 1]] if i > 0 and S[i - 1] in COMP else ('v', 'R', i))
            if len(left) < 2 or len(right) < 2:
                continue
        out[a] = (F, R)
    return out


def fmt(cl):
    if not cl:
        return '-'
    return ';'.join('%s:%s:%s' % (a, ','.join(map(str, F)) or '-', ','.join(map(str, R)) or '-') for a, (F, R) in sorted(cl.items()))


def repeat_match(f, L, forward_only):
    args = ['repeat-match'] + (['-f'] if forward_only else []) + ['-n', str(L), f]
    r = subprocess.run(args, capture_output=True, text=True, timeout=60)
    if r.returncode != 0:  # repeat-match aborts on a genome that is an exact palindrome ('*** Genome is exact palindrome ***')
        return None
    pairs = []
    for ln in r.stdout.splitlines()[2:]:
        p = ln.split()
        if len(p) != 3:
            continue
        i, length = int(p[0]) - 1, int(p[2])
        if p[1].endswith('r'):
            pairs.append((i, int(p[1][:-1]) - length, length, 'r'))
        else:
            pairs.append((i, int(p[1]) - 1, length, 'f'))
    return sorted(pairs)


def planted(alpha, with_n):
    parts = [rand(rnd.randint(0, 6), alpha)]
    for _ in range(rnd.randint(1, 3)):
        u = rand(rnd.randint(3, 8), alpha)
        parts.append(u)
        parts.append(rand(rnd.randint(0, 5), alpha))
        parts.append(rc(u) if rnd.random() < 0.55 else u)
        parts.append(rand(rnd.randint(0, 5), alpha))
    s = ''.join(parts)
    if with_n and len(s) > 3:
        for _ in range(rnd.randint(1, 2)):
            i = rnd.randrange(len(s))
            s = s[:i] + 'N' * rnd.randint(1, 3) + s[i + 1:]
    return s


def main():
    seqs = []
    for alpha in ('ACGT', 'AC', 'AT', 'ACGT', 'ACGT'):
        for _ in range(75):
            seqs.append(planted(alpha, False) if rnd.random() < 0.75 else rand(rnd.randint(4, 28), alpha))
    for _ in range(90):
        seqs.append(planted('ACGT', True))
    seqs += ['GAATTC', 'AAAAAAA', 'ACGTACGT', 'TTTTAAAA', 'GGATCCGGATCC', 'ACGGTTCCGT', 'AAACCCGGGTTT', 'ACGNACGT', 'NNNN', 'GAATTCNGAATTC']
    seqs = [s for s in seqs if s]
    mism = 0
    rows = 0
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, 's.fa')
        for S in seqs:
            open(f, 'w').write('>s\n' + S + 'N\n')  # N sentinel: see B04 F25 (S and revcomp(S) would share the final $ and repeat-match loses pairs / aborts on palindromic genomes)
            has_n = any(c not in COMP for c in S)
            for L in (1, 2, 3, 4):
                try:
                    mx = classes(S, L, True)
                    al = classes(S, L, False) if (len(S) <= 26 and L in (2, 4)) else None
                    if has_n or L < 2:
                        pairs = fpairs = None
                    else:
                        pairs, fpairs = repeat_match(f, L, False), repeat_match(f, L, True)
                except Exception as e:
                    sys.stderr.write('ERR %s %s %s\n' % (S, L, e))
                    print('%s\t%d\tERR\tERR\tERR\tERR' % (S, L))
                    continue
                if pairs is not None and fpairs is not None:
                    cls = {min(S[i:i + ln], rc(S[i:i + ln])) for (i, j, ln, s) in pairs}
                    if cls != set(mx):
                        mism += 1
                        sys.stderr.write('CLASS MISMATCH %s L=%d rm=%s brute=%s\n' % (S, L, sorted(cls), sorted(mx)))
                ps = '-' if pairs is None else (';'.join('%d:%d:%d:%s' % p for p in pairs) or '-')
                fs = '-' if fpairs is None else (';'.join('%d:%d:%d:%s' % p for p in fpairs) or '-')
                print('%s\t%d\t%s\t%s\t%s\t%s' % (S, L, fmt(mx), '?' if al is None else fmt(al), ps, fs))
                rows += 1
    sys.stderr.write('rows=%d class-mismatches=%d\n' % (rows, mism))


if __name__ == '__main__':
    main()
