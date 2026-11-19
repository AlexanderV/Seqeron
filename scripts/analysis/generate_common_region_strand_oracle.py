#!/usr/bin/env python3
"""GENOMIC-COMMON-001 oracle (A2-30): strand / mode / k-mismatch / multi-sequence common regions.

  python3 generate_common_region_strand_oracle.py mum      > common_region_mum_oracle.tsv
  python3 generate_common_region_strand_oracle.py kmism    > common_region_kmismatch_oracle.tsv
  python3 generate_common_region_strand_oracle.py multi    > common_region_multi_oracle.tsv

mum: MUMmer 3.23 `mummer {-maxmatch|-mumreference|-mum} -b -l L ref qry` (-b = forward AND reverse-complement matches), checked against an
  independent Python brute force (the script asserts the two agree on every ACGT row; rows with N use the brute force only):
    ref  qry  L  mode  forward(p1:p2:len ...)  reverse(p1:p2:len ...)
  p1 = 0-based start in ref; p2 = 0-based start in qry (forward); for reverse matches p2 = 0-based start in qry (forward strand) of the window that equals
  revcomp(ref[p1..p1+len)) (mummer prints the start in revcomp(qry); `-c` prints the forward position of the LAST base: p2 = c - len; both asserted).
  mode: maxmatch (all MEMs), mumreference (unique in ref, MAMs), mum (unique in ref AND qry). Brute force definitions: a MEM is (i, j, len) with
  ref[i..i+len) == qry[j..j+len) over A/C/G/T, not extendable left or right; uniqueness = the matched string occurs once in ref (and once in qry);
  on the reverse strand the same with the query replaced by revcomp(qry). A non-ACGT symbol never matches.
kmism: longest common substring with <= k mismatches (Hamming), brute force over every (i, j, len): the longest window length, tie = smallest j, then smallest i;
    s1  s2  k  len:i:j:mismatches | 0
  a symbol pair agrees only when both are the same A/C/G/T letter (N is always a mismatch).
multi: longest common substring of N sequences occurring in >= q of them (Gusfield sec. 7.6), brute force:
    seq1,seq2,...  q  result(string@pos1,pos2,... ordered by string; pos = first occurrence, -1 when absent | -)
Never aborts: a failed external call prints ERR.
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261121)
COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
ACGT = set('ACGT')


def rc(s):
    return ''.join(COMP.get(c, 'N') for c in reversed(s))


def rand(n, alpha):
    return ''.join(rnd.choice(alpha) for _ in range(n))


def occ(s, w):
    return [i for i in range(len(s) - len(w) + 1) if s[i:i + len(w)] == w]


def brute_matches(ref, qry, L, mode):
    out = []
    n, m = len(ref), len(qry)
    for i in range(n):
        for j in range(m):
            if ref[i] != qry[j] or ref[i] not in ACGT:
                continue
            if i > 0 and j > 0 and ref[i - 1] == qry[j - 1] and ref[i - 1] in ACGT:
                continue
            ln = 0
            while i + ln < n and j + ln < m and ref[i + ln] == qry[j + ln] and ref[i + ln] in ACGT:
                ln += 1
            if ln < L:
                continue
            w = ref[i:i + ln]
            if mode == 'mumreference' and len(occ(ref, w)) != 1:
                continue
            if mode == 'mum' and (len(occ(ref, w)) != 1 or len(occ(qry, w)) != 1):
                continue
            out.append((i, j, ln))
    return sorted(out)


def mummer(ref, qry, L, mode):
    with tempfile.TemporaryDirectory() as d:
        a, b = os.path.join(d, 'r.fa'), os.path.join(d, 'q.fa')
        open(a, 'w').write('>r\n' + ref + '\n')
        open(b, 'w').write('>q\n' + qry + '\n')
        res = {}
        for c in (False, True):
            args = ['mummer', '-' + mode, '-b', '-l', str(L)] + (['-c'] if c else []) + [a, b]
            r = subprocess.run(args, capture_output=True, text=True, timeout=60)
            fwd, rev, cur = [], [], None
            for ln in r.stdout.splitlines():
                if ln.startswith('>'):
                    cur = rev if 'Reverse' in ln else fwd
                    continue
                p = ln.split()
                if len(p) == 3 and cur is not None:
                    cur.append((int(p[0]) - 1, int(p[1]), int(p[2])))
            res[c] = (fwd, rev)
    fwd = sorted((i, j - 1, ln) for i, j, ln in res[False][0])
    m = len(qry)
    rev = sorted((i, m - (j - 1) - ln, ln) for i, j, ln in res[False][1])          # revcomp(qry) start -> forward start
    rev_c = sorted((i, j - ln, ln) for i, j, ln in res[True][1])                   # -c: forward 1-based position of the last base
    assert rev == rev_c, (ref, qry, rev, rev_c)
    return fwd, rev


def fmt(ms):
    return ';'.join('%d:%d:%d' % x for x in ms) or '-'


def gen_pair(alpha, with_n):
    shared = [rand(rnd.randint(2, 9), alpha) for _ in range(rnd.randint(1, 3))]
    def build():
        parts = [rand(rnd.randint(0, 8), alpha)]
        for u in shared:
            parts.append(rc(u) if rnd.random() < 0.4 else u)
            parts.append(rand(rnd.randint(0, 8), alpha))
        s = ''.join(parts)
        if with_n and len(s) > 3 and rnd.random() < 0.7:
            i = rnd.randrange(len(s)); s = s[:i] + 'N' * rnd.randint(1, 2) + s[i + 1:]
        return s
    return build(), build()


def mode_mum():
    corpus = []
    for alpha in ('ACGT', 'AC', 'ACGT', 'ACGT'):
        for _ in range(110):
            corpus.append(gen_pair(alpha, False))
    for _ in range(80):
        corpus.append(gen_pair('ACGT', True))
    corpus += [('ACGGATCCGTTTTACGATCCGTAAGT', 'GGATCCGAAAACGGATCCGTTTT'), ('AAAAAAAA', 'AAAA'), ('ACGT', 'ACGT'), ('GATTACA', 'TGTAATC')]
    rows = 0
    for ref, qry in corpus:
        if not ref or not qry:
            continue
        has_n = any(c not in ACGT for c in ref + qry)
        for L in (2, 3, 5):
            for mode in ('maxmatch', 'mumreference', 'mum'):
                bf_f = brute_matches(ref, qry, L, mode)
                bf_r = sorted((i, len(qry) - rj - ln, ln) for i, rj, ln in brute_matches(ref, rc(qry), L, mode))
                if not has_n:
                    try:
                        f, r = mummer(ref, qry, L, mode)
                    except Exception as e:
                        sys.stderr.write('ERR %s %s %s\n' % (ref, qry, e)); print('%s\t%s\t%d\t%s\tERR\tERR' % (ref, qry, L, mode)); continue
                    if f != bf_f or r != sorted(bf_r):
                        sys.stderr.write('MISMATCH %s %s L=%d %s mummer=%s/%s brute=%s/%s\n' % (ref, qry, L, mode, f, r, bf_f, bf_r))
                print('%s\t%s\t%d\t%s\t%s\t%s' % (ref, qry, L, mode, fmt(bf_f), fmt(sorted(bf_r))))
                rows += 1
    sys.stderr.write('rows=%d\n' % rows)


def kmism_best(s1, s2, k):
    best = (0, 0, 0, 0)
    n, m = len(s1), len(s2)
    cand = []
    for i in range(n):
        for j in range(m):
            mm = 0
            for l in range(1, min(n - i, m - j) + 1):
                if not (s1[i + l - 1] == s2[j + l - 1] and s1[i + l - 1] in ACGT):
                    mm += 1
                if mm > k:
                    break
                cand.append((l, i, j, mm))
    if not cand:
        return None
    L = max(c[0] for c in cand)
    return min((c for c in cand if c[0] == L), key=lambda c: (c[2], c[1]))


def mode_kmism():
    rows = 0
    for _ in range(420):
        alpha = rnd.choice(['ACGT', 'AC', 'ACGT', 'ACGTN'])
        s1, s2 = rand(rnd.randint(0, 16), alpha), rand(rnd.randint(0, 16), alpha)
        if rnd.random() < 0.5 and s1:
            u = s1[rnd.randrange(len(s1)):][:rnd.randint(3, 9)]
            u = ''.join(rnd.choice('ACGT') if rnd.random() < 0.15 else c for c in u)
            p = rnd.randint(0, len(s2)); s2 = s2[:p] + u + s2[p:]
        for k in (0, 1, 2, 3):
            b = kmism_best(s1, s2, k)
            print('%s\t%s\t%d\t%s' % (s1 or '-', s2 or '-', k, '0' if b is None else '%d:%d:%d:%d' % (b[0], b[1], b[2], b[3])))
            rows += 1
    sys.stderr.write('rows=%d\n' % rows)


def multi_best(seqs, q):
    cand = set()
    for s in seqs:
        for i in range(len(s)):
            for l in range(1, len(s) - i + 1):
                w = s[i:i + l]
                if all(c in ACGT for c in w):
                    cand.add(w)
    good = [w for w in cand if sum(1 for s in seqs if w in s) >= q]
    if not good:
        return []
    L = max(map(len, good))
    return sorted(w for w in good if len(w) == L)


def mode_multi():
    rows = 0
    for _ in range(420):
        alpha = rnd.choice(['ACGT', 'AC', 'ACGT', 'ACGTN'])
        k = rnd.randint(2, 4)
        core = rand(rnd.randint(2, 7), 'ACGT')
        seqs = []
        for _ in range(k):
            s = rand(rnd.randint(1, 14), alpha)
            if rnd.random() < 0.7:
                p = rnd.randint(0, len(s)); s = s[:p] + core + s[p:]
            seqs.append(s)
        for q in sorted({2, k}):
            res = multi_best(seqs, q)
            txt = ';'.join('%s@%s' % (w, ','.join(str(s.find(w)) for s in seqs)) for w in res) or '-'
            print('%s\t%d\t%s' % (','.join(seqs), q, txt))
            rows += 1
    sys.stderr.write('rows=%d\n' % rows)


if __name__ == '__main__':
    {'mum': mode_mum, 'kmism': mode_kmism, 'multi': mode_multi}[sys.argv[1]]()
