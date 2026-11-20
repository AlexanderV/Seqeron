#!/usr/bin/env python3
"""ALIGN-MULTI-001 oracle (audit A3-9a): MUMmer 3.23 `mgaps` (clustering of exact matches by diagonal difference and separation) vs NucmerClusters.Cluster.

Rows: a list of exact matches `start1:start2:length` (1-based starts, the format `mummer` prints) and the clusters `mgaps -l L -s S -d D -f F [-e]` prints for them.
  (a) 'mummer' rows: the forward-strand matches of `mummer -maxmatch -n -l 6..12` on a random related pair (shared segments with substitutions / indels / repeats / N, so overlapping matches on
      different diagonals, shared starts and tandem repeats occur), fed unchanged to mgaps;
  (b) 'synthetic' rows: 4-45 random matches (lengths 3-70, same-diagonal overlaps, shared start1 / start2, near diagonals), shuffled.
Options per row: -l 1..65, -s 20..200, -d 0..10, -f 0..0.3, optional -e (extent). Rows with two identical (start1, start2) are not generated (std::sort of equal elements is unspecified).
TSV: kind  l  s  d  f  e(0|1)  matches  clusters   where clusters = '|'-joined clusters, a cluster = ';'-joined  start1:start2:length:adj  (adj = the overlap the match is shortened by, 0 for none; starts / lengths as printed)
Never aborts: a failing call is skipped.
Usage: python3 scripts/alignment/generate_mgaps_oracle.py [n=700] > mgaps_oracle.tsv
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261128)


def rseq(n, alpha='ACGT'):
    return ''.join(rnd.choice(alpha) for _ in range(n))


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)


def related_pair():
    unit = rseq(rnd.randint(5, 30))
    core = rseq(rnd.randint(100, 400))
    seg = lambda: rnd.choice([core, core[rnd.randint(0, 50):], unit * rnd.randint(2, 6), rseq(rnd.randint(20, 100))])
    ref = rseq(rnd.randint(0, 40)) + mutate(seg(), 0.03) + rseq(rnd.randint(5, 60)) + mutate(seg(), 0.05) + rseq(rnd.randint(0, 40)) + (unit * 3 if rnd.random() < 0.4 else '')
    qry = rseq(rnd.randint(0, 40)) + mutate(seg(), 0.03) + rseq(rnd.randint(5, 60)) + mutate(seg(), 0.05) + rseq(rnd.randint(0, 40))
    if rnd.random() < 0.5: qry = mutate(ref, 0.04) + rseq(rnd.randint(0, 30))
    return ref, qry


def run_mgaps(text, l, s, d, f, e):
    cmd = ['mgaps', '-l', str(l), '-s', str(s), '-d', str(d), '-f', str(f)] + (['-e'] if e else [])
    r = subprocess.run(cmd, input=text, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    clusters = []; cur = None
    for line in r.stdout.split('\n'):
        if not line.strip(): continue
        if line.startswith('>') or line.startswith('#'):
            cur = []; clusters.append(cur); continue
        f6 = line.split()
        s1, s2, ln = int(f6[0]), int(f6[1]), int(f6[2])
        adj = 0 if f6[3] == 'none' else -int(f6[3])
        cur.append('%d:%d:%d:%d' % (s1, s2, ln, adj))
    return '|'.join(';'.join(c) for c in clusters if c)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 700
    rows = 0
    with tempfile.TemporaryDirectory() as d:
        for idx in range(n):
            if idx % 2 == 0:
                ref, qry = related_pair()
                open(d + '/r.fa', 'w').write('>ref\n' + ref + '\n'); open(d + '/q.fa', 'w').write('>qry\n' + qry + '\n')
                r = subprocess.run(['mummer', '-maxmatch', '-n', '-l', str(rnd.choice([6, 8, 10, 12])), d + '/r.fa', d + '/q.fa'], capture_output=True, text=True)
                ms = []
                for line in r.stdout.split('\n'):
                    if line.startswith('>'):
                        if 'Reverse' in line: break
                        continue
                    p = line.split()
                    if len(p) == 3: ms.append(tuple(map(int, p)))
                kind = 'mummer'
            else:
                ms = []; seen = set()
                for _ in range(rnd.randint(4, 45)):
                    s1 = rnd.randint(1, 400)
                    if ms and rnd.random() < 0.35:
                        b = rnd.choice(ms); diag = b[1] - b[0] + rnd.choice([0, 0, 1, -1, 3, -4, 12])
                        s1 = max(1, b[0] + rnd.randint(-20, 60)); s2 = s1 + diag
                        if rnd.random() < 0.15: s1 = b[0]; s2 = b[1] + rnd.randint(-8, 30)
                        if rnd.random() < 0.15: s2 = b[1]; s1 = b[0] + rnd.randint(-8, 30)
                    else:
                        s2 = rnd.randint(1, 400)
                    if s2 < 1 or s1 < 1: continue
                    ln = rnd.randint(3, 70)
                    if (s1, s2) in seen: continue
                    seen.add((s1, s2)); ms.append((s1, s2, ln))
                rnd.shuffle(ms)
                kind = 'synthetic'
            if not ms or len(ms) > 400: continue
            if len({(a, b) for a, b, _ in ms}) != len(ms): continue
            l = rnd.choice([1, 10, 30, 65, 65]); s = rnd.choice([20, 90, 90, 200]); dd = rnd.choice([0, 5, 5, 10]); f = rnd.choice([0, 0.12, 0.12, 0.3]); e = rnd.random() < 0.2
            text = '> qry\n' + ''.join('%8d %8d %6d\n' % m for m in ms)
            out = run_mgaps(text, l, s, dd, f, e)
            if out is None:
                sys.stderr.write('FAILED %s\n' % text[:50]); continue
            print('\t'.join([kind, str(l), str(s), str(dd), repr(f), '1' if e else '0', ','.join('%d:%d:%d' % m for m in ms), out]))
            rows += 1
    sys.stderr.write('rows=%d\n' % rows)


main()
