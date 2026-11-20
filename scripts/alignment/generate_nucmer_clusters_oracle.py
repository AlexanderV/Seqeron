#!/usr/bin/env python3
"""ALIGN-MULTI-001 oracle (audit A3-9a): the cluster step of `nucmer` 3.23 (`mummer -mum|-mumreference|-maxmatch -b -l L -n | mgaps -l c -s g -d D -f F`, exactly the pipeline of the nucmer script)
vs GenomicAnalyzer.FindNucmerClusters.

Related random pairs (shared segments with substitutions / indels, a repeat, an inverted segment, N, lower case) - the pair is written to FASTA and the installed mummer + mgaps are piped.
TSV: reference  query  mode(mum|mumreference|maxmatch)  l  c  s  D  d   clusters
  clusters = '|'-joined clusters; a cluster = strand letter F|R, '/', then ';'-joined  start1:start2:length:adj  (adj = overlap cut off the match; the starts are those mgaps prints: 1-based, the second one in the
  coordinates of the reverse-complemented query for R).
Never aborts: a failing pipeline is skipped.
Usage: python3 scripts/alignment/generate_nucmer_clusters_oracle.py [n=600] > nucmer_clusters_oracle.tsv
"""
import os, random, subprocess, sys, tempfile

rnd = random.Random(20261129)
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def rseq(n):
    return ''.join(rnd.choice('ACGT') for _ in range(n))


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)


def revcomp(s):
    return s.translate(COMP)[::-1]


def pair():
    seg = [rseq(rnd.randint(40, 250)) for _ in range(rnd.randint(2, 4))]
    unit = rseq(rnd.randint(6, 25))
    ref = rseq(rnd.randint(0, 60))
    qry = rseq(rnd.randint(0, 60))
    for s in seg:
        ref += mutate(s, 0.02) + rseq(rnd.randint(3, 70)) + (unit * rnd.randint(2, 4) if rnd.random() < 0.3 else '')
        piece = mutate(s, rnd.choice([0.0, 0.02, 0.05, 0.09]))
        qry += (revcomp(piece) if rnd.random() < 0.2 else piece) + rseq(rnd.randint(3, 70))
        if rnd.random() < 0.3: qry += unit * rnd.randint(1, 3)
    if rnd.random() < 0.2:
        i = rnd.randrange(len(ref)); ref = ref[:i] + 'N' * rnd.randint(1, 5) + ref[i:]
    return ref, qry


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 600
    rows = 0
    with tempfile.TemporaryDirectory() as d:
        for _ in range(n):
            ref, qry = pair()
            mode = rnd.choice(['mum', 'mum', 'mumreference', 'maxmatch'])
            l = rnd.choice([8, 12, 15, 20]); c = rnd.choice([20, 45, 65, 65]); g = rnd.choice([40, 90, 90, 200]); D = rnd.choice([2, 5, 5, 10]); f = rnd.choice([0.05, 0.12, 0.12, 0.3])
            open(d + '/r.fa', 'w').write('>ref\n' + ref + '\n'); open(d + '/q.fa', 'w').write('>qry\n' + qry + '\n')
            m = subprocess.run(['mummer', '-' + mode, '-b', '-l', str(l), '-n', d + '/r.fa', d + '/q.fa'], capture_output=True, text=True)
            if m.returncode != 0: continue
            g_ = subprocess.run(['mgaps', '-l', str(c), '-s', str(g), '-d', str(D), '-f', str(f)], input=m.stdout, capture_output=True, text=True)
            if g_.returncode != 0: continue
            clusters = []; cur = None; strand = 'F'
            for line in g_.stdout.split('\n'):
                if not line.strip(): continue
                if line.startswith('>'):
                    strand = 'R' if 'Reverse' in line else 'F'
                    cur = []; clusters.append((strand, cur)); continue
                if line.startswith('#'):
                    cur = []; clusters.append((strand, cur)); continue
                p = line.split()
                cur.append('%s:%s:%s:%d' % (p[0], p[1], p[2], 0 if p[3] == 'none' else -int(p[3])))
            out = '|'.join(s + '/' + ';'.join(cl) for s, cl in clusters if cl)
            print('\t'.join([ref, qry, mode, str(l), str(c), str(g), str(D), repr(f), out]))
            rows += 1
    sys.stderr.write('rows=%d\n' % rows)


main()
