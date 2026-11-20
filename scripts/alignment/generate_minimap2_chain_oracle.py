#!/usr/bin/env python3
"""ALIGN-MULTI-001 oracle (audit A3-9c): minimap2 2.26 chaining (mg_lchain_dp + backtrack + compact_a, `lchain.c`) vs Minimap2Chaining.LChainDp.

minimap2 2.26 `--print-seeds --print-chains` prints its own anchor array (SD lines: target, end position on the target, strand, end position on the query - reverse strand in the reverse-complemented
query -, span; sorted as minimap2 sorts them) and the chains it builds from it (CN lines: chain index, then the anchors of the chain). This script runs the installed minimap2 on related random
pairs (shared segments with errors, inversions, tandem repeats, two reference sequences) with random -k -w -g -r b,b (so that the long-join rechaining of minimap2 stays off) -m -n and records, per case, the seed array, the options and
  chains     the anchors of every chain as seed indices (anchors of the chain in the order minimap2 prints them), chains sorted by their first seed
  paf        for every line of the plain-PAF output (primary / secondary): (cm, s1) = number of anchors and chain score
JSON lines: {"k":..., "w":..., "g":..., "bw":..., "m":..., "n":..., "targets": [...], "query": "...", "seeds": [[tid, strand, tpos, qpos, span], ...], "chains": [[seed index, ...], ...], "paf": [[cm, s1], ...]}
  (the target and query sequences let the test run the whole seeding: minimizer sketch, index, seed collection, then chaining)
Cases with a repeat-triggered rechaining (RS > 0) or more than 700 seeds are skipped. Never aborts.
Usage: python3 scripts/alignment/generate_minimap2_chain_oracle.py [n=500] > minimap2_chain_oracle.jsonl
"""
import json, os, random, subprocess, sys, tempfile

rnd = random.Random(20261130)
COMP = str.maketrans('ACGT', 'TGCA')


def rseq(n):
    return ''.join(rnd.choice('ACGT') for _ in range(n))


def revcomp(s):
    return s.translate(COMP)[::-1]


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.4: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.7: pass
        elif r < rate: out.append(c); out.append(rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    rows = 0
    with tempfile.TemporaryDirectory() as d:
        for _ in range(n):
            nref = rnd.choice([1, 1, 2])
            refs = []
            segs = []
            for _r in range(nref):
                parts = []
                for _s in range(rnd.randint(2, 4)):
                    seg = rseq(rnd.randint(60, 400)); segs.append(seg)
                    parts.append(seg + rseq(rnd.randint(10, 150)))
                    if rnd.random() < 0.3: unit = rseq(rnd.randint(4, 20)); parts.append(unit * rnd.randint(3, 8))
                refs.append(''.join(parts))
            qparts = []
            err = rnd.choice([0.01, 0.03, 0.08, 0.15])
            for seg in rnd.sample(segs, k=min(len(segs), rnd.randint(2, 4))):
                piece = mutate(seg, err)
                qparts.append(revcomp(piece) if rnd.random() < 0.25 else piece)
                qparts.append(rseq(rnd.randint(0, 120)))
            qry = ''.join(qparts)
            if len(qry) < 60: continue
            k = rnd.choice([9, 11, 13, 15, 15, 15, 17]); w = rnd.choice([3, 5, 10, 10, 10]); g = rnd.choice([100, 400, 1000, 5000, 5000]); bw = rnd.choice([30, 100, 500, 500]); m = rnd.choice([20, 40, 40, 60]); nn = rnd.choice([2, 3, 3, 5])
            open(d + '/r.fa', 'w').write(''.join('>r%d\n%s\n' % (i, s) for i, s in enumerate(refs))); open(d + '/q.fa', 'w').write('>q\n' + qry + '\n')
            base = ['minimap2', '-k', str(k), '-w', str(w), '-g', str(g), '-r', '%d,%d' % (bw, bw), '-m', str(m), '-n', str(nn)]
            r = subprocess.run(base + ['--print-seeds', '--print-chains', d + '/r.fa', d + '/q.fa'], capture_output=True, text=True)
            if r.returncode != 0: continue
            chains = {}; rep = 0
            for line in r.stderr.split('\n'):
                p = line.split('\t')
                if p[0] == 'RS': rep = int(p[1])
            # the SD line is: SD tname tpos strand qpos span gap
            seeds = []; key = {}
            for line in r.stderr.split('\n'):
                p = line.split('\t')
                if p[0] == 'SD':
                    s = (int(p[1][1:]), 1 if p[3] == '-' else 0, int(p[2]), int(p[4]), int(p[5]))
                    key.setdefault(s, len(seeds)); seeds.append(list(s))
                elif p[0] == 'CN':
                    s = (int(p[2][1:]), 1 if p[4] == '-' else 0, int(p[3]), int(p[5]), int(p[6]))
                    chains.setdefault(int(p[1]), []).append(key[s])
            if rep > 0 or not seeds or len(seeds) > 700: continue
            paf = subprocess.run(base + ['-N', '1000', '-p', '0', d + '/r.fa', d + '/q.fa'], capture_output=True, text=True)
            pl = []
            for line in paf.stdout.split('\n'):
                f = line.split('\t')
                if len(f) < 12: continue
                tags = {t.split(':')[0]: t.split(':')[2] for t in f[12:]}
                pl.append([int(tags['cm']), int(tags['s1'])])
            cl = sorted(chains.values(), key=lambda c: c[0])
            print(json.dumps({'k': k, 'w': w, 'g': g, 'bw': bw, 'm': m, 'n': nn, 'targets': refs, 'query': qry, 'seeds': seeds, 'chains': cl, 'paf': pl}))
            rows += 1
    sys.stderr.write('rows=%d\n' % rows)


main()
