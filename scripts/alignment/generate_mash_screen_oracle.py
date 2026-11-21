#!/usr/bin/env python3
"""A5-10b (B09 F185) oracle: Mash 2.3 `screen` (query sketch from `mash sketch [-k -s -S -n -i]`, mixtures of reads / contigs, -w winner-takes-all, -i minimum identity, -v maximum p-value) as TEXT, byte for byte.

JSON lines: {"k", "s", "seed", "noncanonical", "individual", "winning", "identity", "pvalue", "queries": [{"name", "records": [{"id", "comment", "seq"}]}], "mixtures": [...same...], "expected": "<stdout>"}
Families of related genomes as queries; mixtures of error-bearing reads of several of them (and random reads), or whole contigs; k 4 - 32, sketch sizes 5 - 2000.
Usage: python3 scripts/alignment/generate_mash_screen_oracle.py OUT.jsonl.gz [n=220] [seed]      (needs mash 2.3 on PATH)
"""
import gzip, json, os, random, subprocess, sys, tempfile


def main(out, n, seed):
    rnd = random.Random(seed)

    def rseq(m): return ''.join(rnd.choice('ACGT') for _ in range(m))

    def mutate(s, rate):
        o = []
        for c in s:
            r = rnd.random()
            if r < rate * 0.7: o.append(rnd.choice('ACGT'))
            elif r < rate * 0.85: pass
            elif r < rate: o.append(c); o.append(rnd.choice('ACGT'))
            else: o.append(c)
        return ''.join(o)

    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            base = rseq(rnd.randint(800, 3500))
            nq = rnd.choice([1, 2, 3, 4])
            genomes = []
            for i in range(nq):
                genomes.append(mutate(base, rnd.choice([0.005, 0.03, 0.1, 0.2])) if rnd.random() < 0.8 else rseq(rnd.randint(800, 3500)))
            queries = []
            for i, g in enumerate(genomes):
                recs = []
                parts = rnd.choice([1, 1, 2, 3])
                cuts = sorted(rnd.sample(range(1, len(g)), parts - 1)) if parts > 1 else []
                for j, (a, b) in enumerate(zip([0] + cuts, cuts + [len(g)])):
                    recs.append({'id': 'q%d_%d' % (i, j), 'comment': rnd.choice(['', 'desc %d' % rnd.randint(0, 99), 'strain X']), 'seq': g[a:b]})
                queries.append({'name': 'q%d.fa' % (i + 1), 'records': recs})
            mixtures = []
            for m in range(rnd.choice([1, 1, 2])):
                recs = []
                if rnd.random() < 0.75:      # reads
                    rl = rnd.choice([60, 100, 150, 250]); err = rnd.choice([0.0, 0.01, 0.05])
                    for r in range(rnd.randint(20, 160)):
                        src = rnd.choice(genomes) if rnd.random() < 0.85 else rseq(rnd.randint(300, 1000))
                        st = rnd.randint(0, max(0, len(src) - rl)); read = src[st:st + rl]
                        if rnd.random() < 0.5: read = read.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
                        recs.append({'id': 'r%d_%d' % (m, r), 'comment': '', 'seq': mutate(read, err) if err else read})
                else:                         # contigs
                    for r in range(rnd.randint(1, 4)):
                        recs.append({'id': 'c%d_%d' % (m, r), 'comment': 'contig', 'seq': rnd.choice(genomes) if rnd.random() < 0.7 else rseq(rnd.randint(500, 2000))})
                mixtures.append({'name': 'm%d.fa' % (m + 1), 'records': recs})
            k = rnd.choice([4, 8, 12, 16, 17, 21, 21, 21, 24, 32]); s = rnd.choice([5, 20, 100, 500, 1000, 2000]); seed_ = rnd.choice([42, 42, 42, 0, 7, 1234567])
            nonc = rnd.random() < 0.15; indiv = rnd.random() < 0.2; win = rnd.random() < 0.35
            ident = rnd.choice([0.0, 0.0, 0.0, -1, 0.5, 0.9, 0.99]); pv = rnd.choice([1.0, 1.0, 1.0, 0.05, 1e-6])
            for grp in (queries, mixtures):
                for f in grp:
                    with open(os.path.join(d, f['name']), 'w') as h:
                        for r in f['records']: h.write('>%s%s\n%s\n' % (r['id'], (' ' + r['comment']) if r['comment'] else '', r['seq']))
            sk = ['sketch', '-o', 'q', '-k', str(k), '-s', str(s), '-S', str(seed_)] + (['-n'] if nonc else []) + (['-i'] if indiv else []) + [f['name'] for f in queries]
            p = subprocess.run(['mash'] + sk, cwd=d, capture_output=True, text=True)
            if p.returncode != 0 or not os.path.exists(os.path.join(d, 'q.msh')): continue
            sc = ['screen'] + (['-w'] if win else []) + ['-i', str(ident), '-v', str(pv), 'q.msh'] + [f['name'] for f in mixtures]
            p = subprocess.run(['mash'] + sc, cwd=d, capture_output=True, text=True)
            if p.returncode != 0: continue
            rows.append({'k': k, 's': s, 'seed': seed_, 'noncanonical': nonc, 'individual': indiv, 'winning': win, 'identity': ident, 'pvalue': pv, 'queries': queries, 'mixtures': mixtures, 'expected': p.stdout})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 220, int(sys.argv[3]) if len(sys.argv) > 3 else 20261223)
