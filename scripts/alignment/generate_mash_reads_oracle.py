#!/usr/bin/env python3
"""A5-10b (B09 F185) oracle: the READ-SET mode and the comment option of Mash 2.3 `dist` / `triangle`: `-r` with `-m` (minimum copies: the noise filter), `-c` (target coverage), `-g` (genome size) and `-C`
(comments with the names), as TEXT. Read sets of 1 - 12x coverage with substitution / indel errors from related genomes (and genomes themselves for `-C`), k 4 - 32, sketch sizes 5 - 2000, seeds, `-n`, `-v`, `-d`, `-t`.
The Bloom filter (`-b`) is not part of the port (registered).

JSON lines: {"cmd": "dist|triangle", "k", "s", "seed", "noncanonical", "reads": bool, "minCopies", "targetCoverage", "genomeSize", "comments": bool, "table", "edge", "v", "d",
             "refs": [{"name", "records": [{"id", "comment", "seq"}]}], "queries": [...], "expected": "<stdout>"}   (triangle: only "refs" holds the inputs)
Usage: python3 scripts/alignment/generate_mash_reads_oracle.py OUT.jsonl.gz [n=240] [seed]      (needs mash 2.3 on PATH)
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

    def readset(name, genome, cov, rl, err):
        recs = []
        for i in range(max(3, int(cov * len(genome) / rl))):
            st = rnd.randint(0, max(0, len(genome) - rl)); r = genome[st:st + rl]
            if rnd.random() < 0.5: r = r.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
            recs.append({'id': '%s_r%d' % (name[:-3], i), 'comment': rnd.choice(['', 'lane %d' % rnd.randint(1, 4)]), 'seq': mutate(r, err) if err else r})
        return {'name': name, 'records': recs}

    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            base = rseq(rnd.randint(800, 3000))
            nfiles = rnd.choice([2, 3, 4])
            cmd = rnd.choice(['dist', 'dist', 'triangle'])
            files = []
            for i in range(nfiles):
                g = mutate(base, rnd.choice([0.005, 0.03, 0.1])) if rnd.random() < 0.85 else rseq(rnd.randint(800, 3000))
                files.append(readset('f%d.fa' % (i + 1), g, rnd.choice([1, 2, 3, 6, 12]), rnd.choice([50, 100, 150]), rnd.choice([0.0, 0.01, 0.03])))
            k = rnd.choice([4, 8, 12, 16, 17, 21, 21, 24, 32]); s = rnd.choice([5, 20, 100, 500, 1000, 2000]); seed_ = rnd.choice([42, 42, 0, 7]); nonc = rnd.random() < 0.15
            reads = rnd.random() < 0.85; m = rnd.choice([1, 1, 2, 3, 5]) if reads else 1; c = rnd.choice([0, 0, 0, 1, 2, 4, 8]) if reads else 0
            gsz = rnd.choice([0, 0, 0, 1000, 5000, 100000]) if reads else 0; comments = rnd.random() < 0.3
            table = cmd == 'dist' and rnd.random() < 0.2; edge = cmd == 'triangle' and rnd.random() < 0.3
            v = rnd.choice([1.0, 1.0, 0.05]); dd = rnd.choice([1.0, 1.0, 0.5])
            for f in files:
                with open(os.path.join(d, f['name']), 'w') as h:
                    for r in f['records']: h.write('>%s%s\n%s\n' % (r['id'], (' ' + r['comment']) if r['comment'] else '', r['seq']))
            opts = ['-k', str(k), '-s', str(s), '-S', str(seed_)] + (['-v', str(v)] if v < 1 else []) + (['-d', str(dd)] if dd < 1 else []) + (['-n'] if nonc else []) + (['-r', '-m', str(m)] + (['-c', str(c)] if c else []) + (['-g', str(gsz)] if gsz else []) if reads else []) + (['-C'] if comments else [])
            if cmd == 'dist':
                refs, queries = files[:1], files[1:]
                args = ['dist'] + (['-t'] if table else []) + opts + [f['name'] for f in refs] + [f['name'] for f in queries]
            else:
                refs, queries = files, []
                args = ['triangle'] + (['-E'] if edge else []) + opts + [f['name'] for f in refs]
            p = subprocess.run(['mash'] + args, cwd=d, capture_output=True, text=True)
            if p.returncode != 0: continue
            rows.append({'cmd': cmd, 'k': k, 's': s, 'seed': seed_, 'noncanonical': nonc, 'reads': reads, 'minCopies': m, 'targetCoverage': c, 'genomeSize': gsz, 'comments': comments, 'table': table, 'edge': edge, 'v': v, 'd': dd,
                         'refs': refs, 'queries': queries, 'expected': p.stdout})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 240, int(sys.argv[3]) if len(sys.argv) > 3 else 20261224)
