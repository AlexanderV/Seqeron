#!/usr/bin/env python3
"""A5-10b3 (B09 F194) oracle: the .msh file of READ SETS (`mash sketch -r -o x.msh`), whose `counts32` are listed in the iteration order of Mash's bundled robin_hood::unordered_map.
Same JSON schema as generate_mash_msh_oracle.py (without `info`). Read sets of 10 - 1 200 reads of 50 / 100 bp drawn from a 500 - 20 000 bp source, k 4 - 32, sketch sizes 1 - 3 000 (tables of 8 - 4 096 buckets
with evictions: erase = backward shift), -m 1 - 3, -c, -g, -n, seeds.
Usage: python3 scripts/alignment/generate_mash_reads_msh_oracle.py OUT.jsonl.gz [n=140] [seed]      (needs mash 2.3 on PATH)
"""
import base64, gzip, json, os, random, subprocess, sys, tempfile


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            k = rnd.choice([4, 8, 12, 16, 17, 21, 21, 24, 32]); s = rnd.choice([1, 5, 50, 500, 1000, 2000, 3000]); sd = rnd.choice([42, 1, 12345, 0, 777])
            glen = rnd.choice([500, 2000, 6000, 20000])
            g = ''.join(rnd.choice('ACGT') for _ in range(glen))
            files = []
            for i in range(rnd.choice([1, 1, 2])):
                recs = []
                for j in range(rnd.randint(10, 1200)):
                    st = rnd.randint(0, max(0, len(g) - 100)); r = g[st:st + rnd.choice([50, 100])]
                    if rnd.random() < 0.3: r = ''.join(rnd.choice('ACGT') if rnd.random() < 0.02 else c for c in r)
                    recs.append(('r%d_%d' % (i, j), rnd.choice(['', 'lane 1']), r))
                files.append({'name': 'reads%d.fq.fa' % i, 'fasta': ''.join('>%s%s\n%s\n' % (a, (' ' + b) if b else '', c) for a, b, c in recs)})
            args = ['-k', str(k), '-s', str(s), '-S', str(sd), '-r', '-m', str(rnd.choice([1, 1, 2, 3]))]
            if rnd.random() < 0.2: args.append('-n')
            if rnd.random() < 0.3: args += ['-c', str(rnd.choice([1, 3, 8]))]
            if rnd.random() < 0.3: args += ['-g', str(rnd.choice([1000, 100000]))]
            names = [f['name'] for f in files]
            for f in files:
                open(os.path.join(d, f['name']), 'w', newline='\n').write(f['fasta'])
            o = os.path.join(d, 'o')
            if os.path.exists(o + '.msh'): os.remove(o + '.msh')
            p = subprocess.run(['mash', 'sketch', '-o', 'o'] + args + names, cwd=d, capture_output=True, text=True, timeout=120)
            if p.returncode != 0 or not os.path.exists(o + '.msh'):
                continue
            rows.append({'id': len(rows), 'args': args + names, 'files': files, 'list': None, 'msh': base64.b64encode(open(o + '.msh', 'rb').read()).decode(), 'info': {}})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 140, int(sys.argv[3]) if len(sys.argv) > 3 else 20261121)
