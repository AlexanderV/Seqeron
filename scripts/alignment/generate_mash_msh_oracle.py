#!/usr/bin/env python3
"""A5-10b (B09 F193) oracle: the Mash 2.3 sketch FILE (`mash sketch -o x.msh`, a Cap'n Proto message of MinHash.capnp built by capnp::MallocMessageBuilder) and the text of `mash info` on it
(default, -t, -H, -d, -c). Genomes (1 - 3 files of 1 - 4 records, 20 bp - 60 kb with N / lower case) and read sets (-r -m -c -g: counts32), k 4 - 32, sketch sizes 1 - 3 000 (64-bit sketches of 1 000 and more
hashes need several segments and far pointers), seeds (42 = referenceListOld, others = referenceList), -n, -i, -l (a list file), -Z is not used here (alphabets are registered separately).
JSON lines: {"id", "args": [...mash sketch arguments after the program name, files as typed], "files": [{"name", "fasta"}], "list": <text of the -l list file or null>, "msh": "<base64>", "info": {"default",
             "t", "H", "d", "c"} (stdout, or the text of stderr: "ERROR: ...")}
Usage: python3 scripts/alignment/generate_mash_msh_oracle.py OUT.jsonl.gz [n=120] [seed]      (needs mash 2.3 on PATH)
"""
import base64, gzip, json, os, random, subprocess, sys, tempfile


def main(out, n, seed):
    rnd = random.Random(seed)

    def rseq(m):
        s = ''.join(rnd.choice('ACGT') for _ in range(m))
        if rnd.random() < 0.2:
            s = ''.join(c.lower() if rnd.random() < 0.3 else c for c in s)
        if rnd.random() < 0.15:
            s = ''.join('N' if rnd.random() < 0.02 else c for c in s)
        return s

    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            reads = rnd.random() < 0.3
            k = rnd.choice([4, 8, 12, 16, 17, 21, 21, 24, 32]); s = rnd.choice([1, 5, 50, 500, 500, 1000, 1000, 2000, 3000]); sd = rnd.choice([42, 42, 42, 1, 12345, 4294967295, 0])
            files = []
            for i in range(rnd.choice([1, 1, 2, 3]) if s < 1000 else rnd.choice([1, 1, 2])):
                recs = []
                if reads:
                    g = rseq(rnd.randint(500, 3000))
                    for j in range(rnd.randint(10, 150)):
                        st = rnd.randint(0, max(0, len(g) - 100)); recs.append(('r%d_%d' % (i, j), rnd.choice(['', 'lane 1']), g[st:st + rnd.choice([50, 100])]))
                else:
                    for j in range(rnd.choice([1, 1, 2, 3]) if s < 1000 else 1):
                        L = rnd.choice([3000, 6000, 12000]) if s >= 1000 else rnd.choice([20, 200, 2000, 6000])
                        recs.append(('seq%d_%d' % (i, j), rnd.choice(['', '', 'desc %d' % j, 'a longer comment with words']), rseq(L)))
                files.append({'name': rnd.choice(['f%d.fa' % i, 'genome_%d.fasta' % i, 'x%d' % i]), 'fasta': ''.join('>%s%s\n%s\n' % (a, (' ' + b) if b else '', c) for a, b, c in recs)})
            names = [f['name'] for f in files]
            if len(set(names)) != len(names): continue
            args = ['-k', str(k), '-s', str(s), '-S', str(sd)]
            if rnd.random() < 0.2: args.append('-n')
            if not reads and rnd.random() < 0.25: args.append('-i')
            if rnd.random() < 0.15: args += ['-I', rnd.choice(['myid', 'sample 1', 'x'])]
            if rnd.random() < 0.15: args += ['-C', rnd.choice(['a comment', 'c', 'two words here'])]
            if reads:
                args += ['-r', '-m', str(rnd.choice([1, 1, 2])), ] + (['-c', str(rnd.choice([1, 3, 8]))] if rnd.random() < 0.3 else []) + (['-g', str(rnd.choice([1000, 100000]))] if rnd.random() < 0.3 else [])
            for f in files:
                open(os.path.join(d, f['name']), 'w', newline='\n').write(f['fasta'])
            listtext = None
            if rnd.random() < 0.15:
                listtext = ''.join(nm + '\n' for nm in names); open(os.path.join(d, 'list.txt'), 'w').write(listtext)
                cmdfiles = ['-l', 'list.txt']
            else:
                cmdfiles = names
            o = os.path.join(d, 'o')
            if os.path.exists(o + '.msh'): os.remove(o + '.msh')
            p = subprocess.run(['mash', 'sketch', '-o', 'o'] + args + cmdfiles, cwd=d, capture_output=True, text=True, timeout=120)
            if p.returncode != 0 or not os.path.exists(o + '.msh'):
                continue
            msh = open(o + '.msh', 'rb').read()
            info = {}
            for key, opt in (('default', []), ('t', ['-t']), ('H', ['-H']), ('d', ['-d']), ('c', ['-c'])):
                q = subprocess.run(['mash', 'info'] + opt + ['o.msh'], cwd=d, capture_output=True, text=True, timeout=120)
                info[key] = q.stdout if q.returncode == 0 else q.stderr.strip()
            rows.append({'id': len(rows), 'args': args + cmdfiles, 'files': files, 'list': listtext, 'msh': base64.b64encode(msh).decode(), 'info': info})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 120, int(sys.argv[3]) if len(sys.argv) > 3 else 20261222)
