#!/usr/bin/env python3
"""A5-10b (B09 F193) oracle: `mash paste` and `mash dist` / `triangle` / `screen` with SKETCH FILES (.msh) as inputs, byte for byte.
Related genomes are sketched with `mash sketch` (random -k -s -S -n, -i); the commands then take the sketch files and/or the sequence files: paste (file lists with equal / different / larger / smaller sketch
sizes, other seeds / k / strand modes: the program skips the incompatible ones and reduces the larger ones), dist (reference .msh, queries .msh and / or sequences; reference a sequence file and queries .msh), triangle
(.msh inputs, possibly with several sketches per file), screen (query .msh, mixtures of reads).
JSON lines: {"cmd", "args": [...options typed], "inputs": [{"name", "kind": "seq"|"msh", "fasta": text | "msh": base64}], "expected": "<stdout>" | "paste": "<base64 of the written file>"}
Usage: python3 scripts/alignment/generate_mash_msh_commands_oracle.py OUT.jsonl.gz [n=160] [seed]      (needs mash 2.3 on PATH)
"""
import base64, gzip, json, os, random, subprocess, sys, tempfile


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

    def run(d, *a):
        return subprocess.run(['mash'] + list(a), cwd=d, capture_output=True, text=True, timeout=300)

    rows = []
    with tempfile.TemporaryDirectory() as d:
        while len(rows) < n:
            base = rseq(rnd.randint(800, 3000))
            nfiles = rnd.choice([2, 3, 4])
            fastas = []
            for i in range(nfiles):
                g = mutate(base, rnd.choice([0.005, 0.03, 0.1])) if rnd.random() < 0.85 else rseq(rnd.randint(800, 3000))
                recs = [('s%d_%d' % (i, j), rnd.choice(['', 'desc %d' % j]), g[j * len(g) // 2:(j + 1) * len(g) // 2]) for j in range(rnd.choice([1, 2]))]
                fastas.append(''.join('>%s%s\n%s\n' % (a, (' ' + b) if b else '', c) for a, b, c in recs))
            k = rnd.choice([8, 12, 16, 17, 21, 24]); s = rnd.choice([20, 100, 500, 1000]); sd = rnd.choice([42, 42, 42, 7])
            for i, f in enumerate(fastas): open(os.path.join(d, 'g%d.fa' % i), 'w').write(f)
            cmd = rnd.choice(['paste', 'dist', 'dist', 'triangle', 'screen'])
            # sketch files: mostly the common parameters, sometimes another size / seed / k / strand mode
            def sketch(i, kk, ss, sdd, nn, indiv):
                o = 'g%d_%d_%d_%d_%d_%d' % (i, kk, ss, sdd, nn, indiv)
                p = run(d, 'sketch', '-o', o, '-k', str(kk), '-s', str(ss), '-S', str(sdd), *(['-n'] if nn else []), *(['-i'] if indiv else []), 'g%d.fa' % i)
                return o + '.msh' if p.returncode == 0 else None
            def variant(first):
                if first or rnd.random() < 0.6: return (k, s, sd, 0)
                r = rnd.random()
                if r < 0.3: return (k, rnd.choice([10, 20, 50, 2000]), sd, 0)
                if r < 0.5: return (k, s, rnd.choice([1, 99]), 0)
                if r < 0.7: return (rnd.choice([11, 14, 19]), s, sd, 0)
                return (k, s, sd, 1)
            if cmd == 'paste':
                names = []
                for i in range(nfiles):
                    kk, ss, sdd, nn = variant(i == 0)
                    m = sketch(i, kk, ss, sdd, nn, rnd.random() < 0.2 and i > 0)
                    if m: names.append(m)
                if len(names) < 2: continue
                if os.path.exists(os.path.join(d, 'pasted.msh')): os.remove(os.path.join(d, 'pasted.msh'))
                p = run(d, 'paste', 'pasted', *names)
                if p.returncode != 0 or not os.path.exists(os.path.join(d, 'pasted.msh')): continue
                rows.append({'cmd': 'paste', 'args': [], 'inputs': [{'name': m, 'kind': 'msh', 'msh': base64.b64encode(open(os.path.join(d, m), 'rb').read()).decode()} for m in names],
                             'paste': base64.b64encode(open(os.path.join(d, 'pasted.msh'), 'rb').read()).decode()})
                continue
            inputs = []; names = []
            for i in range(nfiles):
                use_msh = rnd.random() < (0.9 if i == 0 and cmd != 'dist' else 0.6) or (i == 0 and cmd in ('triangle', 'screen'))
                if use_msh:
                    kk, ss, sdd, nn = variant(i == 0)
                    m = sketch(i, kk, ss, sdd, nn, rnd.random() < 0.15)
                    if m is None: continue
                    names.append(m); inputs.append({'name': m, 'kind': 'msh', 'msh': base64.b64encode(open(os.path.join(d, m), 'rb').read()).decode()})
                else:
                    names.append('g%d.fa' % i); inputs.append({'name': 'g%d.fa' % i, 'kind': 'seq', 'fasta': fastas[i]})
            if len(inputs) < 2: continue
            opts = []
            first_is_msh = inputs[0]['kind'] == 'msh'
            if not first_is_msh or cmd == 'screen':
                opts += ['-k', str(k), '-s', str(s), '-S', str(sd)]
            if cmd == 'dist':
                if rnd.random() < 0.2: opts += ['-t']
                args = ['dist'] + opts + names[:1] + names[1:]
                if not first_is_msh and False: pass
            elif cmd == 'triangle':
                if rnd.random() < 0.3: opts += ['-E']
                args = ['triangle'] + opts + names
            else:
                # screen: query sketch first, then the mixture (reads of the genomes)
                mix = 'mix.fq'
                with open(os.path.join(d, 'mix.fa'), 'w') as f:
                    for r in range(rnd.randint(30, 100)):
                        src = rnd.choice(fastas).split('\n')[1]; st = rnd.randint(0, max(0, len(src) - 100)); f.write('>r%d\n%s\n' % (r, src[st:st + 100]))
                inputs = inputs[:1]; names = names[:1]
                if inputs[0]['kind'] != 'msh': continue
                args = ['screen'] + [] + names + ['mix.fa']
                inputs.append({'name': 'mix.fa', 'kind': 'seq', 'fasta': open(os.path.join(d, 'mix.fa')).read()})
            p = run(d, *args)
            if p.returncode != 0: continue
            rows.append({'cmd': cmd, 'args': args[1:], 'inputs': inputs, 'expected': p.stdout})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    c = {}
    for r in rows: c[r['cmd']] = c.get(r['cmd'], 0) + 1
    print(len(rows), 'rows', c)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 160, int(sys.argv[3]) if len(sys.argv) > 3 else 20261223)
