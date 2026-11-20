#!/usr/bin/env python3
"""A5-10 (B09 F174) oracle: Mash 2.3 `bounds`, `dist` (-t, -v, -d, -i, -n, -k, -s, -S) and `triangle` (-E, -v, -d, -i) as TEXT, byte for byte.

JSON lines: {"cmd": "bounds|dist|triangle", "k", "s", "seed", "noncanonical", "individual", "table", "edge", "v", "d",
             "refs": [{"name", "records": [{"id", "comment", "seq"}]}], "queries": [...], "expected": "<stdout>"}
(bounds: only "k" and "p".)
Group "error" (audit A6-6, written LAST so the rows above are unchanged for the committed seed): {"cmd": "error", "sub": "dist|triangle", options as above, "refs", "queries", "status": <exit status or -signal>, "message": "<the WARNING / ERROR line of stderr>"} -
runs of `mash dist` / `mash triangle` on FASTA files (no sketch file) where an input file has no record of at least k symbols (status 1, or a crash: SIGFPE / SIGABRT with -i), which the generator used to drop (`if r.returncode: continue`). The files are written under the names a1.fa, a2.fa, ... so that the names the program prints are those of the inputs; dist references are sketched first with the same -k -s -S -n -i
(`mash sketch -o refs.msh`: the reference file of `mash dist` is one sketch file or one FASTA file).
Families of related genomes (0.5 - 25 % substitutions / indels), unrelated ones, several records per file, k 4 - 32, sketch sizes 5 - 2000 (smaller and larger than the genome).
Usage: python3 scripts/alignment/generate_mash_commands_oracle.py [per_command=220] > mash_commands_oracle.jsonl      (needs mash 2.3 on PATH)
"""
import json, os, random, subprocess, sys, tempfile

rnd = random.Random(20261211)


def rseq(n): return ''.join(rnd.choice('ACGT') for _ in range(n))


def mutate(s, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.7: out.append(rnd.choice('ACGT'))
        elif r < rate * 0.85: pass
        elif r < rate: out.append(c); out.append(rnd.choice('ACGT'))
        else: out.append(c)
    return ''.join(out)


def make_files(n, prefix):
    base = rseq(rnd.randint(1200, 5000))
    files = []
    for i in range(n):
        rate = rnd.choice([0.005, 0.02, 0.06, 0.12, 0.25]) if rnd.random() < 0.85 else None
        genome = mutate(base, rate) if rate is not None else rseq(rnd.randint(1200, 5000))
        recs = []
        pos = 0
        parts = rnd.choice([1, 1, 2, 3])
        cuts = sorted(rnd.sample(range(1, len(genome)), parts - 1)) if parts > 1 else []
        for j, (a, b) in enumerate(zip([0] + cuts, cuts + [len(genome)])):
            recs.append({'id': 'rec%d_%d' % (i, j), 'comment': rnd.choice(['', 'desc %d' % rnd.randint(0, 99), 'strain X']), 'seq': genome[a:b]})
        files.append({'name': '%s%d.fa' % (prefix, i + 1), 'records': recs})
    return files


def write(d, files):
    for f in files:
        with open(os.path.join(d, f['name']), 'w') as h:
            for r in f['records']:
                h.write('>%s%s\n%s\n' % (r['id'], (' ' + r['comment']) if r['comment'] else '', r['seq']))


def mash(d, args):
    p = subprocess.run(['mash'] + args, cwd=d, capture_output=True, text=True)
    return p


def common(o):
    a = ['-k', str(o['k']), '-s', str(o['s']), '-S', str(o['seed'])]
    if o['noncanonical']: a.append('-n')
    return a


def rand_opts():
    k = rnd.choice([4, 6, 9, 12, 15, 16, 17, 21, 21, 21, 24, 28, 32])
    s = rnd.choice([5, 20, 100, 500, 1000, 1000, 2000])
    return {'k': k, 's': s, 'seed': rnd.choice([42, 42, 42, 1, 7, 123456]), 'noncanonical': rnd.random() < 0.25, 'individual': False, 'table': False, 'edge': False, 'v': None, 'd': None}


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 220
    out = []
    for _ in range(per):                       # bounds: k 1..32 x p
        k = rnd.randint(1, 32); p = rnd.choice([0.5, 0.9, 0.95, 0.99, 0.999, 0.75, 0.6827, rnd.random()])
        r = mash(None, ['bounds', '-k', str(k), '-p', str(p)])
        if r.returncode == 0:
            print(json.dumps({'cmd': 'bounds', 'k': k, 'p': p, 'expected': r.stdout}))
    for cmd in ('dist', 'triangle'):
        made = 0; tries = 0
        while made < per and tries < per * 4:
            tries += 1
            o = rand_opts()
            o['individual'] = rnd.random() < 0.2
            with tempfile.TemporaryDirectory() as d:
                if cmd == 'dist':
                    refs = make_files(rnd.randint(1, 4), 'r'); qs = make_files(rnd.randint(1, 3), 'q')
                    o['table'] = rnd.random() < 0.3
                    if rnd.random() < 0.35: o['v'] = rnd.choice([0.001, 0.05, 0.5, 1e-10])
                    if rnd.random() < 0.35: o['d'] = rnd.choice([0.02, 0.1, 0.3, 0.6])
                    write(d, refs); write(d, qs)
                    sk = mash(d, ['sketch', '-o', 'refs'] + common(o) + (['-i'] if o['individual'] else []) + [f['name'] for f in refs])
                    if sk.returncode: continue
                    args = ['dist'] + (['-t'] if o['table'] else []) + (['-v', str(o['v'])] if o['v'] is not None else []) + (['-d', str(o['d'])] if o['d'] is not None else []) + (['-i'] if o['individual'] else []) + ['refs.msh'] + [f['name'] for f in qs]
                else:
                    refs = make_files(rnd.randint(2, 6), 'a'); qs = []
                    o['edge'] = rnd.random() < 0.3
                    if rnd.random() < 0.25: o['v'] = rnd.choice([0.001, 0.05, 0.5])
                    if rnd.random() < 0.25: o['d'] = rnd.choice([0.05, 0.2, 0.5])
                    write(d, refs)
                    args = ['triangle'] + common(o) + (['-E'] if o['edge'] else []) + (['-v', str(o['v'])] if o['v'] is not None else []) + (['-d', str(o['d'])] if o['d'] is not None else []) + (['-i'] if o['individual'] else []) + [f['name'] for f in refs]
                r = mash(d, args)
                if r.returncode: continue
                made += 1
                print(json.dumps(dict(o, cmd=cmd, refs=refs, queries=qs, expected=r.stdout)))


def error_rows(per):
    """Runs the program refuses (A6-6): every non-zero status is KEPT with its message."""
    made = 0
    while made < per:
        sub = 'dist' if made % 2 == 0 else 'triangle'
        o = rand_opts()
        o['individual'] = rnd.random() < 0.25
        k = o['k']
        def bad_file(name):
            if rnd.random() < 0.3:
                return {'name': name, 'records': []}                                     # an empty file: "Did not find fasta records"
            return {'name': name, 'records': [{'id': 'short%d' % j, 'comment': rnd.choice(['', 'tiny']), 'seq': rseq(rnd.randint(1, k - 1))} for j in range(rnd.randint(1, 3))]}
        if k == 1:
            continue
        if sub == 'dist':
            refs = make_files(1, 'r'); qs = make_files(rnd.randint(1, 2), 'q')
            o['table'] = rnd.random() < 0.3
            if rnd.random() < 0.5: refs[0] = bad_file('r1.fa')
            else: qs[rnd.randrange(len(qs))] = bad_file(qs[0]['name'] if len(qs) == 1 else 'q%d.fa' % (1 + rnd.randrange(len(qs))))
            qs = [dict(q, name='q%d.fa' % (i + 1)) for i, q in enumerate(qs)]
        else:
            refs = make_files(rnd.randint(2, 4), 'a'); qs = []
            o['edge'] = rnd.random() < 0.3
            i = rnd.randrange(len(refs)); refs[i] = bad_file(refs[i]['name'])
        with tempfile.TemporaryDirectory() as d:
            write(d, refs); write(d, qs)
            if sub == 'dist':
                args = ['dist'] + (['-t'] if o['table'] else []) + common(o) + (['-i'] if o['individual'] else []) + [refs[0]['name']] + [f['name'] for f in qs]
            else:
                args = ['triangle'] + common(o) + (['-E'] if o['edge'] else []) + (['-i'] if o['individual'] else []) + [f['name'] for f in refs]
            r = mash(d, args)
        if r.returncode == 0:
            continue                                                                       # -i drops the short records: no error (the main groups cover it)
        lines = [l for l in r.stderr.splitlines() if l.startswith(('WARNING', 'ERROR'))]
        made += 1
        print(json.dumps(dict(o, cmd='error', sub=sub, refs=refs, queries=qs, status=r.returncode, message=lines[0] if lines else '')))


if __name__ == '__main__':
    main()
    error_rows(int(sys.argv[2]) if len(sys.argv) > 2 else 80)
