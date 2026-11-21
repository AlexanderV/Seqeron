#!/usr/bin/env python3
"""A3-9f oracle: the file naming of `nucmer -p|--prefix` (MUMmer 3.23): which files the script leaves behind for every option family (<prefix>.delta | .cluster | .coords; the temporaries .ntref / .mgaps are deleted), their texts,
the stderr of the whole pipeline (progress lines, prenuc warnings, mummer's reading lines without its timing lines, postnuc warnings) and the failing runs (a prefix in a directory that does not exist, two matching algorithms).
Rows: {"ref","qry","args","files":{name:text},"err","rc"}. Usage: python3 generate_nucmer_prefix_oracle.py [n=170] [seed=20261206] | gzip > nucmer_prefix_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09pfx'

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 170
    rnd = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 20261206)
    def rseq(k): return ''.join(rnd.choice('ACGT') for _ in range(k))
    def mutate(s, rate):
        return ''.join(rnd.choice('ACGT') if rnd.random() < rate else c for c in s)
    comp = str.maketrans('ACGT', 'TGCA')
    rows = 0
    for case in range(n):
        base = rseq(rnd.randint(150, 700))
        ref = '>ref1 desc\n' + mutate(base, .01) + '\n' + ('>ref2\n' + rseq(rnd.randint(30, 200)) + '\n' if rnd.random() < .4 else '')
        q = mutate(base[rnd.randint(0, 40):], rnd.choice([0, .03, .08]))
        if rnd.random() < .3: q = q.translate(comp)[::-1]
        if rnd.random() < .15: q = q[:len(q) // 2] + 'NN' + 'e' + q[len(q) // 2:]
        qry = '>q1\n' + q + '\n' + ('>q2\n' + rseq(rnd.randint(20, 100)) + '\n' if rnd.random() < .3 else '')
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK + '/sub')
        open(WORK + '/ref.fa', 'w').write(ref); open(WORK + '/qry.fa', 'w').write(qry)
        args = []
        r = rnd.random()
        pfx = rnd.choice(['out', 'pre', 'sub/pfx', 'a.b', 'nodir/x']) if r < .6 else None
        if pfx: args += ['-p', pfx] if rnd.random() < .7 else ['--prefix', pfx]
        if rnd.random() < .5: args += ['-l', str(rnd.choice([10, 15, 20]))]
        if rnd.random() < .5: args += ['-c', str(rnd.choice([20, 40]))]
        if rnd.random() < .2: args.append('--nodelta')
        if rnd.random() < .3: args.append('-o' if rnd.random() < .5 else '--coords')
        if rnd.random() < .15: args.append(rnd.choice(['--mum', '--maxmatch', '--mumreference']))
        if rnd.random() < .05: args += ['--mum', '--maxmatch']
        if rnd.random() < .1: args += ['-b', rnd.choice(['0', '20000', '100'])]
        if rnd.random() < .1: args.append('-f')
        p = subprocess.run(['nucmer'] + args + [WORK + '/ref.fa', WORK + '/qry.fa'], cwd=WORK, capture_output=True, text=True)
        files = {}
        for dp, dn, fn in os.walk(WORK):
            for f in fn:
                rel = os.path.relpath(os.path.join(dp, f), WORK)
                if rel not in ('ref.fa', 'qry.fa'): files[rel] = open(os.path.join(dp, f)).read()
        err = '\n'.join(l.replace(WORK + '/', '/W/') for l in p.stderr.split('\n') if not l.startswith(('# CONSTRUCTIONTIME', '# COMPLETETIME', '# SPACE')))
        print(json.dumps({'ref': ref, 'qry': qry, 'args': args, 'files': files, 'err': err, 'rc': p.returncode})); rows += 1
    sys.stderr.write('rows=%d\n' % rows)

main()
