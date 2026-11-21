#!/usr/bin/env python3
"""A3-9e oracle: `delta2maf` (Debian MUMmer 3.23 package) on nucmer delta files of multi-record, both-strand inputs: no ids / a reference id / a query id / both, -q -r, -w (also values the program rejects), -x,
ids without alignments (the program prints only the maf header and exits with 1).
Rows: {"ref","qry","delta","args","out","rc"}. Usage: python3 scripts/alignment/generate_delta2maf_oracle.py [n=140] [seed=20261204] | gzip > delta2maf_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09maf'
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 140
    rnd = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 20261204)
    def rseq(k): return ''.join(rnd.choice('acgt') for _ in range(k))
    def rc(s): return s.translate(COMP)[::-1]
    def mutate(s, rate):
        o = []
        for c in s:
            r = rnd.random()
            if r < rate * .6: o.append(rnd.choice('acgt'))
            elif r < rate * .75: pass
            elif r < rate: o.append(c + rnd.choice('acgt') * rnd.choice([1, 1, 3]))
            else: o.append(c)
        return ''.join(o)
    def fasta(recs, w):
        out = []
        for nm, sq in recs:
            out.append('>' + nm + '\n')
            for i in range(0, len(sq), w): out.append(sq[i:i + w] + '\n')
        return ''.join(out)
    rows = 0
    for case in range(n):
        blocks = [rseq(rnd.randint(120, 700)) for _ in range(rnd.randint(2, 5))]
        nref = rnd.choice([1, 2, 3]); nq = rnd.choice([1, 2, 3, 4])
        refs = [[] for _ in range(nref)]; qs = [[] for _ in range(nq)]
        for b in blocks:
            refs[rnd.randrange(nref)].append(mutate(b, .01)); refs[rnd.randrange(nref)].append(rseq(rnd.randint(0, 60)))
            for _ in range(rnd.choice([1, 1, 2])):
                piece = mutate(b, rnd.choice([0, .02, .05, .1]))
                if rnd.random() < .35: piece = rc(piece)
                qs[rnd.randrange(nq)].append(piece); qs[rnd.randrange(nq)].append(rseq(rnd.randint(0, 50)))
        rrec = [('ref%d' % (i + 1), ''.join(p) or rseq(30)) for i, p in enumerate(refs)]
        qrec = [('qry%d' % (i + 1), ''.join(p) or rseq(30)) for i, p in enumerate(qs)]
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
        ref = fasta(rrec, 60); qry = fasta(qrec, 70)
        open(WORK + '/ref.fa', 'w').write(ref); open(WORK + '/qry.fa', 'w').write(qry)
        p = subprocess.run(['nucmer', '-p', 'out', '-l', str(rnd.choice([12, 15, 20])), '-c', str(rnd.choice([20, 40, 65])), '--maxmatch' if rnd.random() < .5 else '--mumreference', 'ref.fa', 'qry.fa'], cwd=WORK, capture_output=True)
        if p.returncode != 0 or not os.path.exists(WORK + '/out.delta'): continue
        delta = open(WORK + '/out.delta').read()
        recs = [tuple(l[1:].split()[:2]) for l in delta.split('\n') if l.startswith('>')]
        if not recs: continue
        for _ in range(3):
            args = []
            if rnd.random() < .25: args.append(rnd.choice(['-q', '-r', '-q', '-r']))
            if rnd.random() < .1: args += ['-q', '-r']
            if rnd.random() < .2: args += ['-w', rnd.choice(['5', '12', '30', '100'])]
            if rnd.random() < .15: args += ['-x', rnd.choice(['1', '7'])]
            r = rnd.random(); ids = []
            if r < .3: pass
            elif r < .5: ids = [rnd.choice(recs)[0]]
            elif r < .75: ids = [rnd.choice(recs)[0], rnd.choice(recs)[1]] if rnd.random() < .5 else [rnd.choice(recs)[0], rnd.choice(recs)[1]]
            elif r < .85:
                ids = ['', rnd.choice(recs)[1]]          # the query id only: an empty reference id argument
            else: ids = ['nonexistent']
            pr = subprocess.run(['delta2maf'] + args + ['out.delta'] + ids, cwd=WORK, capture_output=True)
            print(json.dumps({'ref': ref, 'qry': qry, 'delta': delta, 'args': args + ['out.delta'] + ids, 'out': pr.stdout.decode('latin-1'), 'rc': pr.returncode})); rows += 1
    sys.stderr.write('rows=%d\n' % rows)

main()
