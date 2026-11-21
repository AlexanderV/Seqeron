#!/usr/bin/env python3
"""A3-9e oracle: `gaps`, `mummer-annotate` (annotate), `delta2blocks` and the script `run-mummer1` of MUMmer 3.23 (Debian 3.23+dfsg-8).
Rows (JSON lines):
  {"kind":"gaps","file","rflag","input","out","rc"}                                   gaps <file> [-r] < input (random and mummer-made match lists, ties in the second start, overlaps; the empty list = SIGSEGV)
  {"kind":"delta2blocks","ref","qry","delta","args","out","rc"}                         delta2blocks on nucmer deltas (and promer deltas: SIGSEGV)
  {"kind":"annotate","files":{name:text},"gaps","data","out","errors","ub"}              mummer-annotate gaps data; "ub": the output changes when malloc fills fresh memory with junk (MALLOC_PERTURB_ = 165): the program read
                                                                                         memory it never wrote (a window before / past a sequence end) - the genuine non-deterministic class, not compared
  {"kind":"runmummer1","ref","qry","prefix","args","files":{ext:text},"rc"}             run-mummer1 ref qry prefix [-r]
Usage: python3 scripts/alignment/generate_mummer_gaps_oracle.py [seed=20261203] | gzip > mummer_gaps_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09gaps'
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')

def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261203
    rnd = random.Random(seed)
    def rseq(k): return ''.join(rnd.choice('acgt') for _ in range(k))
    def rc(s): return s.translate(COMP)[::-1]
    def mutate(s, rate):
        o = []
        for c in s:
            r = rnd.random()
            if r < rate * .6: o.append(rnd.choice('acgt'))
            elif r < rate * .75: pass
            elif r < rate: o.append(c + rnd.choice('acgt'))
            else: o.append(c)
        return ''.join(o)
    def fasta(recs, width=60):
        out = []
        for nm, sq in recs:
            out.append('>' + nm + '\n')
            for i in range(0, len(sq), width): out.append(sq[i:i + width] + '\n')
        return ''.join(out)
    def run(cmd, **kw): return subprocess.run(cmd, cwd=WORK, capture_output=True, **kw)
    def reset():
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)

    rows = []
    # ---- gaps
    for k in range(260):
        r = rnd.random()
        if k == 0: inp = ''
        else:
            n = rnd.randint(1, 14)
            ms = []
            if r < .5:      # a diagonal band of matches with errors between them
                p1, p2 = rnd.randint(1, 200), rnd.randint(1, 200)
                for _ in range(n):
                    l = rnd.randint(5, 60); ms.append((p1, p2, l)); p1 += l + rnd.randint(-8, 40); p2 += l + rnd.randint(-8, 40)
                if rnd.random() < .5:
                    for _ in range(rnd.randint(1, 4)): ms.append((rnd.randint(1, 600), rnd.randint(1, 600), rnd.randint(5, 50)))
                rnd.shuffle(ms) if rnd.random() < .5 else None
            elif r < .8:    # ties in the second start, overlaps
                for _ in range(n):
                    ms.append((rnd.randint(1, 80), rnd.choice([10, 10, 25, 40, 60, rnd.randint(1, 100)]), rnd.randint(5, 40)))
            else:
                for _ in range(n): ms.append((rnd.randint(1, 1000), rnd.randint(1, 1000), rnd.randint(1, 100)))
            inp = ''.join('%8d  %8d  %8d\n' % m for m in ms)
            if rnd.random() < .08: inp = inp.rstrip('\n') + ' 7\n'          # a partial triple at the end
        rflag = rnd.random() < .4
        reset()
        p = run(['gaps', 'g1.fa'] + (['-r'] if rflag else []), input=inp.encode())
        rows.append({'kind': 'gaps', 'file': 'g1.fa', 'rflag': rflag, 'input': inp, 'out': p.stdout.decode('latin-1'), 'rc': p.returncode})

    # ---- delta2blocks (nucmer deltas; promer deltas crash)
    for k in range(125):
        base = rseq(rnd.randint(300, 1500))
        ref = [('ref%d' % i, mutate(base[rnd.randint(0, 50):], .01) + rseq(rnd.randint(0, 80))) for i in range(rnd.choice([1, 1, 2]))]
        qry = []
        for i in range(rnd.choice([1, 2, 3])):
            piece = mutate(base[rnd.randint(0, 100):rnd.randint(len(base) // 2, len(base))], rnd.choice([0, .02, .06, .12]))
            if rnd.random() < .35: piece = rc(piece)
            qry.append(('q%d' % i, rseq(rnd.randint(0, 30)) + piece + rseq(rnd.randint(0, 30))))
        reset()
        open(WORK + '/ref.fa', 'w').write(fasta(ref)); open(WORK + '/qry.fa', 'w').write(fasta(qry, 70))
        promer = k % 12 == 11
        cmd = ['promer', '-p', 'out', '-l', '5', '-c', '10'] if promer else ['nucmer', '-p', 'out', '-l', str(rnd.choice([12, 15, 20])), '-c', str(rnd.choice([20, 40, 65]))]
        if run(cmd + ['ref.fa', 'qry.fa']).returncode != 0 or not os.path.exists(WORK + '/out.delta'): continue
        delta = open(WORK + '/out.delta').read()
        pairs = [tuple(l[1:].split()[:2]) for l in delta.split('\n') if l.startswith('>')]
        if not pairs: continue
        for _ in range(2):
            rid, qid = rnd.choice(pairs)
            args = []
            if rnd.random() < .3: args.append(rnd.choice(['-q', '-r']))
            p = run(['delta2blocks'] + args + ['out.delta', rid, qid])
            rows.append({'kind': 'delta2blocks', 'ref': open(WORK + '/ref.fa').read(), 'qry': open(WORK + '/qry.fa').read(), 'delta': delta, 'args': args + [rid, qid], 'out': p.stdout.decode('latin-1'), 'rc': p.returncode})

    # ---- mummer-annotate and run-mummer1
    for k in range(345):
        base = rseq(rnd.randint(500, 2500))
        ref = base
        qpiece = mutate(base[rnd.randint(0, 150):len(base) - rnd.randint(0, 150)], rnd.choice([.005, .01, .02, .04]))
        if rnd.random() < .3: qpiece = qpiece[:len(qpiece) // 2] + rseq(rnd.randint(1, 40)) + qpiece[len(qpiece) // 2:]
        rev = rnd.random() < .5
        qry = rc(qpiece) if rev else qpiece
        if rnd.random() < .4: qry = rseq(rnd.randint(0, 120)) + qry + rseq(rnd.randint(0, 120))     # windows near the ends: sometimes undefined
        reset()
        open(WORK + '/ref.fa', 'w').write('>ref\n' + '\n'.join(ref[i:i + 60] for i in range(0, len(ref), 60)) + '\n')
        open(WORK + '/qry.fa', 'w').write('>qry\n' + '\n'.join(qry[i:i + 60] for i in range(0, len(qry), 60)) + '\n')
        if k % 3 == 2:
            args = ['-r'] if rev or rnd.random() < .5 else []
            p = run(['run-mummer1', 'ref.fa', 'qry.fa', 'pf'] + args)
            files = {}
            for ext in ('out', 'gaps', 'align', 'errorsgaps'):
                fp = '%s/pf.%s' % (WORK, ext)
                files[ext] = open(fp, encoding='latin-1').read() if os.path.exists(fp) else None
            if p.returncode == 0:      # the second run perturbs fresh malloc memory: a different result = a read of memory the program never wrote
                env = dict(os.environ, MALLOC_PERTURB_='165')
                run(['run-mummer1', 'ref.fa', 'qry.fa', 'pg'] + args, env=env)
                ub = any((open('%s/pg.%s' % (WORK, e), encoding='latin-1').read() if os.path.exists('%s/pg.%s' % (WORK, e)) else None) != files[e] for e in ('out', 'gaps', 'align', 'errorsgaps'))
            else: ub = False
            rows.append({'kind': 'runmummer1', 'ref': open(WORK + '/ref.fa').read(), 'qry': open(WORK + '/qry.fa').read(), 'prefix': 'pf', 'args': args, 'files': files, 'rc': p.returncode, 'ub': ub})
            continue
        mm = subprocess.run('mummer -mum -l %d %s ref.fa qry.fa 2>/dev/null | tail +2 | gaps ref.fa %s' % (rnd.choice([12, 15, 20]), '-r' if rev else '', '-r' if rev else ''), shell=True, cwd=WORK, capture_output=True)
        gaps = mm.stdout.decode('latin-1')
        if not gaps.strip(): continue
        open(WORK + '/g.gaps', 'w').write(gaps)
        def annotate(env=None):
            if os.path.exists(WORK + '/witherrors.gaps'): os.remove(WORK + '/witherrors.gaps')
            p = run(['mummer-annotate', 'g.gaps', 'qry.fa'], env=env)
            return p, (open(WORK + '/witherrors.gaps', encoding='latin-1').read() if os.path.exists(WORK + '/witherrors.gaps') else None)
        p, we = annotate()
        if p.returncode != 0: continue
        p2, we2 = annotate(dict(os.environ, MALLOC_PERTURB_='165'))
        ub = p2.stdout != p.stdout or we2 != we
        rows.append({'kind': 'annotate', 'files': {'ref.fa': open(WORK + '/ref.fa').read()}, 'gaps': gaps, 'data': open(WORK + '/qry.fa').read(), 'out': p.stdout.decode('latin-1'), 'errors': we, 'ub': ub})
    # the genuine non-deterministic class: a gap line whose window reaches before the first base (the program reads the heap in front of its buffer)
    for k in range(14):
        base = rseq(rnd.randint(300, 600))
        reset()
        open(WORK + '/ref.fa', 'w').write('>ref\n' + base + '\n'); open(WORK + '/qry.fa', 'w').write('>qry\n' + base + '\n')
        g = '> ref.fa  Consistent matches\n%8d %8d %6d    none      -      -\n%8d %8d %6d    none      %d      %d\n' % (1, 1, 30, rnd.randint(1, 12), rnd.randint(1, 12), 30, rnd.randint(1, 6), rnd.randint(1, 6))
        open(WORK + '/g.gaps', 'w').write(g)
        def annotate(env=None):
            if os.path.exists(WORK + '/witherrors.gaps'): os.remove(WORK + '/witherrors.gaps')
            p = run(['mummer-annotate', 'g.gaps', 'qry.fa'], env=env)
            return p, (open(WORK + '/witherrors.gaps', encoding='latin-1').read() if os.path.exists(WORK + '/witherrors.gaps') else None)
        p, we = annotate(); p2, we2 = annotate(dict(os.environ, MALLOC_PERTURB_='165'))
        if p.returncode != 0: continue
        vg = run(['valgrind', '-q', '--error-exitcode=9', 'mummer-annotate', 'g.gaps', 'qry.fa'])      # proof: an invalid read in front of the buffer
        rows.append({'kind': 'annotate', 'files': {'ref.fa': open(WORK + '/ref.fa').read()}, 'gaps': g, 'data': open(WORK + '/qry.fa').read(), 'out': p.stdout.decode('latin-1'), 'errors': we,
                     'ub': p2.stdout != p.stdout or we2 != we or vg.returncode == 9, 'valgrind': vg.returncode == 9, 'perturb': p2.stdout != p.stdout or we2 != we})
    for r in rows: print(json.dumps(r))
    kinds = {}
    for r in rows: kinds[r['kind']] = kinds.get(r['kind'], 0) + 1
    sys.stderr.write('rows=%d %s ub=%d\n' % (len(rows), kinds, sum(1 for r in rows if r.get('ub'))))

main()
