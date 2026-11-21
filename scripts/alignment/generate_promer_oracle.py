#!/usr/bin/env python3
"""A3-9d oracle: `promer` (MUMmer 3.23 perl script 3.07 + prepro / mummer / mgaps / postpro) on related protein-coding DNA (random proteins back-translated with random codons, the copies mutated at the amino-acid level:
substitutions, in-frame indels, frame shifts, both strands, several records, N / IUPAC characters, lower case, tiny and unrelated records). Every run of promer is recorded with every option family
(--mum --mumreference --maxmatch -b -c -d -g -l -m -x --noextend --nooptimize --nodelta -o -p) together with the intermediate files computed by the sub-programs run by hand (prepro -r / -q, mummer | mgaps) so
that every stage can be checked on its own. Output: JSON lines {"ref","qry","args","delta","cluster","coords","err","rc","aaref","aaqry","mgaps"} (stderr without the timing lines, work paths replaced).
Usage: python3 scripts/alignment/generate_promer_oracle.py [n=420] [seed=20261201] | gzip > promer_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09pro'
CODON = {}
bases = 'TCAG'
aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
for i, a in enumerate(aas):
    CODON.setdefault(a, []).append(bases[i // 16] + bases[(i // 4) % 4] + bases[i % 4])
AAS = 'ACDEFGHIKLMNPQRSTVWY'
SIM = {'I': 'LVM', 'L': 'IVM', 'V': 'ILM', 'M': 'ILV', 'F': 'YW', 'Y': 'FW', 'W': 'FY', 'K': 'RQ', 'R': 'KQ', 'D': 'EN', 'E': 'DQ', 'S': 'TA', 'T': 'SA', 'A': 'ST', 'N': 'DQ', 'Q': 'EK', 'H': 'NQ'}
COMP = str.maketrans('ACGTacgtRYKMBDHVrykmbdhv', 'TGCAtgcaYRMKVHDByrmkvhdb')

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 420
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261201
    rnd = random.Random(seed)

    def prot(k): return ''.join(rnd.choice(AAS) for _ in range(k))

    def mutate_prot(p, rate):
        o = []
        for a in p:
            r = rnd.random()
            if r < rate * .6: o.append(rnd.choice(SIM.get(a, AAS)) if rnd.random() < .6 else rnd.choice(AAS))
            elif r < rate * .75: pass
            elif r < rate * .9: o.append(a); o.append(prot(rnd.choice([1, 1, 2, 4])))
            else: o.append(a)
        return ''.join(o)

    def back(p, stops=False):
        return ''.join(rnd.choice(CODON[a]) for a in p)

    def rcomp(s): return s.translate(COMP)[::-1]

    def spacer(k): return ''.join(rnd.choice('ACGT') for _ in range(k))

    def decorate(s):
        if rnd.random() < .12:
            i = rnd.randrange(len(s) + 1); s = s[:i] + 'N' * rnd.choice([1, 2, 3, 6, 12]) + s[i:]
        if rnd.random() < .1:
            i = rnd.randrange(len(s) + 1); s = s[:i] + rnd.choice('RYKMSWBDHV') + s[i:]
        if rnd.random() < .06:
            i = rnd.randrange(len(s) + 1); s = s[:i] + rnd.choice('efzq*-') + s[i:]      # characters promer cannot translate / read
        if rnd.random() < .1: s = s.lower()
        return s

    def fasta(recs, width):
        out = []
        for nm, sq in recs:
            out.append('>' + nm + '\n')
            if width <= 0: out.append(sq + '\n')
            else:
                for i in range(0, len(sq), width): out.append(sq[i:i + width] + '\n')
        return ''.join(out)

    def make_case():
        genes = [prot(rnd.randint(25, 220)) for _ in range(rnd.randint(1, 5))]
        nref = rnd.choice([1, 1, 2, 3]); nq = rnd.choice([1, 1, 2, 3])
        refs = [[] for _ in range(nref)]; qs = [[] for _ in range(nq)]
        for g in genes:
            t = refs[rnd.randrange(nref)]
            t.append(spacer(rnd.randint(0, 60)))
            gd = back(g)
            t.append(rcomp(gd) if rnd.random() < .25 else gd)
            for _ in range(rnd.choice([1, 1, 1, 2])):
                p = mutate_prot(g, rnd.choice([0, .05, .15, .3]))
                d = back(p)
                if rnd.random() < .12 and len(d) > 30:
                    i = rnd.randrange(3, len(d) - 3); d = d[:i] + spacer(rnd.choice([1, 2, 4])) + d[i:]       # frame shift
                if rnd.random() < .3: d = rcomp(d)
                if rnd.random() < .3: d = d[rnd.randint(0, 20):len(d) - rnd.randint(0, 20)]
                q = qs[rnd.randrange(nq)]
                q.append(spacer(rnd.randint(0, 60))); q.append(d)
        rrec = [('ref%d desc' % (i + 1), decorate(''.join(p) or spacer(rnd.randint(0, 20)))) for i, p in enumerate(refs)]
        qrec = [('qry%d' % (i + 1), decorate(''.join(p) or spacer(rnd.randint(0, 20)))) for i, p in enumerate(qs)]
        if rnd.random() < .08: rrec.append(('tiny', spacer(rnd.randint(0, 5))))
        if rnd.random() < .08: qrec.append(('tinyq', spacer(rnd.randint(0, 5))))
        return rrec, qrec

    def sh(cmd, **kw):
        return subprocess.run(cmd, cwd=WORK, capture_output=True, text=True, timeout=300, **kw)

    def clean_err(e):
        out = []
        for l in e.split('\n'):
            if l.startswith(('# CONSTRUCTIONTIME', '# COMPLETETIME', '# SPACE')): continue
            out.append(l.replace(WORK + '/', ''))
        return '\n'.join(out)

    rows = 0; stats = {}
    for case in range(n):
        rrec, qrec = make_case()
        ref = fasta(rrec, rnd.choice([0, 60, 70])); qry = fasta(qrec, rnd.choice([0, 60, 70]))
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
        open(WORK + '/ref.fa', 'w').write(ref); open(WORK + '/qry.fa', 'w').write(qry)
        args = []
        r = rnd.random()
        if r < .2: args.append('--mum')
        elif r < .35: args.append('--maxmatch')
        elif r < .45: args.append(rnd.choice(['--mumreference', '--mumcand']))
        if rnd.random() < .5: args += ['-l', str(rnd.choice([4, 5, 6, 7, 8, 10]))]
        if rnd.random() < .4: args += ['-c', str(rnd.choice([5, 10, 20, 30, 50]))]
        if rnd.random() < .35: args += ['-g', str(rnd.choice([5, 15, 30, 60, 100]))]
        if rnd.random() < .3: args += ['-b', str(rnd.choice([10, 30, 60, 120, 200, 20000, 0]))]
        if rnd.random() < .3: args += ['-d', str(rnd.choice(['0.05', '.11', '0.3', '0']))]
        if rnd.random() < .3: args += ['-m', str(rnd.choice([0, 3, 8, 15, 30, -2]))]
        if rnd.random() < .35: args += ['-x', str(rnd.choice([1, 2, 3, 3, 1, 0, 9]))]
        if rnd.random() < .12: args.append('--noextend')
        if rnd.random() < .15: args.append('--nooptimize')
        if rnd.random() < .08: args.append('--nodelta')
        if rnd.random() < .25: args.append('-o')
        if rnd.random() < .1: args += ['-p', 'pre']
        pfx = args[args.index('-p') + 1] if '-p' in args else 'out'
        p = sh(['promer'] + args + [WORK + '/ref.fa', WORK + '/qry.fa'])
        files = {}
        for ext in ('delta', 'cluster', 'coords', 'aaref', 'aaqry', 'mgaps'):
            fp = '%s/%s.%s' % (WORK, pfx, ext)
            files[ext] = open(fp).read() if os.path.exists(fp) else None
        # the intermediate files by hand, with the settings the script derived from the options
        def opt(name, default):
            return args[args.index(name) + 1] if name in args else str(default)
        mask = opt('-m', 8); size = opt('-l', 6); clus = opt('-c', 20); gap = opt('-g', 30); diff = opt('-d', '.11')
        algo = '-mum' if '--mum' in args else '-maxmatch' if '--maxmatch' in args else '-mumreference'
        pr = sh(['/usr/lib/mummer/prepro', '-m', mask, '-r', WORK + '/ref.fa']); aaref = pr.stdout
        pq = sh(['/usr/lib/mummer/prepro', '-m', mask, '-q', WORK + '/qry.fa']); aaqry = pq.stdout
        mg = None
        if pr.returncode == 0 and pq.returncode == 0:
            open(WORK + '/h.aaref', 'w').write(aaref); open(WORK + '/h.aaqry', 'w').write(aaqry)
            mm = subprocess.run('mummer %s -l %s h.aaref h.aaqry 2>/dev/null | mgaps -l %s -s %s -f %s' % (algo, size, clus, gap, diff), shell=True, cwd=WORK, capture_output=True, text=True)
            if mm.returncode == 0: mg = mm.stdout
        row = {'ref': ref, 'qry': qry, 'args': args, 'rc': p.returncode, 'err': clean_err(p.stderr), 'delta': files['delta'], 'cluster': files['cluster'], 'coords': files['coords'],
               'aaref': aaref, 'aaqry': aaqry, 'mgaps': mg, 'prepro_err': (pr.stderr, pq.stderr)}
        stats[p.returncode] = stats.get(p.returncode, 0) + 1
        print(json.dumps(row))
        rows += 1
    sys.stderr.write('rows=%d rc=%s\n' % (rows, stats))

main()
