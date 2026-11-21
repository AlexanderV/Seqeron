#!/usr/bin/env python3
"""A3-9d oracle (downstream): the MUMmer 3.23 delta programs `delta-filter`, `show-coords`, `show-aligns`, `show-snps`, `show-diff` and `show-tiling` on PROMER (amino-acid) delta files written by `promer` from
related protein-coding DNA (genes back-translated with random codons, mutated at the amino-acid level, both strands, several records; genome-like query contigs cut from a mutated copy of the reference for show-tiling).
Every program is run with random option sets (LD_PRELOAD shim: time() = 1700000000, TZ=UTC). Same row format as generate_mummer_show_oracle.py: {"ref","qry","delta","runs":[{"prog","args","stdin","out","rc","files"}]}.
Usage: python3 scripts/alignment/generate_promer_show_oracle.py [n=170] [seed=20261202] [faketime.so] | gzip > promer_show_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09pro2'
COMP = str.maketrans('ACGTacgtRYKMBDHVrykmbdhv', 'TGCAtgcaYRMKVHDByrmkvhdb')


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 170
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261202
    shim = sys.argv[3] if len(sys.argv) > 3 else './faketime.so'
    rnd = random.Random(seed)
    env = dict(os.environ, LD_PRELOAD=os.path.abspath(shim), TZ='UTC')

    CODON = {}
    bases = 'TCAG'
    aas = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
    for i, a in enumerate(aas):
        CODON.setdefault(a, []).append(bases[i // 16] + bases[(i // 4) % 4] + bases[i % 4])
    AAS = 'ACDEFGHIKLMNPQRSTVWY'
    SIM = {'I': 'LVM', 'L': 'IVM', 'V': 'ILM', 'M': 'ILV', 'F': 'YW', 'Y': 'FW', 'W': 'FY', 'K': 'RQ', 'R': 'KQ', 'D': 'EN', 'E': 'DQ', 'S': 'TA', 'T': 'SA', 'A': 'ST', 'N': 'DQ', 'Q': 'EK', 'H': 'NQ'}

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

    def back(p): return ''.join(rnd.choice(CODON[a]) for a in p)

    def rc(s): return s.translate(COMP)[::-1]

    def spacer(k): return ''.join(rnd.choice('ACGT') for _ in range(k))

    def fasta(records, width):
        out = []
        for name, seq in records:
            out.append('>' + name + '\n')
            for i in range(0, len(seq), width): out.append(seq[i:i + width] + '\n')
        return ''.join(out)

    def decorate(seq):
        if rnd.random() < 0.1:
            i = rnd.randrange(len(seq) + 1); seq = seq[:i] + 'N' * rnd.choice([1, 3, 12]) + seq[i:]
        if rnd.random() < 0.08:
            i = rnd.randrange(len(seq) + 1); seq = seq[:i] + rnd.choice('RYKMSWBDHV') + seq[i:]
        if rnd.random() < 0.1: seq = seq.lower()
        return seq

    def make_case():
        genes = [prot(rnd.randint(30, 200)) for _ in range(rnd.randint(2, 6))]
        nref = rnd.choice([1, 1, 2, 3]); nq = rnd.choice([1, 2, 3, 4])
        refs = [[] for _ in range(nref)]; qs = [[] for _ in range(nq)]
        for g in genes:
            gd = back(g)
            refs[rnd.randrange(nref)].append(spacer(rnd.randint(0, 80))); refs[rnd.randrange(nref)].append(rc(gd) if rnd.random() < .25 else gd)
            for _ in range(rnd.choice([1, 1, 1, 2])):
                d = back(mutate_prot(g, rnd.choice([0, .05, .15, .3])))
                if rnd.random() < .12 and len(d) > 30:
                    i = rnd.randrange(3, len(d) - 3); d = d[:i] + spacer(rnd.choice([1, 2, 4])) + d[i:]
                if rnd.random() < .3: d = rc(d)
                if rnd.random() < .3: d = d[rnd.randint(0, 20):len(d) - rnd.randint(0, 20)]
                qs[rnd.randrange(nq)].append(spacer(rnd.randint(0, 80))); qs[rnd.randrange(nq)].append(d)
        rrec = [('ref%d' % (i + 1), decorate(''.join(p) or spacer(50))) for i, p in enumerate(refs)]
        qrec = [('qry%d' % (i + 1), decorate(''.join(p) or spacer(50))) for i, p in enumerate(qs)]
        return rrec, qrec

    def make_tiling_case():
        nref = rnd.choice([1, 1, 2])
        refprots = [[prot(rnd.randint(60, 220)) for _ in range(rnd.randint(4, 9))] for _ in range(nref)]
        refs = []; mut = []
        for ps in refprots:
            r = ''.join(spacer(rnd.randint(0, 60)) + back(p) for p in ps); refs.append(r)
            mut.append(''.join(spacer(rnd.randint(0, 60)) + back(mutate_prot(p, rnd.choice([0, .03, .1]))) for p in ps))
        qs = []
        for m in mut:
            pos = rnd.randint(0, 100)
            while pos < len(m) - 200:
                ln = rnd.randint(300, 900)
                piece = m[pos:pos + ln]
                if rnd.random() < .3: piece = rc(piece)
                qs.append(piece)
                pos += max(120, ln - rnd.randint(-100, 300))
        for _ in range(rnd.choice([0, 1, 2])): qs.append(spacer(rnd.randint(150, 500)))
        rnd.shuffle(qs)
        return [('ref%d' % (i + 1), x) for i, x in enumerate(refs)], [('c%d' % (i + 1), x) for i, x in enumerate(qs)]

    def run(cmd, stdin=None):
        p = subprocess.run(cmd, cwd=WORK, env=env, input=stdin, capture_output=True, text=True, timeout=120)
        return p

    rows = 0; runs_total = {}
    for case in range(n):
        tiling = rnd.random() < 0.4
        rrec, qrec = make_tiling_case() if tiling else make_case()
        ref = fasta(rrec, 60); qry = fasta(qrec, 70)
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
        open(WORK + '/ref.fa', 'w').write(ref); open(WORK + '/qry.fa', 'w').write(qry)
        nuc = ['promer', '-p', 'out', '-l', str(rnd.choice([4, 5, 6, 8])), '-c', str(rnd.choice([10, 20, 30]))]
        nuc += [rnd.choice(['--maxmatch', '--mum', '--mumreference', '--maxmatch'])]
        if rnd.random() < 0.3: nuc += ['-x', str(rnd.choice([1, 3]))]
        p = subprocess.run(nuc + ['ref.fa', 'qry.fa'], cwd=WORK, capture_output=True, text=True, timeout=300)
        if p.returncode != 0 or not os.path.exists(WORK + '/out.delta'): continue
        delta = open(WORK + '/out.delta').read()
        pairs = []
        cur = None
        for line in delta.split('\n'):
            if line.startswith('>'):
                f = line[1:].split(); pairs.append((f[0], f[1]))
        if not pairs: continue
        runs = []

        def record(prog, args, stdin=None):
            p = run([prog] + args + ['out.delta'] + (extra if False else []), stdin)
            runs.append({'prog': prog, 'args': args, 'stdin': stdin, 'out': p.stdout, 'rc': p.returncode})
            runs_total[prog] = runs_total.get(prog, 0) + 1
            return p

        # delta-filter
        for _ in range(2):
            args = []
            for opt in rnd.sample(['-q', '-r', '-g', '-1', '-m'], rnd.choice([0, 1, 1, 2])): args.append(opt)
            if rnd.random() < 0.3: args += ['-i', rnd.choice(['80', '90', '95', '99.5'])]
            if rnd.random() < 0.3: args += ['-l', rnd.choice(['20', '50', '100', '200'])]
            if rnd.random() < 0.2: args += ['-u', rnd.choice(['50', '90', '100'])]
            if rnd.random() < 0.2: args += ['-o', rnd.choice(['0', '10', '50', '100'])]
            if rnd.random() < 0.15: args += ['-e', rnd.choice(['0', '5', '20'])]
            record('delta-filter', args)
        # show-coords
        for _ in range(3):
            args = rnd.sample(['-b', '-c', '-d', '-H', '-k', '-l', '-o', '-T', '-k'], rnd.choice([0, 1, 2, 3]))
            if rnd.random() < 0.25: args += ['-I', rnd.choice(['80', '95', '99'])]
            if rnd.random() < 0.25: args += ['-L', rnd.choice(['30', '100'])]
            r = rnd.random()
            if r < 0.25: args += ['-r']
            elif r < 0.4: args += ['-q']
            if rnd.random() < 0.08: args += ['-g', '-r'] if '-r' not in args and '-q' not in args else ['-g']
            if rnd.random() < 0.04: args += ['-B']
            if rnd.random() < 0.04: args += ['-w']
            record('show-coords', args)
        # show-aligns
        for _ in range(2):
            rid, qid = rnd.choice(pairs)
            args = []
            if rnd.random() < 0.3: args += [rnd.choice(['-q', '-r'])]
            if rnd.random() < 0.4: args += ['-w', rnd.choice(['12', '30', '60', '100'])]
            if rnd.random() < 0.4: args += ['-x', rnd.choice(['1', '2', '3', '3', '7'])]
            p = run(['show-aligns'] + args + ['out.delta', rid, qid])
            runs.append({'prog': 'show-aligns', 'args': args + [rid, qid], 'stdin': None, 'out': p.stdout, 'rc': p.returncode})
            runs_total['show-aligns'] = runs_total.get('show-aligns', 0) + 1
        # show-snps
        for _ in range(3):
            args = rnd.sample(['-C', '-H', '-I', '-l', '-T'], rnd.choice([0, 1, 2]))
            r = rnd.random()
            if r < 0.3: args += ['-r']
            elif r < 0.6: args += ['-q']
            if rnd.random() < 0.5: args += ['-x', rnd.choice(['1', '3', '8'])]
            stdin = None
            if rnd.random() < 0.12:
                sc = run(['show-coords', '-H'] + rnd.choice([[], ['-r'], ['-q']]) + ['out.delta']).stdout
                lines = [l for l in sc.split('\n') if l.strip()]
                if lines:
                    sel = rnd.sample(lines, rnd.randint(1, min(3, len(lines))))
                    stdin = '\n'.join(sel) + '\n'; args += ['-S']
            record('show-snps', args, stdin)
        # show-diff
        for _ in range(2):
            args = rnd.sample(['-f', '-H'], rnd.choice([0, 1, 2])) + [rnd.choice(['-r', '-q', '-r'])]
            record('show-diff', args)
        # show-tiling (the nucmer data the program is made for: contigs of one genome against a reference)
        for _ in range(4 if tiling else 1):
            args = []
            if rnd.random() < 0.2: args += ['-a']
            elif rnd.random() < 0.15: args += ['-x']
            if rnd.random() < 0.2: args += ['-c']
            if rnd.random() < 0.3: args += ['-g', rnd.choice(['-1', '0', '100', '1000', '5000'])]
            if rnd.random() < 0.3: args += ['-i', rnd.choice(['0', '80', '90', '99'])]
            if rnd.random() < 0.3: args += ['-l', rnd.choice(['1', '100', '400', '800'])]
            if rnd.random() < 0.3: args += ['-v', rnd.choice(['0', '50', '80', '95', '99'])]
            if rnd.random() < 0.3: args += ['-V', rnd.choice(['0', '5', '10', '30'])]
            if rnd.random() < 0.15: args += ['-R']
            files = {}
            for opt, name in (('-t', 'tile.t'), ('-u', 'tile.u'), ('-p', 'tile.p')):
                if rnd.random() < 0.25:
                    args += [opt, name]; files[opt] = name
                    try: os.remove(WORK + '/' + name)
                    except OSError: pass
            p = run(['show-tiling'] + args + ['out.delta'])
            got = {opt: (open(WORK + '/' + nm).read() if os.path.exists(WORK + '/' + nm) else None) for opt, nm in files.items()}
            runs.append({'prog': 'show-tiling', 'args': args, 'stdin': None, 'out': p.stdout, 'rc': p.returncode, 'files': got})
            runs_total['show-tiling'] = runs_total.get('show-tiling', 0) + 1
        print(json.dumps({'ref': ref, 'qry': qry, 'delta': delta, 'runs': runs}))
        rows += 1
    sys.stderr.write('cases=%d runs=%s\n' % (rows, runs_total))


main()
