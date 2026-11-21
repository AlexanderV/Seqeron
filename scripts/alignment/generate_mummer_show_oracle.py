#!/usr/bin/env python3
"""A3-9b oracle (stage 2): the MUMmer 3.23 delta programs `delta-filter`, `show-coords`, `show-aligns`, `show-snps`, `show-diff` (and `show-tiling`) on delta files written by `nucmer` from rearranged,
duplicated, inverted, mutated multi-contig inputs. Every program is run (LD_PRELOAD shim that makes time() = 1700000000 and TZ=UTC: the programs call srand (time (NULL)) and print the date of a btab line) with
random option combinations; the work directory is the fixed /tmp/b09nuc so the header paths are reproducible.
Output: JSON lines {"ref", "qry", "delta", "runs": [{"prog", "args", "stdin", "out", "rc"}]}.
Usage: python3 scripts/alignment/generate_mummer_show_oracle.py [n=150] [seed=20261122] [faketime.so] | gzip > mummer_show_oracle.jsonl.gz
(build the shim: printf '#include <time.h>\\ntime_t time(time_t *t){ if (t) *t = 1700000000; return 1700000000; }\\n' > ft.c; gcc -shared -fPIC -o faketime.so ft.c)
"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09nuc'
COMP = str.maketrans('ACGTacgtRYKMBDHVrykmbdhv', 'TGCAtgcaYRMKVHDByrmkvhdb')


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261122
    shim = sys.argv[3] if len(sys.argv) > 3 else './faketime.so'
    rnd = random.Random(seed)
    env = dict(os.environ, LD_PRELOAD=os.path.abspath(shim), TZ='UTC')

    def rseq(k, gc=0.5):
        return ''.join(rnd.choice('GC') if rnd.random() < gc else rnd.choice('AT') for _ in range(k))

    def mutate(s, rate):
        out = []
        for c in s:
            r = rnd.random()
            if r < rate * 0.6: out.append(rnd.choice('ACGT'))
            elif r < rate * 0.75: pass
            elif r < rate * 0.9: out.append(c); out.append(rnd.choice('ACGT') * rnd.choice([1, 1, 2, 4]))
            else: out.append(c)
        return ''.join(out)

    def rc(s):
        return s.translate(COMP)[::-1]

    def fasta(records, width):
        out = []
        for name, seq in records:
            out.append('>' + name + '\n')
            for i in range(0, len(seq), width): out.append(seq[i:i + width] + '\n')
        return ''.join(out)

    def decorate(seq):
        if rnd.random() < 0.15:
            i = rnd.randrange(len(seq) + 1); seq = seq[:i] + 'N' * rnd.choice([1, 3, 12]) + seq[i:]
        if rnd.random() < 0.1:
            i = rnd.randrange(len(seq) + 1); seq = seq[:i] + rnd.choice('RYKMSWBDHV') + seq[i:]
        if rnd.random() < 0.1:
            seq = seq.lower()
        return seq

    def make_case():
        blocks = [rseq(rnd.randint(80, 600), rnd.choice([0.35, 0.5, 0.6])) for _ in range(rnd.randint(3, 8))]
        nref = rnd.choice([1, 1, 2, 3]); nq = rnd.choice([1, 2, 3, 4, 5])
        refs = [[] for _ in range(nref)]
        for b in blocks:
            refs[rnd.randrange(nref)].append(mutate(b, 0.01)); refs[rnd.randrange(nref)].append(rseq(rnd.randint(0, 120)))
        qs = [[] for _ in range(nq)]
        for b in blocks:
            for _ in range(rnd.choice([1, 1, 1, 2])):               # duplicated blocks give overlapping alignments
                piece = mutate(b, rnd.choice([0.0, 0.02, 0.05, 0.1]))
                if rnd.random() < 0.3: piece = rc(piece)
                if rnd.random() < 0.3 and len(piece) > 120:
                    a = rnd.randint(30, len(piece) - 30); b2 = a + rnd.randint(1, 60); piece = piece[:a] + rseq(rnd.randint(0, 40)) + piece[b2:]
                if rnd.random() < 0.3: piece = piece[rnd.randint(0, 40):len(piece) - rnd.randint(0, 40)]
                qs[rnd.randrange(nq)].append(piece); qs[rnd.randrange(nq)].append(rseq(rnd.randint(0, 100)))
        rrec = [('ref%d' % (i + 1), decorate(''.join(p) or rseq(50))) for i, p in enumerate(refs)]
        qrec = [('qry%d' % (i + 1), decorate(''.join(p) or rseq(50))) for i, p in enumerate(qs)]
        rnd.shuffle(qrec) if rnd.random() < 0.3 else None
        return rrec, qrec

    def make_tiling_case():
        nref = rnd.choice([1, 1, 2])
        refs = [rseq(rnd.randint(1500, 4000), rnd.choice([0.4, 0.5, 0.6])) for _ in range(nref)]
        # a repeat copied into the reference, so that some contigs align twice
        if rnd.random() < 0.4:
            seg = refs[0][200:200 + rnd.randint(80, 300)]
            k = rnd.randrange(len(refs)); i = rnd.randrange(len(refs[k])); refs[k] = refs[k][:i] + seg + refs[k][i:]
        qs = []
        for i, r in enumerate(refs):
            pos = rnd.randint(0, 100)
            while pos < len(r) - 150:
                ln = rnd.randint(200, 900)
                piece = r[pos:pos + ln]
                piece = mutate(piece, rnd.choice([0.0, 0.005, 0.02, 0.06]))
                if rnd.random() < 0.3: piece = rc(piece)
                if rnd.random() < 0.1: piece = piece[:len(piece) // 2] + rseq(rnd.randint(5, 60)) + piece[len(piece) // 2:]
                qs.append(piece)
                pos += max(60, ln - rnd.randint(-100, 300))          # overlap or gap between contigs
        for _ in range(rnd.choice([0, 1, 2])): qs.append(rseq(rnd.randint(100, 500)))   # unalignable contig
        if rnd.random() < 0.3: qs.append(rseq(rnd.randint(40, 100)))                      # tiny contig
        rnd.shuffle(qs)
        return [('ref%d' % (i + 1), x) for i, x in enumerate(refs)], [('c%d' % (i + 1), decorate(x) if rnd.random() < 0.1 else x) for i, x in enumerate(qs)]

    def run(cmd, stdin=None):
        p = subprocess.run(cmd, cwd=WORK, env=env, input=stdin, capture_output=True, text=True, timeout=120)
        return p

    rows = 0; runs_total = {}
    for case in range(n):
        tiling = rnd.random() < 0.45
        rrec, qrec = make_tiling_case() if tiling else make_case()
        ref = fasta(rrec, 60); qry = fasta(qrec, 70)
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
        open(WORK + '/ref.fa', 'w').write(ref); open(WORK + '/qry.fa', 'w').write(qry)
        nuc = ['nucmer', '-p', 'out', '-l', str(rnd.choice([10, 12, 15, 20])), '-c', str(rnd.choice([20, 30, 65]))]
        nuc += [rnd.choice(['--maxmatch', '--mum', '--mumreference', '--maxmatch'])]
        if rnd.random() < 0.3: nuc += ['--nosimplify']
        p = subprocess.run(nuc + ['ref.fa', 'qry.fa'], cwd=WORK, capture_output=True, text=True, timeout=120)
        if p.returncode != 0 or not os.path.exists(WORK + '/out.delta'): continue
        delta = open(WORK + '/out.delta').read()
        pairs = []
        cur = None
        for line in delta.split('\n'):
            if line.startswith('>'):
                f = line[1:].split(); pairs.append((f[0], f[1]))
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
            args = rnd.sample(['-b', '-c', '-d', '-H', '-k', '-l', '-o', '-T'], rnd.choice([0, 1, 2, 3]))
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
        # dnadiff: the perl script (nucmer --maxmatch + delta-filter + show-* + the report); the order of the hash-driven lines of the report and of .unref / .unqry is random (perl hash randomisation)
        for k, dargs in enumerate([[], ['-d', 'out.delta']]):
            pfx = 'dd%d' % k
            for suffix in ('report', 'delta', '1delta', 'mdelta', '1coords', 'mcoords', 'snps', 'rdiff', 'qdiff', 'unref', 'unqry'):
                try: os.remove('%s/%s.%s' % (WORK, pfx, suffix))
                except OSError: pass
            p = subprocess.run(['dnadiff', '-p', pfx] + ([] if dargs else ['ref.fa', 'qry.fa']) + dargs[0:0] + (['-d', 'out.delta'] if dargs else []), cwd=WORK, env=env if False else None, capture_output=True, text=True, timeout=300)
            files = {sfx: (open('%s/%s.%s' % (WORK, pfx, sfx)).read() if os.path.exists('%s/%s.%s' % (WORK, pfx, sfx)) else None) for sfx in ('report', 'delta', '1delta', 'mdelta', '1coords', 'mcoords', 'snps', 'rdiff', 'qdiff', 'unref', 'unqry')}
            runs.append({'prog': 'dnadiff', 'args': ['-d'] if dargs else [], 'stdin': None, 'out': '', 'rc': p.returncode, 'files': files})
            runs_total['dnadiff'] = runs_total.get('dnadiff', 0) + 1
        print(json.dumps({'ref': ref, 'qry': qry, 'delta': delta, 'runs': runs}))
        rows += 1
    sys.stderr.write('cases=%d runs=%s\n' % (rows, runs_total))


main()
