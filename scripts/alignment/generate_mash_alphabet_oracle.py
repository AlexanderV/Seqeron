#!/usr/bin/env python3
"""A5-10b2 (B09 F214) oracle: the alphabet options of Mash 2.3 `-a` (amino acids, implies -n -k 9), `-z TEXT` (own alphabet, implies -n), `-Z` (case kept) for
  sketch  (the whole .msh file, genomes and read sets, -i, -r -m -c -g -b),
  dist / triangle (the stdout text, -t -v -d -E -i -C ...),
  screen  (query sketch made with the alphabet; protein queries translate the mixture in 6 frames: reads / contigs of several files, the 2^20-character chunking, -w -i -v).
JSON lines: {"id", "kind", "args": [the mash command line after the program name, files as given], "files": [{"name", "fasta"}], "qmsh": base64 of the query sketch of a screen row, "out": stdout, "msh": base64 .msh of a sketch row}.
Usage: python3 scripts/alignment/generate_mash_alphabet_oracle.py OUT_PREFIX   (writes OUT_PREFIX_sketch.jsonl.gz, _dist.jsonl.gz, _screen.jsonl.gz) [n_sketch=300] [n_dist=300] [n_screen=300] [seed]     (needs mash 2.3 on PATH)
"""
import base64, gzip, json, os, random, subprocess, sys, tempfile

AA = 'ACDEFGHIKLMNPQRSTVWY'
CODON = {}
_b = 'ACGT'
_code = 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF'
for i, x in enumerate(_b):
    for j, y in enumerate(_b):
        for k, z in enumerate(_b):
            CODON.setdefault(_code[16 * i + 4 * j + k], []).append(x + y + z)


def main(out, n_sketch, n_dist, n_screen, seed):
    rnd = random.Random(seed)

    def prot(m, extra=0.0):
        s = []
        for _ in range(m):
            r = rnd.random()
            s.append(rnd.choice('XBZ*') if r < extra * 0.5 else rnd.choice(AA).lower() if r < extra else rnd.choice(AA))
        return ''.join(s)

    def dna(m, extra=0.0):
        s = []
        for _ in range(m):
            r = rnd.random()
            s.append('N' if r < extra * 0.5 else rnd.choice('acgt') if r < extra else rnd.choice('ACGT'))
        return ''.join(s)

    def mutate_p(s, rate):
        return ''.join(rnd.choice(AA) if rnd.random() < rate else c for c in s)

    def back(p):                                   # a DNA coding for the protein p (the stop of the code for '*')
        return ''.join(rnd.choice(CODON.get(c.upper(), ['NNN'])) for c in p)

    def revcomp(s): return s.translate(str.maketrans('ACGTacgtN', 'TGCAtgcaN'))[::-1]

    ALPHABETS = ['ACGT', 'ACGTN', 'ACDEFGHIKLMNPQRSTVWY', 'AB', 'ACGU', 'abc', 'XYZ', 'ACGTacgt', 'AC', 'ACDEFGHIKLMNPQRSTVWYX', 'ACDEFGHIKLMNPQRSTVWY*', 'RYKMSWN', '01', 'A-C', 'AACCGGTT']

    def fasta(recs):
        return ''.join('>%s%s\n%s\n' % (i, (' ' + c) if c else '', s) for i, c, s in recs)

    def run(d, args, ok_codes=(0,)):
        p = subprocess.run(['mash'] + args, cwd=d, capture_output=True, text=True, timeout=300)
        return p

    def write_files(d, files):
        for f in files:
            with open(os.path.join(d, f['name']), 'w', newline='\n') as h:
                h.write(f['fasta'])

    def options():
        """a random choice of the alphabet options: (args, kind)"""
        r = rnd.random()
        a = []
        if r < 0.4:
            a.append('-a')
            if rnd.random() < 0.5: a += ['-k', str(rnd.choice([3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 20, 32]))]
        elif r < 0.8:
            a += ['-z', rnd.choice(ALPHABETS)]
            a += ['-k', str(rnd.choice([1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 17, 21, 32]))] if rnd.random() < 0.85 else []
        if rnd.random() < 0.35: a.append('-Z')
        if rnd.random() < 0.1: a.append('-n')
        return a

    def seq_for(alpha_args):
        """a pool of sequences that matches the alphabet chosen by the options"""
        if '-a' in alpha_args:
            return lambda m, e=0.03: prot(m, e)
        if '-z' in alpha_args:
            letters = alpha_args[alpha_args.index('-z') + 1]
            def f(m, e=0.03):
                pool = letters if rnd.random() < 0.9 else letters + 'N'
                return ''.join(rnd.choice(pool) if rnd.random() > e else rnd.choice('NacgtXx') for _ in range(m))
            return f
        return lambda m, e=0.03: dna(m, e)

    rows = []
    with tempfile.TemporaryDirectory() as d:
        # ------------------------------------------------------------------ sketch files
        while sum(1 for r in rows if r['kind'] == 'sketch') < n_sketch:
            al = options()
            gen = seq_for(al)
            nfiles = rnd.choice([1, 1, 2, 3])
            files = []
            base = gen(rnd.randint(200, 3000), 0.0)
            for i in range(nfiles):
                recs = []
                for j in range(rnd.choice([1, 1, 2, 3, 5])):
                    sq = ''.join(rnd.choice(base) if rnd.random() < 0.1 else c for c in base[rnd.randint(0, 50):rnd.randint(150, len(base))]) if rnd.random() < 0.7 else gen(rnd.randint(2, 900))
                    recs.append(('s%d_%d' % (i, j), rnd.choice(['', 'a comment', 'taxid 1234']), sq))
                files.append({'name': 'f%d.fa' % i, 'fasta': fasta(recs)})
            args = ['sketch', '-o', 'o'] + al
            if '-k' not in al and '-a' not in al: args += ['-k', str(rnd.choice([4, 8, 12, 16, 17, 21, 24, 32]))]
            args += ['-s', str(rnd.choice([1, 5, 20, 50, 200, 1000])), '-S', str(rnd.choice([42, 42, 1, 12345, 0]))]
            mode = rnd.random()
            if mode < 0.12: args.append('-i')
            elif mode < 0.4: args += ['-r'] + (['-m', str(rnd.choice([1, 2, 3]))] if rnd.random() < 0.5 else []) + (['-c', str(rnd.choice([1, 3, 8]))] if rnd.random() < 0.3 else []) + (['-g', str(rnd.choice([1000, 100000]))] if rnd.random() < 0.3 else [])
            elif mode < 0.5: args += ['-b', rnd.choice(['40', '1K', '100K', '2M'])]
            if rnd.random() < 0.2: args += ['-I', 'myid']
            names = [f['name'] for f in files]
            write_files(d, files)
            o = os.path.join(d, 'o.msh')
            if os.path.exists(o): os.remove(o)
            p = run(d, args + names)
            if p.returncode != 0 or not os.path.exists(o): continue
            info = run(d, ['info', 'o.msh']).stdout
            rows.append({'id': len(rows), 'kind': 'sketch', 'args': args + names, 'files': files, 'qmsh': None, 'out': info, 'msh': base64.b64encode(open(o, 'rb').read()).decode()})

        # ------------------------------------------------------------------ dist / triangle
        while sum(1 for r in rows if r['kind'] in ('dist', 'triangle')) < n_dist:
            al = options()
            gen = seq_for(al)
            fam = gen(rnd.randint(300, 2500), 0.0)
            def member(rate):
                if '-a' in al:
                    return mutate_p(fam, rate)
                return ''.join(rnd.choice(fam) if rnd.random() < rate else c for c in fam)
            nfiles = rnd.randint(2, 5)
            files = []
            for i in range(nfiles):
                recs = []
                for j in range(rnd.choice([1, 1, 2, 3])):
                    recs.append(('d%d_%d' % (i, j), rnd.choice(['', 'c%d' % i]), member(rnd.choice([0.01, 0.1, 0.3, 0.6])) if rnd.random() < 0.85 else gen(rnd.randint(30, 1500))))
                files.append({'name': 'g%d.fa' % i, 'fasta': fasta(recs)})
            write_files(d, files)
            kind = rnd.choice(['dist', 'triangle'])
            args = [kind] + al
            if '-k' not in al and '-a' not in al: args += ['-k', str(rnd.choice([4, 8, 12, 16, 17, 21, 24, 32]))]
            args += ['-s', str(rnd.choice([5, 50, 500, 1000, 2000]))]
            if rnd.random() < 0.15: args.append('-i')
            elif rnd.random() < 0.15: args += ['-r']
            if kind == 'dist':
                if rnd.random() < 0.3: args.append('-t')
                if rnd.random() < 0.2: args += ['-v', rnd.choice(['0.05', '1e-5'])]
                if rnd.random() < 0.2: args += ['-d', rnd.choice(['0.3', '0.6'])]
                if rnd.random() < 0.2: args.append('-C')
            else:
                if rnd.random() < 0.25: args.append('-E')
                if rnd.random() < 0.15: args += ['-d', rnd.choice(['0.3', '0.6'])]
                if rnd.random() < 0.15: args.append('-C')
            names = [f['name'] for f in files]
            p = run(d, args + names)
            if p.returncode != 0: continue
            rows.append({'id': len(rows), 'kind': kind, 'args': args + names, 'files': files, 'qmsh': None, 'out': p.stdout, 'msh': None})

        # ------------------------------------------------------------------ screen
        while sum(1 for r in rows if r['kind'] == 'screen') < n_screen:
            protein = rnd.random() < 0.65
            al = ['-a'] if protein else options()
            if protein and rnd.random() < 0.5: al += ['-k', str(rnd.choice([4, 5, 6, 7, 8, 9, 10, 12]))]
            if protein and rnd.random() < 0.3: al.append('-Z')
            if not protein and '-a' in al: al.remove('-a')
            gen = seq_for(al)
            nq = rnd.choice([1, 2, 3, 4])
            fam = gen(rnd.randint(200, 1500), 0.0)
            qfiles = []
            genomes = []
            for i in range(nq):
                g = (mutate_p(fam, rnd.choice([0.0, 0.05, 0.3])) if protein else ''.join(rnd.choice(fam) if rnd.random() < 0.05 else c for c in fam)) if rnd.random() < 0.8 else gen(rnd.randint(200, 1500))
                g = g.replace('*', 'L') if protein else g
                genomes.append(g)
                qfiles.append({'name': 'q%d.fa' % (i + 1), 'fasta': fasta([('q%d' % i, rnd.choice(['', 'desc']), g)])})
            mfiles = []
            big = sum(1 for r in rows if r['kind'] == 'screen' and r.get('big')) < 6 and rnd.random() < 0.04
            for m in range(rnd.choice([1, 1, 2, 3]) if not big else rnd.choice([1, 2])):
                recs = []
                for r in range(rnd.randint(1, 30) if not big else rnd.randint(60, 140)):
                    src = rnd.choice(genomes) if rnd.random() < 0.8 else gen(rnd.randint(100, 600))
                    if protein:
                        st = rnd.randint(0, max(0, len(src) - 60)); pr = src[st:st + rnd.choice([20, 40, 60, 120])]
                        dn = back(pr.upper() if rnd.random() < 0.9 else pr)
                        if rnd.random() < 0.3: dn = dn[rnd.randint(0, 2):]
                        if rnd.random() < 0.5: dn = revcomp(dn)
                        if rnd.random() < 0.15: dn = dn.lower()
                        if rnd.random() < 0.1: dn = dn[:rnd.randint(0, len(dn))] + 'N' + dn[rnd.randint(0, len(dn)):]
                    else:
                        st = rnd.randint(0, max(0, len(src) - 80)); dn = src[st:st + rnd.choice([30, 60, 100])]
                    recs.append(('r%d_%d' % (m, r), '', dn))
                if big:                      # > 2^20 characters: the program hashes the mixture in chunks (a repeated block keeps the file small)
                    total = sum(len(x[2]) + 1 for x in recs)
                    recs = recs * (1100000 // total + 1 + rnd.randint(0, 2)) if m == 0 else recs * rnd.randint(1, 5)
                mfiles.append({'name': 'm%d.fa' % (m + 1), 'fasta': fasta(recs)})
            if rnd.random() < 0.06:
                mfiles.append({'name': 'empty.fa', 'fasta': ''})
            write_files(d, qfiles + mfiles)
            sk = ['sketch', '-o', 'q'] + al
            if '-k' not in al and '-a' not in al: sk += ['-k', str(rnd.choice([4, 8, 12, 16, 21]))]
            sk += ['-s', str(rnd.choice([5, 50, 200, 1000]))] + (['-i'] if rnd.random() < 0.15 else [])
            qn = [f['name'] for f in qfiles]
            p = run(d, sk + qn)
            if p.returncode != 0 or not os.path.exists(os.path.join(d, 'q.msh')): continue
            sc = ['screen'] + (['-w'] if rnd.random() < 0.35 else []) + ['-i', str(rnd.choice([0, 0, -1, 0.3, 0.9])), '-v', str(rnd.choice([1.0, 1.0, 0.05, 1e-6])), 'q.msh']
            mn = [f['name'] for f in mfiles]
            p = run(d, sc + mn)
            if p.returncode != 0: continue
            rows.append({'id': len(rows), 'kind': 'screen', 'args': sc + mn, 'files': mfiles, 'qmsh': base64.b64encode(open(os.path.join(d, 'q.msh'), 'rb').read()).decode(), 'out': p.stdout, 'msh': None, 'sketchargs': sk + qn, 'qfiles': qfiles, 'big': big})
            os.remove(os.path.join(d, 'q.msh'))
    for name, kinds in (('sketch', ('sketch',)), ('dist', ('dist', 'triangle')), ('screen', ('screen',))):
        with gzip.open(out + '_' + name + '.jsonl.gz', 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
            for r in rows:
                if r['kind'] in kinds: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows', {k: sum(1 for r in rows if r['kind'] == k) for k in ('sketch', 'dist', 'triangle', 'screen')})


if __name__ == '__main__':
    a = sys.argv
    main(a[1], int(a[2]) if len(a) > 2 else 300, int(a[3]) if len(a) > 3 else 300, int(a[4]) if len(a) > 4 else 300, int(a[5]) if len(a) > 5 else 20261221)
