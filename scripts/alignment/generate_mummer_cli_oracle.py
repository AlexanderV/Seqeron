#!/usr/bin/env python3
"""A3-9e oracle: the `mummer` command (MUMmer 3.23, Kurtz maxmat3) on random reference / query files: every matching mode (-mum -mumcand -mumreference -maxmatch), -n, -l, -b, -r, -s, -c, -F, -L, several query
files, multi-record references with white space in the headers, empty records, text before the first header, upper / lower case, CR LF line ends, N / IUPAC / odd characters, and the error runs.
Output: JSON lines {"ref", "qrys", "args", "out", "rc", "err"} (stdout and stderr decoded as latin-1; stderr without the three timing lines).
Usage: python3 scripts/alignment/generate_mummer_cli_oracle.py [n=700] [seed=20261130] | gzip > mummer_cli_oracle.jsonl.gz"""
import json, os, random, subprocess, sys, tempfile

def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 700
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261130
    rnd = random.Random(seed)
    work = tempfile.mkdtemp(prefix='b09mm')

    def rs(k): return ''.join(rnd.choice('acgt') for _ in range(k))

    def mutate(s, rate):
        o = []
        for c in s:
            r = rnd.random()
            if r < rate * .6: o.append(rnd.choice('acgt'))
            elif r < rate * .8: pass
            elif r < rate: o.append(c + rnd.choice('acgt'))
            else: o.append(c)
        return ''.join(o)

    def rcomp(s): return s.translate(str.maketrans('acgtACGT', 'tgcaTGCA'))[::-1]

    def noise(t):
        o = list(t)
        for _ in range(rnd.choice([0, 0, 0, 2, 6])):
            if o: o[rnd.randrange(len(o))] = rnd.choice('nNRYKMSWBDHVx-*XJ0.')
        return ''.join(o)

    def fasta(recs, eol, width, junk_first=False):
        out = []
        if junk_first: out.append('acgt' + eol)
        for name, seq in recs:
            out.append('>' + name + eol)
            if width <= 0: out.append(seq + eol)
            else:
                for i in range(0, len(seq), width): out.append(seq[i:i + width] + eol)
        return ''.join(out)

    rows = 0
    stats = {}
    for case in range(n):
        base = rs(rnd.randint(40, 400))
        nref = rnd.choice([1, 1, 1, 2, 3, 5])
        refrecs = []
        for i in range(nref):
            seq = mutate(base[rnd.randint(0, 20):rnd.randint(len(base) // 2, len(base))], 0.02) + rs(rnd.randint(0, 40))
            if rnd.random() < 0.3: seq = mutate(base, 0.01)
            if rnd.random() < 0.15: seq = seq + seq[:rnd.randint(5, 30)]
            if rnd.random() < 0.3: seq = noise(seq)
            if rnd.random() < 0.2: seq = seq.upper()
            if rnd.random() < 0.05: seq = ''
            nm = rnd.choice(['r%d' % i, 'ref_%d description text' % i, 'longer_reference_name_%d x' % (i * 7), 'a', 'chr%d\tTAB' % i])
            refrecs.append((nm, seq))
        nq = rnd.choice([1, 1, 1, 2, 3])
        qfiles = []
        for f in range(rnd.choice([1, 1, 1, 2])):
            qrecs = []
            for i in range(nq):
                parts = [rs(rnd.randint(0, 20))]
                for _ in range(rnd.randint(1, 3)):
                    piece = mutate(base[rnd.randint(0, len(base) // 2):rnd.randint(len(base) // 2, len(base))], rnd.choice([0, .01, .05, .1]))
                    if rnd.random() < .4: piece = rcomp(piece)
                    parts.append(piece); parts.append(rs(rnd.randint(0, 20)))
                seq = ''.join(parts)
                if rnd.random() < 0.3: seq = noise(seq)
                if rnd.random() < 0.15: seq = seq.upper()
                if rnd.random() < 0.05: seq = ''
                if rnd.random() < 0.05: seq = seq[:rnd.randint(0, 12)]
                qrecs.append(('q%d_%d some words' % (f, i), seq))
            qfiles.append(qrecs)
        eol = rnd.choice(['\n', '\n', '\n', '\r\n'])
        ref = fasta(refrecs, eol, rnd.choice([0, 60, 70, 13]), rnd.random() < 0.04)
        qrys = [fasta(q, eol, rnd.choice([0, 60, 70, 17])) for q in qfiles]
        args = []
        r = rnd.random()
        if r < .45: args.append('-maxmatch')
        elif r < .65: args.append('-mum')
        elif r < .8: args.append('-mumreference')
        elif r < .9: args.append('-mumcand')
        if rnd.random() < .6: args.append('-n')
        args += ['-l', str(rnd.choice([6, 8, 10, 12, 15, 20, 30]) if '-maxmatch' in args else rnd.choice([1, 3, 5, 8, 12, 15, 20, 30]))] if rnd.random() < .9 else []
        r = rnd.random()
        if r < .35: args.append('-b')
        elif r < .5: args.append('-r')
        if rnd.random() < .2: args.append('-s')
        if rnd.random() < .25 and ('-b' in args or '-r' in args): args.append('-c')
        if rnd.random() < .2: args.append('-F')
        if rnd.random() < .2: args.append('-L')
        rnd.shuffle(args) if False else None
        rf = os.path.join(work, 'r.fa'); open(rf, 'w', newline='', encoding='latin-1').write(ref)
        qf = []
        for k, q in enumerate(qrys):
            p = os.path.join(work, 'q%d.fa' % k); open(p, 'w', newline='', encoding='latin-1').write(q); qf.append(p)
        # args are given as flags then files; keep -l value together
        cmd = ['mummer'] + args + [rf] + qf
        p = subprocess.run(cmd, capture_output=True, timeout=120)
        err = '\n'.join(l for l in p.stderr.decode('latin-1').split('\n') if not l.startswith(('# CONSTRUCTIONTIME', '# COMPLETETIME', '# SPACE')))
        # file names in stderr are the temporary paths; map them to the fixed names of the port
        err = err.replace(work + '/r.fa', 'reference.fa')
        for k in range(len(qf)): err = err.replace(work + '/q%d.fa' % k, 'query%d.fa' % k)
        err = err.replace(work, '')
        row = {'ref': ref, 'qrys': qrys, 'args': args, 'out': p.stdout.decode('latin-1'), 'rc': p.returncode, 'err': err}
        stats[p.returncode] = stats.get(p.returncode, 0) + 1
        print(json.dumps(row))
        rows += 1
    sys.stderr.write('rows=%d rc=%s\n' % (rows, stats))

main()
