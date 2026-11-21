#!/usr/bin/env python3
"""A3-9b oracle: the whole `nucmer` 3.23 pipeline (prenuc | mummer | mgaps | postnuc) - the text of out.delta (or out.cluster with --nodelta) - vs GenomicAnalyzer.Nucmer.

Each case: a multi-record reference and query FASTA (related segments with substitutions / indels, inverted segments, repeats, N runs, lower case, IUPAC and non-letter characters,
records that share a segment, self comparison) and the nucmer options of every family (--mum / --mumreference / --maxmatch, -b, -c, -g, -l, -D, -d, --banded, --noextend, --nooptimize,
--nosimplify, -f, -r, --nodelta). The program is run in the fixed directory /tmp/b09nuc, so the header line of the output (the absolute paths) is reproducible.
Output: JSON lines  {"ref": fasta, "qry": fasta, "args": [..], "out": text | null, "error": bool}.
Usage: python3 scripts/alignment/generate_nucmer_oracle.py [n=360] [seed=20261121] | gzip > nucmer_oracle.jsonl.gz
"""
import json, os, random, subprocess, sys, shutil

WORK = '/tmp/b09nuc'
COMP = str.maketrans('ACGTacgtRYKMBDHVrykmbdhv', 'TGCAtgcaYRMKVHDByrmkvhdb')


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 360
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261121
    rnd = random.Random(seed)
    os.makedirs(WORK, exist_ok=True)

    def rseq(k, gc=0.5):
        return ''.join(rnd.choice('GC') if rnd.random() < gc else rnd.choice('AT') for _ in range(k))

    def mutate(s, rate):
        out = []
        for c in s:
            r = rnd.random()
            if r < rate * 0.5: out.append(rnd.choice('ACGT'))
            elif r < rate * 0.7: pass
            elif r < rate * 0.9: out.append(c); out.append(rnd.choice('ACGT') * rnd.choice([1, 1, 2, 5, 12]))
            else: out.append(c)
        return ''.join(out)

    def rc(s):
        return s.translate(COMP)[::-1]

    def fasta(records, width):
        out = []
        for name, seq in records:
            out.append('>' + name + '\n')
            for i in range(0, len(seq), width): out.append(seq[i:i + width] + '\n')
            if not seq: out.append('\n')
        return ''.join(out)

    def decorate(seq):
        r = rnd.random()
        if r < 0.15:
            i = rnd.randrange(len(seq) + 1); seq = seq[:i] + 'N' * rnd.choice([1, 3, 10, 40]) + seq[i:]
        if rnd.random() < 0.08:
            i = rnd.randrange(len(seq) + 1); seq = seq[:i] + rnd.choice('RYKMSWBDHV') + seq[i:]
        if rnd.random() < 0.05:
            i = rnd.randrange(len(seq) + 1); seq = seq[:i] + rnd.choice(['-', '*', '.']) + seq[i:]
        if rnd.random() < 0.12:
            seq = seq.lower()
        elif rnd.random() < 0.1:
            i = rnd.randrange(len(seq) + 1); j = min(len(seq), i + rnd.randint(5, 80)); seq = seq[:i] + seq[i:j].lower() + seq[j:]
        return seq

    def make_case():
        nseg = rnd.randint(2, 6)
        segs = [rseq(rnd.randint(60, 900), rnd.choice([0.3, 0.5, 0.5, 0.65])) for _ in range(nseg)]
        unit = rseq(rnd.randint(5, 40))
        nref = rnd.choice([1, 1, 1, 2, 3, 4]); nq = rnd.choice([1, 1, 2, 2, 3])
        refs = [[] for _ in range(nref)]; qs = [[] for _ in range(nq)]
        for s in segs:
            ri = rnd.randrange(nref); qi = rnd.randrange(nq)
            refs[ri].append(mutate(s, rnd.choice([0.0, 0.01, 0.02])))
            refs[ri].append(rseq(rnd.randint(0, 150)))
            if rnd.random() < 0.25: refs[rnd.randrange(nref)].append(mutate(s, 0.03))        # second copy
            if rnd.random() < 0.3: refs[ri].append(unit * rnd.randint(2, 6))
            piece = mutate(s, rnd.choice([0.0, 0.01, 0.03, 0.06, 0.1, 0.18]))
            if rnd.random() < 0.3: piece = rc(piece)
            if rnd.random() < 0.15 and len(piece) > 100:
                a = rnd.randint(20, len(piece) - 20); piece = piece[:a] + rseq(rnd.randint(1, 400)) + piece[a:]   # insertion between halves
            if rnd.random() < 0.1 and len(piece) > 100:
                a = rnd.randint(20, len(piece) - 60); piece = piece[:a] + piece[a + rnd.randint(5, 150):]      # deletion
            qs[qi].append(piece)
            qs[qi].append(rseq(rnd.randint(0, 150)))
            if rnd.random() < 0.25: qs[rnd.randrange(nq)].append(unit * rnd.randint(1, 4))
        rrec = []; qrec = []
        for i, parts in enumerate(refs):
            rnd.shuffle(parts) if rnd.random() < 0.2 else None
            rrec.append(('ref%d%s' % (i + 1, rnd.choice(['', ' desc text', '\tx y'])), decorate(''.join(parts))))
        for i, parts in enumerate(qs):
            qrec.append(('qry%d%s' % (i + 1, rnd.choice(['', ' sample', ' length=5'])), decorate(''.join(parts))))
        if rnd.random() < 0.06:
            qrec.append(('qempty', ''))
        if rnd.random() < 0.04:
            rrec.append(('refshort', rseq(rnd.randint(1, 12))))
        return rrec, qrec

    rows = 0; errors = 0
    edge = []
    s1 = rseq(300); s2 = rseq(300)
    edge.append(('', fasta([('q', s1)], 60), []))                                  # empty reference file
    edge.append((fasta([('r', '')], 60), fasta([('q', s1)], 60), []))              # reference record without sequence
    edge.append((fasta([('r', s1)], 60), '', []))                                  # empty query file
    edge.append((fasta([('r', s1)], 60), fasta([('q', '')], 60), []))              # query record without sequence
    edge.append((fasta([('r', s1)], 60), fasta([('q', s1)], 60), ['--mum', '--maxmatch']))   # two matching algorithms
    edge.append((fasta([('r', s1)], 60), fasta([('q', s1)], 60), ['--mumreference', '--mumcand']))
    edge.append((fasta([('r', 'N' * 200)], 60), fasta([('q', 'N' * 200)], 60), []))
    edge.append((fasta([('r', s1), ('r', s2)], 60), fasta([('q', s2)], 60), []))   # duplicate reference ids
    edge.append((fasta([('r', s1)], 60), fasta([('q', s1), ('q', s2)], 60), []))   # duplicate query ids
    edge.append((fasta([('r1', s1), ('r2', s1)], 60), fasta([('q', s1)], 60), ['--mum']))     # shared segment: no MUM
    edge.append((fasta([('r1', s1), ('r2', s1)], 60), fasta([('q', s1)], 60), ['--maxmatch']))
    edge.append((s1 + '\n', fasta([('q', s1)], 60), []))                           # reference without a header
    edge.append((fasta([('r', s1)], 60), fasta([('q', s1[:100] + s2[:100] + s1[100:])], 60), ['-b', '30']))
    edge.append((fasta([('r', s1)], 60), fasta([('q', rc(s1))], 60), ['-r', '--nooptimize']))
    for ref, qry, args in edge:
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
        open(WORK + '/ref.fa', 'w', newline='').write(ref); open(WORK + '/qry.fa', 'w', newline='').write(qry)
        p = subprocess.run(['nucmer'] + args + ['-p', 'out', 'ref.fa', 'qry.fa'], cwd=WORK, capture_output=True, text=True, timeout=120)
        out = open(WORK + '/out.delta').read() if p.returncode == 0 and os.path.exists(WORK + '/out.delta') else None
        errors += out is None
        print(json.dumps({'ref': ref, 'qry': qry, 'args': args, 'out': out, 'error': out is None}))
        rows += 1
    for case in range(n):
        rrec, qrec = make_case()
        w = rnd.choice([60, 70, 80, 1000000])
        ref = fasta(rrec, w); qry = fasta(qrec, rnd.choice([60, 70, 80]))
        if rnd.random() < 0.07:  # self comparison
            qry = ref
        args = []
        fam = rnd.random()
        args += [rnd.choice(['--mum', '--mumreference', '--maxmatch', '--mumcand', '--mumreference', '--maxmatch', '--mum'])] if rnd.random() < 0.8 else []
        if rnd.random() < 0.55: args += ['-b', str(rnd.choice([1, 20, 50, 100, 200, 500, 1000, 5000, 10000, 10001, 0, -5]))]
        if rnd.random() < 0.4: args += ['-c', str(rnd.choice([1, 20, 30, 65, 100, 200]))]
        if rnd.random() < 0.35: args += ['-g', str(rnd.choice([1, 20, 90, 200, 1000]))]
        if rnd.random() < 0.45: args += ['-l', str(rnd.choice([8, 10, 12, 15, 20, 30, 50]))]
        if rnd.random() < 0.3: args += ['-D', str(rnd.choice([0, 2, 5, 10, 30]))]
        if rnd.random() < 0.3: args += ['-d', str(rnd.choice([0.0, 0.05, 0.12, 0.3, 1.0]))]
        if rnd.random() < 0.2: args += ['--banded']
        if rnd.random() < 0.2: args += ['--noextend']
        if rnd.random() < 0.2: args += ['--nooptimize']
        if rnd.random() < 0.2: args += ['--nosimplify']
        r = rnd.random()
        if r < 0.1: args += ['-f']
        elif r < 0.2: args += ['-r']
        elif r < 0.25: args += ['-f', '-r']
        if rnd.random() < 0.05: args += ['--nodelta']
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
        open(WORK + '/ref.fa', 'w', newline='').write(ref); open(WORK + '/qry.fa', 'w', newline='').write(qry)
        p = subprocess.run(['nucmer'] + args + ['-p', 'out', 'ref.fa', 'qry.fa'], cwd=WORK, capture_output=True, text=True, timeout=120)
        out = None
        fn = 'out.cluster' if '--nodelta' in args else 'out.delta'
        if p.returncode == 0 and os.path.exists(WORK + '/' + fn):
            out = open(WORK + '/' + fn).read()
        err = out is None
        errors += err
        print(json.dumps({'ref': ref, 'qry': qry, 'args': args, 'out': out, 'error': err}))
        rows += 1
    sys.stderr.write('rows=%d errors=%d\n' % (rows, errors))


main()
