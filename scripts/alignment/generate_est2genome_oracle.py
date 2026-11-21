#!/usr/bin/env python3
"""A2-28d oracle: EMBOSS 6.6.0 `est2genome` (spliced alignment of ESTs to a genomic sequence) vs EmbossEst2genomePrograms.

Cases: a genomic sequence of 2 - 6 exons separated by introns (GT..AG, CT..AC = reverse gene, or non-canonical), ESTs built from the exons (substitutions, small indels, a poly-A tail, partial, reverse-complemented,
N / IUPAC / lower case, unrelated ESTs), 1 - 3 ESTs per file, a description line, and the options: -match -mismatch -gappenalty -intronpenalty -splicepenalty -minscore -reverse -nousesplice -mode -nobest -space
(small values force the linear-space recursion) -shuffle -seed -align -width, plus sequence qualifiers (-sbegin1 / -send1 / -sreverse1 / -slower1 / -supper1 for the ESTs, -sbegin2 / -send2 / ... for the genome).
Output: JSON lines {"est", "genome", "args", "out", "rc", "err"} (out = the -outfile text, null when the program did not finish; rc as returned, -11 = SIGSEGV).
Usage: python3 scripts/alignment/generate_est2genome_oracle.py [n=300] [seed=20261123] | gzip > est2genome_oracle.jsonl.gz
"""
import json, os, random, subprocess, sys, tempfile

COMP = str.maketrans('ACGTacgtRYKMBDHVrykmbdhvNn', 'TGCAtgcaYRMKVHDByrmkvhdbNn')


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    seed = int(sys.argv[2]) if len(sys.argv) > 2 else 20261123
    rnd = random.Random(seed)
    work = tempfile.mkdtemp()

    def rseq(k, gc=0.5):
        return ''.join(rnd.choice('GC') if rnd.random() < gc else rnd.choice('AT') for _ in range(k))

    def rc(s):
        return s.translate(COMP)[::-1]

    def mutate(s, rate):
        out = []
        for c in s:
            r = rnd.random()
            if r < rate * 0.6: out.append(rnd.choice('ACGT'))
            elif r < rate * 0.8: pass
            elif r < rate: out.append(c); out.append(rnd.choice('ACGT') * rnd.choice([1, 1, 2, 3]))
            else: out.append(c)
        return ''.join(out)

    def fasta(name, desc, seq, width=60):
        return '>' + name + (' ' + desc if desc else '') + '\n' + ''.join(seq[i:i + width] + '\n' for i in range(0, len(seq), width))

    def make_case():
        nex = rnd.randint(2, 6)
        gc = rnd.choice([0.4, 0.5, 0.6])
        exons = [rseq(rnd.randint(25, 300), gc) for _ in range(nex)]
        kind = rnd.choice(['fwd', 'fwd', 'rev', 'plain'])
        introns = []
        for _ in range(nex - 1):
            body = rseq(rnd.randint(20, 400), 0.4)
            if kind == 'fwd': introns.append('GT' + body + 'AG')
            elif kind == 'rev': introns.append('CT' + body + 'AC')
            else: introns.append(body)
        genome = rseq(rnd.randint(0, 200), 0.45)
        for i, e in enumerate(exons):
            genome += e
            if i < nex - 1: genome += introns[i]
        genome += rseq(rnd.randint(0, 200), 0.45)
        ests = []
        for _ in range(rnd.choice([1, 1, 1, 2, 3])):
            sel = list(range(nex))
            if rnd.random() < 0.25 and nex > 2: sel.pop(rnd.randrange(1, nex - 1))                  # skipped exon
            est = ''.join(exons[i] for i in sel)
            est = mutate(est, rnd.choice([0, 0, 0.01, 0.03, 0.08]))
            if rnd.random() < 0.3: est += 'A' * rnd.randint(5, 30)
            if rnd.random() < 0.3: est = est[rnd.randint(0, len(est) // 3):]
            if rnd.random() < 0.3: est = est[:len(est) - rnd.randint(0, len(est) // 3)]
            if rnd.random() < 0.35: est = rc(est)
            if rnd.random() < 0.08: est = rseq(rnd.randint(40, 200))                                  # unrelated
            if rnd.random() < 0.12:
                i = rnd.randrange(len(est) + 1); est = est[:i] + rnd.choice(['N', 'NN', 'R', 'Y', 'NNNNN']) + est[i:]
            if rnd.random() < 0.1: est = est.lower()
            ests.append(est or 'ACGT')
        if rnd.random() < 0.08: genome = genome.lower()
        if rnd.random() < 0.08:
            i = rnd.randrange(len(genome)); genome = genome[:i] + rnd.choice(['N', 'NNN', 'R']) + genome[i:]
        return genome, ests

    rows = 0; errors = 0
    for case in range(n):
        genome, ests = make_case()
        g = fasta(rnd.choice(['genomic', 'chr1', 'gene_with_a_long_name', 'g']), rnd.choice(['', 'a genomic region', 'locus 1 length']), genome, rnd.choice([60, 70, 1000000]))
        e = ''.join(fasta(rnd.choice(['est%d' % (i + 1), 'EST_%d_longer' % (i + 1), 'e%d' % (i + 1)]), rnd.choice(['', 'an est', 'cDNA clone ' + str(i)]), s, rnd.choice([50, 60, 1000000])) for i, s in enumerate(ests))
        args = []
        if rnd.random() < 0.25: args += ['-match', str(rnd.choice([1, 2, 3, 5]))]
        if rnd.random() < 0.25: args += ['-mismatch', str(rnd.choice([1, 2, 4, -3, 0]))]
        if rnd.random() < 0.25: args += ['-gappenalty', str(rnd.choice([1, 2, 3, 6, -2]))]
        if rnd.random() < 0.25: args += ['-intronpenalty', str(rnd.choice([10, 20, 40, 60, 100]))]
        if rnd.random() < 0.25: args += ['-splicepenalty', str(rnd.choice([5, 10, 20, 30, 50]))]
        if rnd.random() < 0.3: args += ['-minscore', str(rnd.choice([0, 10, 30, 60, 150]))]
        if rnd.random() < 0.2: args += ['-reverse']
        if rnd.random() < 0.15: args += ['-nousesplice']
        r = rnd.random()
        if r < 0.15: args += ['-mode', 'forward']
        elif r < 0.3: args += ['-mode', 'reverse']
        elif r < 0.35: args += ['-mode', rnd.choice(['f', 'REV', 'b', 'Both'])]
        if rnd.random() < 0.2: args += ['-nobest']
        if rnd.random() < 0.35: args += ['-space', rnd.choice(['0.0005', '0.001', '0.002', '0.01', '0.05', '1', '10'])]
        if rnd.random() < 0.08: args += ['-shuffle', str(rnd.choice([1, 2, 5]))]
        if rnd.random() < 0.05: args += ['-seed', str(rnd.choice([0, 1, 777, -5]))]
        if rnd.random() < 0.4: args += ['-align']
        if rnd.random() < 0.3: args += ['-width', str(rnd.choice([10, 30, 50, 60, 100, 5000]))]
        if rnd.random() < 0.12: args += ['-sbegin1', str(rnd.choice([1, 5, 20])), '-send1', str(rnd.choice([-1, -10, 100, 500]))]
        if rnd.random() < 0.1: args += ['-sbegin2', str(rnd.choice([1, 50, 100])), '-send2', str(rnd.choice([-1, -50, 600, 2000]))]
        if rnd.random() < 0.06: args += [rnd.choice(['-sreverse1', '-slower1', '-supper1', '-sreverse2', '-slower2', '-supper2'])]
        open(work + '/e.fa', 'w').write(e); open(work + '/g.fa', 'w').write(g)
        if os.path.exists(work + '/o.txt'): os.remove(work + '/o.txt')
        try:
            p = subprocess.run(['est2genome', '-estsequence', 'e.fa', '-genomesequence', 'g.fa', '-outfile', 'o.txt', '-auto'] + args, cwd=work, capture_output=True, text=True, timeout=60)
        except subprocess.TimeoutExpired:
            continue
        out = open(work + '/o.txt').read() if p.returncode == 0 and os.path.exists(work + '/o.txt') else None
        errors += out is None
        print(json.dumps({'est': e, 'genome': g, 'args': args, 'out': out, 'rc': p.returncode, 'err': p.stderr}))
        rows += 1
    sys.stderr.write('rows=%d errors=%d\n' % (rows, errors))


main()
