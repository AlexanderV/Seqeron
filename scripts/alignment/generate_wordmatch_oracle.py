#!/usr/bin/env python3
"""A2-28d (B09 F198) oracle: the files of EMBOSS 6.6.0 `wordmatch` (exact matches of a pattern set against large sequences, Rabin-Karp): the alignments (-outfile, default -aformat match), the statistics log
(-logfile), the GFF3 feature files -aoutfeat / -boutfeat; -wordsize 2 - 12, -nodumpalign, -nodumpfeat, pattern sets of 1 - 4 sequences (some with the same name, with N / X), several scanned sequences, nucleotide and
protein. The run date is faked (scripts/alignment/faketime.c).
JSONL.GZ rows: id, nuc, fasta1, fasta2, wordsize, dumpalign, dumpfeat, aformat, cmd, epoch, align, log, feat1, feat2, error, rc.
Usage: python3 scripts/alignment/generate_wordmatch_oracle.py OUT.jsonl.gz [n=260] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

NUC = 'ACGT'; AA = 'ACDEFGHIKLMNPQRSTVWY'


def fasta(name, s, desc=''):
    return '>%s%s\n%s\n' % (name, (' ' + desc) if desc else '', '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)))


def case(rnd):
    nuc = rnd.random() < 0.75
    alpha = NUC if nuc else AA
    genomes = [''.join(rnd.choice(alpha) for _ in range(rnd.randint(150, 700))) for _ in range(rnd.choice([1, 1, 2, 3]))]
    pats = []
    for i in range(rnd.choice([1, 2, 2, 3, 4])):
        g = rnd.choice(genomes); a = rnd.randrange(len(g) - 30); s = g[a:a + rnd.randint(25, 120)]
        if rnd.random() < 0.7:
            s = list(s)
            for _ in range(rnd.randint(0, 4)): s[rnd.randrange(len(s))] = rnd.choice(alpha)
            s = ''.join(s)
        if rnd.random() < 0.3: s = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 20))) + s
        if nuc and rnd.random() < 0.15:
            p = rnd.randrange(len(s)); s = s[:p] + 'N' * rnd.randint(1, 3) + s[p + 3:]
        if not nuc and rnd.random() < 0.1:
            p = rnd.randrange(len(s)); s = s[:p] + 'X' + s[p + 1:]
        if nuc and rnd.random() < 0.2: s = ''.join(c.lower() if rnd.random() < 0.4 else c for c in s)
        pats.append(s)
    if rnd.random() < 0.08: pats = [''.join(rnd.choice(alpha) for _ in range(rnd.randint(25, 80))) for _ in pats]      # unrelated: few or no matches
    names = ['pat%d' % i for i in range(len(pats))]
    if len(names) > 1 and rnd.random() < 0.08: names[1] = names[0]
    gnames = ['genome%d' % i for i in range(len(genomes))]
    return nuc, [(n, p) for n, p in zip(names, pats)], [(n, g) for n, g in zip(gnames, genomes)]


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            nuc, pats, genomes = case(rnd)
            f1 = ''.join(fasta(a, b, 'x' if rnd.random() < 0.2 else '') for a, b in pats); f2 = ''.join(fasta(a, b) for a, b in genomes)
            ws = rnd.choice([2, 3, 4, 4, 5, 6, 7, 8, 8, 10, 12]); da = rnd.random() > 0.12; dfeat = rnd.random() > 0.12
            aformat = rnd.choice(['match'] * 6 + ['pair', 'simple', 'fasta', 'srspair', 'markx3', 'score', 'srs'])
            args = ['a.fa', 'b.fa', '-wordsize', str(ws)] + ([] if da else ['-nodumpalign']) + ([] if dfeat else ['-nodumpfeat']) + (['-aformat', aformat] if aformat != 'match' else []) \
                + (['-outfile', 'o.match'] if da else []) + (['-aoutfeat', 'a.gff', '-boutfeat', 'b.gff'] if dfeat else []) + ['-logfile', 'log.txt', '-auto']
            for nme in ('o.match', 'log.txt', 'a.gff', 'b.gff'):
                if os.path.exists(os.path.join(d, nme)): os.remove(os.path.join(d, nme))
            open(os.path.join(d, 'a.fa'), 'w', newline='\n').write(f1); open(os.path.join(d, 'b.fa'), 'w', newline='\n').write(f2)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            try:
                p = subprocess.run(['wordmatch'] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=30)
            except subprocess.TimeoutExpired:
                continue
            def read(nme):
                q = os.path.join(d, nme)
                return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) else ''
            entries = ['[-asequence] a.fa', '[-bsequence] b.fa', '-wordsize %d' % ws] + ([] if da else ['-nodumpalign']) + ([] if dfeat else ['-nodumpfeat']) + (['-aformat %s' % aformat] if aformat != 'match' else []) \
                + (['-outfile o.match'] if da else []) + (['-aoutfeat a.gff', '-boutfeat b.gff'] if dfeat else []) + ['-logfile log.txt', '-auto']
            rows.append({'id': len(rows), 'nuc': nuc, 'fasta1': f1, 'fasta2': f2, 'wordsize': ws, 'dumpalign': da, 'dumpfeat': dfeat, 'aformat': aformat, 'cmd': entries, 'epoch': epoch,
                         'align': read('o.match'), 'log': read('log.txt'), 'feat1': read('a.gff'), 'feat2': read('b.gff'), 'error': p.stderr.decode('latin-1').strip() if isinstance(p.stderr, bytes) else (p.stderr.strip() or None), 'rc': p.returncode})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 260, int(sys.argv[3]) if len(sys.argv) > 3 else 20261201)
