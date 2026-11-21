#!/usr/bin/env python3
"""A5-8c3 (B09 F196) oracle: the files of EMBOSS 6.6.0 `diffseq` (differences between two nearly identical sequences): the report (7 report tags) in EVERY -rformat, the GFF3 feature files -aoutfeat / -boutfeat,
for nucleotide and protein pairs (SNPs, runs of mismatches, insertions and deletions, flanks that make the ends differ), -wordsize 2 - 14, -globaldifferences, the report options -raccshow -rdesshow -rscoreshow
-rstrandshow -rusashow -rmaxall -rmaxseq. The run date is faked (scripts/alignment/faketime.c).
JSONL.GZ rows: id, nuc, fasta1, fasta2, wordsize, glob, rq (the report qualifiers text incl. -rformat), cmd (header entries), epoch, report, feat1, feat2 (file texts; '' when the program wrote none), error, rc.
Usage: python3 scripts/alignment/generate_diffseq_oracle.py OUT.jsonl.gz [per_format=14] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_emboss_report_formats_oracle as G  # noqa: E402

NUC = 'ACGT'; AA = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(rnd, s, alpha, n):
    s = list(s)
    for _ in range(n):
        k = rnd.random(); p = rnd.randrange(len(s))
        if k < 0.45: s[p] = rnd.choice(alpha)
        elif k < 0.7: del s[p:p + rnd.randint(1, 4)]
        elif k < 0.9: s[p:p] = [rnd.choice(alpha) for _ in range(rnd.randint(1, 5))]
        else:
            for q in range(p, min(len(s), p + rnd.randint(2, 6))): s[q] = rnd.choice(alpha)
    return ''.join(s)


def case(rnd):
    nuc = rnd.random() < 0.7
    alpha = NUC if nuc else AA
    base = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(40, 400)))
    other = mutate(rnd, base, alpha, rnd.randint(0, 8))
    if rnd.random() < 0.4: other = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 25))) + other
    if rnd.random() < 0.4: other = other + ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 25)))
    if rnd.random() < 0.3: base = base[rnd.randint(0, 20):len(base) - rnd.randint(0, 20)]
    if nuc and rnd.random() < 0.15: other = ''.join(c.lower() if rnd.random() < 0.3 else c for c in other)
    if nuc and rnd.random() < 0.1:
        p = rnd.randrange(len(other)); other = other[:p] + 'N' * rnd.randint(1, 3) + other[p + 3:]
    if rnd.random() < 0.05: other = ''.join(rnd.choice(alpha) for _ in range(len(other)))        # unrelated: no match
    return nuc, base, other


def fasta(name, desc, s):
    return '>%s%s\n%s\n' % (name, (' ' + desc) if desc else '', '\n'.join(s[i:i + 60] for i in range(0, len(s), 60)))


def main(out, per, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        for fmt in G.FORMATS:
            made = 0; tries = 0
            while made < per and tries < per * 5:
                tries += 1
                nuc, a, b = case(rnd)
                if not a or not b: continue
                n1 = rnd.choice(['seqa', 'X56789.2', 'gb_a', 'Q12345']); n2 = rnd.choice(['seqb', 'sample2', 'chr7', 'P99999'])
                f1 = fasta(n1, rnd.choice(['', 'first sequence', 'OS=Homo sapiens']), a); f2 = fasta(n2, rnd.choice(['', 'second one']), b)
                ws = rnd.choice([2, 3, 4, 5, 6, 8, 10, 10, 12, 14]); glob = rnd.random() < 0.4
                rq = ['-rformat', fmt]
                if rnd.random() < 0.25: rq.append('-raccshow')
                if rnd.random() < 0.25: rq.append('-rdesshow')
                if rnd.random() < 0.15: rq.append(rnd.choice(['-norscoreshow', '-rscoreshow']))
                if rnd.random() < 0.15: rq.append(rnd.choice(['-norstrandshow', '-rstrandshow']))
                if rnd.random() < 0.15: rq.append('-rusashow')
                if rnd.random() < 0.15: rq += ['-rmaxall', str(rnd.randint(1, 4))]
                if rnd.random() < 0.15: rq += ['-rmaxseq', str(rnd.randint(1, 3))]
                args = ['a.fa', 'b.fa', '-wordsize', str(ws)] + (['-globaldifferences'] if glob else []) + rq + ['-outfile', 'o.rep', '-aoutfeat', 'a.gff', '-boutfeat', 'b.gff', '-auto']
                for n in ('o.rep', 'a.gff', 'b.gff'):
                    if os.path.exists(os.path.join(d, n)): os.remove(os.path.join(d, n))
                open(os.path.join(d, 'a.fa'), 'w', newline='\n').write(f1); open(os.path.join(d, 'b.fa'), 'w', newline='\n').write(f2)
                epoch = rnd.randint(978307200, 2051222400)
                env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
                try:
                    p = subprocess.run(['diffseq'] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=20)
                except subprocess.TimeoutExpired:
                    continue
                def read(n):
                    q = os.path.join(d, n)
                    return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) else ''
                entries = ['[-asequence] a.fa', '[-bsequence] b.fa', '-wordsize %d' % ws] + (['-globaldifferences'] if glob else [])
                i = 0
                while i < len(rq):
                    if i + 1 < len(rq) and not rq[i + 1].startswith('-'): entries.append(rq[i] + ' ' + rq[i + 1]); i += 2
                    else: entries.append(rq[i]); i += 1
                entries += ['-outfile o.rep', '-aoutfeat a.gff', '-boutfeat b.gff', '-auto']
                rows.append({'id': len(rows), 'nuc': nuc, 'fasta1': f1, 'fasta2': f2, 'wordsize': ws, 'glob': glob, 'rq': ' '.join(rq), 'cmd': entries, 'epoch': epoch,
                             'report': read('o.rep'), 'feat1': read('a.gff'), 'feat2': read('b.gff'), 'error': p.stderr.strip() or None, 'rc': p.returncode})
                made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 14, int(sys.argv[3]) if len(sys.argv) > 3 else 20261129)
