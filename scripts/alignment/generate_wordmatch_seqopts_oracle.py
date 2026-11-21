#!/usr/bin/env python3
"""A5-8c3e (B09 F239) oracle: EMBOSS 6.6.0 `wordmatch` with the SEQUENCE qualifiers of its two inputs (-sbegin1/2 -send1/2 -sreverse1/2 -slower1/2 -supper1/2 -snucleotide1/2 -sprotein1/2 -scircular1/2; the first input is a
seqset, the second a seqall): the alignments (-outfile, several -aformat), the statistics (-logfile), the GFF3 feature files -aoutfeat / -boutfeat the program wrote, byte for byte, or its message when it refuses.
Same row schema as generate_wordmatch_featout_oracle.py (`qual` = the qualifiers typed, `aformat`).
Usage: python3 scripts/alignment/generate_wordmatch_seqopts_oracle.py OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_wordmatch_oracle as W  # noqa: E402


def quals(rnd, k, lens, nuc):
    q = []
    r = rnd.random()
    if r < 0.5:
        q += ['-sbegin%d' % k, str(rnd.choice([1, rnd.randint(1, max(1, min(lens))), rnd.randint(1, max(lens)), -rnd.randint(1, max(1, min(lens))), 0]))]
    r = rnd.random()
    if r < 0.5:
        q += ['-send%d' % k, str(rnd.choice([rnd.randint(1, max(lens)), max(lens) + 3, rnd.randint(max(1, min(lens) // 2), max(1, min(lens))), -rnd.randint(1, max(1, min(lens))), 0]))]
    if nuc and rnd.random() < 0.3:
        q += ['-sreverse%d' % k]
    r = rnd.random()
    if r < 0.1:
        q += ['-slower%d' % k]
    elif r < 0.2:
        q += ['-supper%d' % k]
    r = rnd.random()
    if r < 0.04:
        q += ['-snucleotide%d' % k]
    elif r < 0.08:
        q += ['-sprotein%d' % k]
    if rnd.random() < 0.03:
        q += ['-scircular%d' % k]
    return q


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            nuc, pats, genomes = W.case(rnd)
            f1 = ''.join(W.fasta(a, b, 'x' if rnd.random() < 0.2 else '') for a, b in pats); f2 = ''.join(W.fasta(a, b) for a, b in genomes)
            ws = rnd.choice([2, 3, 4, 4, 5, 6, 7, 8, 8, 10, 12]); da = rnd.random() > 0.1; dfeat = rnd.random() > 0.1
            aformat = rnd.choice(['match'] * 6 + ['pair', 'simple', 'fasta', 'srspair', 'markx3', 'score', 'srs'])
            rq = quals(rnd, 1, [len(p) for _, p in pats], nuc) + quals(rnd, 2, [len(g) for _, g in genomes], nuc)
            args = ['a.fa', 'b.fa', '-wordsize', str(ws)] + ([] if da else ['-nodumpalign']) + ([] if dfeat else ['-nodumpfeat']) + rq + (['-aformat', aformat] if aformat != 'match' else []) \
                + (['-outfile', 'o.match'] if da else []) + ['-aoutfeat', 'a.gff', '-boutfeat', 'b.gff'] + ['-logfile', 'log.txt', '-auto']
            for nme in os.listdir(d):
                if nme.endswith(('.out', '.match', '.fa', '.txt', '.gff')): os.remove(os.path.join(d, nme))
            open(os.path.join(d, 'a.fa'), 'w', newline='\n').write(f1); open(os.path.join(d, 'b.fa'), 'w', newline='\n').write(f2)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            try:
                p = subprocess.run(['wordmatch'] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=30)
            except subprocess.TimeoutExpired:
                continue

            def read(nme):
                q = os.path.join(d, nme)
                return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) else None
            written = sorted(x for x in os.listdir(d) if x not in ('a.fa', 'b.fa', 'o.match', 'log.txt', 'faketime.so', 'a.wordmatch'))
            rows.append({'id': len(rows), 'nuc': nuc, 'fasta1': f1, 'fasta2': f2, 'wordsize': ws, 'dumpalign': da, 'dumpfeat': dfeat, 'aformat': aformat, 'qual': ' '.join(rq), 'epoch': epoch,
                         'align': read('o.match') or '', 'log': read('log.txt') or '', 'files': {w: read(w) for w in written}, 'error': p.stderr.strip() or None, 'rc': p.returncode})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if r['error']), 'with a message')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261328)
