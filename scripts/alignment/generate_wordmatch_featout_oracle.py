#!/usr/bin/env python3
"""A5-8c3d (B09 F225) oracle: EMBOSS 6.6.0 `wordmatch` with the FEATURE OUTPUT qualifiers of both feature files (-offormat4/5 of every output feature format and its aliases, other case, an unknown name, the `format::file` form of
-aoutfeat / -boutfeat, -ofopenfile, -ofname / -ofextension / -ofsingle) on nucleotide and protein pattern sets (several tables per file, one per pattern / scanned sequence): the alignments (-outfile), the log and every file the
program wrote byte for byte. Same row schema as generate_wordmatch_oracle.py plus `names` (the -aoutfeat / -boutfeat values), `qual` (the feature qualifiers typed) and `files` (every file written, name -> text).
Usage: python3 scripts/alignment/generate_wordmatch_featout_oracle.py OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_wordmatch_oracle as W  # noqa: E402

FORMATS = ['gff', 'gff2', 'gff3', 'bed', 'embl', 'em', 'genbank', 'gb', 'ddbj', 'refseq', 'refseqp', 'pir', 'nbrf', 'swiss', 'sw', 'swissprot', 'uniprot', 'dasgff', 'draw', 'debug', 'unknown']
ODD = ['GFF3', 'Embl', 'GENBANK', 'bogus', 'fasta', 'gff4']


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            nuc, pats, genomes = W.case(rnd)
            f1 = ''.join(W.fasta(a, b, 'x' if rnd.random() < 0.2 else '') for a, b in pats); f2 = ''.join(W.fasta(a, b) for a, b in genomes)
            ws = rnd.choice([2, 3, 4, 4, 5, 6, 7, 8, 8, 10, 12])
            names = ['a.out', 'b.out']
            rq = []
            for k in (4, 5):
                f = rnd.choice(FORMATS) if rnd.random() < 0.8 else rnd.choice(ODD)
                if nuc and f in ('pir', 'nbrf', 'swiss', 'sw', 'swissprot', 'uniprot', 'refseqp') and rnd.random() < 0.6: f = rnd.choice(['embl', 'genbank', 'gff2', 'bed', 'dasgff'])
                if not nuc and f in ('embl', 'em', 'genbank', 'gb', 'ddbj', 'refseq') and rnd.random() < 0.6: f = rnd.choice(['pir', 'swiss', 'refseqp', 'gff2', 'bed', 'dasgff', 'draw'])
                r = rnd.random()
                if r < 0.07: names[k - 4] = '%s::%s' % (rnd.choice(FORMATS), names[k - 4])
                elif r < 0.1: rq += ['-ofopenfile%d' % k, 'f%d.out' % k]
                elif r < 0.13: rq += ['-ofname%d' % k, 'nm']
                elif r < 0.16: rq += ['-ofextension%d' % k, 'ext']
                elif r < 0.19: rq += ['-ofsingle%d' % k]
                if rnd.random() < 0.93: rq += ['-offormat%d' % k, f]
            dfeat = rnd.random() > 0.08; da = rnd.random() > 0.2
            args = ['a.fa', 'b.fa', '-wordsize', str(ws)] + ([] if da else ['-nodumpalign']) + ([] if dfeat else ['-nodumpfeat']) + rq + (['-outfile', 'o.match'] if da else []) + ['-aoutfeat', names[0], '-boutfeat', names[1], '-logfile', 'log.txt', '-auto']
            for nme in os.listdir(d):
                if nme.endswith(('.out', '.match', '.fa', '.txt', '.gff')) or nme.startswith(('nm', 'ext')): os.remove(os.path.join(d, nme))
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
            if p.returncode != 0 and not read('log.txt'): continue
            written = sorted(x for x in os.listdir(d) if x not in ('a.fa', 'b.fa', 'o.match', 'log.txt', 'faketime.so', 'a.wordmatch'))
            rows.append({'id': len(rows), 'nuc': nuc, 'fasta1': f1, 'fasta2': f2, 'wordsize': ws, 'dumpalign': da, 'dumpfeat': dfeat, 'qual': ' '.join(rq), 'names': names, 'epoch': epoch,
                         'align': read('o.match') or '', 'log': read('log.txt') or '', 'files': {w: read(w) for w in written}, 'error': p.stderr.strip() or None, 'rc': p.returncode})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if r['error']), 'with a message')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261227)
