#!/usr/bin/env python3
"""A5-8c3b (B09 F216) oracle: EMBOSS 6.6.0 `diffseq` with the SEQUENCE qualifiers of both inputs (-sbegin1/2 -send1/2 -sreverse1/2 -slower -supper -snucleotide -sprotein, -sformat, -rusashow) on FASTA inputs:
the report (7-tag `diffseq` format and other -rformat) and both GFF3 feature files, byte for byte. Same row schema as generate_diffseq_oracle.py (fasta1 / fasta2 are the input files; `rq` the qualifier text after the two files).
Usage: python3 scripts/alignment/generate_diffseq_qualifiers_oracle.py OUT.jsonl.gz [n=260] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_diffseq_oracle as D  # noqa: E402
import generate_emboss_report_formats_oracle as G  # noqa: E402


def seq_qualifiers(rnd, n, nuc, length):
    """the qualifier words of input n (1 / 2)"""
    q = []
    if rnd.random() < 0.7:
        b = rnd.choice([rnd.randint(1, max(1, length // 3)), -rnd.randint(1, max(1, length // 2)), rnd.randint(1, length + 5)])
        q += ['-sbegin%d' % n, str(b)]
    if rnd.random() < 0.7:
        e = rnd.choice([rnd.randint(max(1, length // 2), length), -rnd.randint(1, max(1, length // 3)), rnd.randint(1, length + 5), rnd.randint(length // 2, length + 30)])
        q += ['-send%d' % n, str(e)]
    if rnd.random() < 0.3: q.append('-sreverse%d' % n)
    r = rnd.random()
    if r < 0.12: q.append('-slower%d' % n)
    elif r < 0.24: q.append('-supper%d' % n)
    if rnd.random() < 0.05: q.append('-snucleotide%d' % n if nuc else '-sprotein%d' % n)
    if rnd.random() < 0.1: q += ['-sformat%d' % n, 'fasta']
    return q


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            nuc, a, b = D.case(rnd)
            if not a or not b: continue
            n1 = rnd.choice(['seqa', 'X56789.2', 'gb_a', 'Q12345']); n2 = rnd.choice(['seqb', 'sample2', 'chr7', 'P99999'])
            f1 = D.fasta(n1, rnd.choice(['', 'first sequence', 'OS=Homo sapiens']), a); f2 = D.fasta(n2, rnd.choice(['', 'second one']), b)
            ws = rnd.choice([2, 3, 4, 5, 6, 8, 10, 10, 12]); glob = rnd.random() < 0.4
            rq = (seq_qualifiers(rnd, 1, nuc, len(a)) if rnd.random() < 0.8 else []) + (seq_qualifiers(rnd, 2, nuc, len(b)) if rnd.random() < 0.8 else [])
            rq += ['-rformat', rnd.choice(['diffseq', 'diffseq', 'diffseq', 'table', 'gff', 'embl', 'excel', 'simple', 'tagseq', 'srs', 'listfile', 'seqtable'])]
            if rnd.random() < 0.15: rq.append('-rusashow')
            if rnd.random() < 0.1: rq.append('-raccshow')
            args = ['a.fa', 'b.fa', '-wordsize', str(ws)] + (['-globaldifferences'] if glob else []) + rq + ['-outfile', 'o.rep', '-aoutfeat', 'a.gff', '-boutfeat', 'b.gff', '-auto']
            for nme in ('o.rep', 'a.gff', 'b.gff'):
                if os.path.exists(os.path.join(d, nme)): os.remove(os.path.join(d, nme))
            open(os.path.join(d, 'a.fa'), 'w', newline='\n').write(f1); open(os.path.join(d, 'b.fa'), 'w', newline='\n').write(f2)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            try:
                p = subprocess.run(['diffseq'] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=20)
            except subprocess.TimeoutExpired:
                continue
            def read(nme):
                q = os.path.join(d, nme)
                return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) else ''
            entries = ['[-asequence] a.fa', '[-bsequence] b.fa', '-wordsize %d' % ws] + (['-globaldifferences'] if glob else [])
            i = 0
            while i < len(rq):
                if i + 1 < len(rq) and not rq[i + 1].startswith('-') or (i + 1 < len(rq) and rq[i + 1].lstrip('-').isdigit()): entries.append(rq[i] + ' ' + rq[i + 1]); i += 2
                else: entries.append(rq[i]); i += 1
            entries += ['-outfile o.rep', '-aoutfeat a.gff', '-boutfeat b.gff', '-auto']
            rows.append({'id': len(rows), 'nuc': nuc, 'fasta1': f1, 'fasta2': f2, 'wordsize': ws, 'glob': glob, 'rq': ' '.join(rq), 'cmd': entries, 'epoch': epoch,
                         'report': read('o.rep'), 'feat1': read('a.gff'), 'feat2': read('b.gff'), 'error': p.stderr.strip() or None, 'rc': p.returncode})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if r['error']), 'with a message')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 260, int(sys.argv[3]) if len(sys.argv) > 3 else 20261222)
