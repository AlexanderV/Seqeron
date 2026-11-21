#!/usr/bin/env python3
"""A5-8c3c (B09 F220) oracle: EMBOSS 6.6.0 `diffseq` with the FEATURE OUTPUT qualifiers of both feature files (-offormat4/5 of every output feature format and its aliases, other case, an unknown name, the `format::file` form of
-aoutfeat / -boutfeat, -ofopenfile, -ofname / -ofextension / -ofsingle) on nucleotide and protein pairs, with the sequence qualifiers (-sbegin / -send / -sreverse ...): the report and both feature files byte for byte, the message of the
program when a format is refused. Same row schema as generate_diffseq_qualifiers_oracle.py plus `name1` / `name2` (the files named on the command line).
Usage: python3 scripts/alignment/generate_diffseq_featout_oracle.py OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_diffseq_oracle as D  # noqa: E402
import generate_diffseq_qualifiers_oracle as Q  # noqa: E402

FORMATS = ['gff', 'gff2', 'gff3', 'bed', 'embl', 'em', 'genbank', 'gb', 'ddbj', 'refseq', 'refseqp', 'pir', 'nbrf', 'swiss', 'sw', 'swissprot', 'uniprot', 'dasgff', 'draw', 'debug', 'unknown']
ODD = ['GFF3', 'Embl', 'GENBANK', 'bogus', 'fasta', 'gff4']


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            nuc, a, b = D.case(rnd)
            if not a or not b: continue
            if nuc and rnd.random() < 0.2:                       # letters the tag validator of /replace refuses (the program writes a note instead)
                b = ''.join(c if rnd.random() > 0.03 else rnd.choice('ELQPZJFIO') for c in b)
            n1 = rnd.choice(['seqa', 'X56789.2', 'gb_a', 'Q12345', 'a"b', 'long_sequence_name_that_needs_wrapping_in_the_note_of_a_feature_xx']); n2 = rnd.choice(['seqb', 'sample2', 'chr7', 'P99999', 'it\'s', 'x y'])
            f1 = D.fasta(n1, rnd.choice(['', 'first sequence', 'OS=Homo sapiens']), a); f2 = D.fasta(n2, rnd.choice(['', 'second one']), b)
            ws = rnd.choice([2, 3, 4, 5, 6, 8, 10, 10, 12]); glob = rnd.random() < 0.4
            rq = (Q.seq_qualifiers(rnd, 1, nuc, len(a)) if rnd.random() < 0.35 else []) + (Q.seq_qualifiers(rnd, 2, nuc, len(b)) if rnd.random() < 0.35 else [])
            if rnd.random() < 0.25: rq += ['-rformat', rnd.choice(['table', 'gff', 'embl', 'genbank', 'dasgff', 'pir', 'swiss', 'simple', 'tagseq', 'debug'])]
            names = ['a.out', 'b.out']
            for k in (4, 5):
                r = rnd.random()
                f = rnd.choice(FORMATS) if r < 0.8 else rnd.choice(ODD)
                if nuc and f in ('pir', 'nbrf', 'swiss', 'sw', 'swissprot', 'uniprot', 'refseqp') and rnd.random() < 0.6: f = rnd.choice(['embl', 'genbank', 'gff2', 'bed', 'dasgff'])
                if not nuc and f in ('embl', 'em', 'genbank', 'gb', 'ddbj', 'refseq') and rnd.random() < 0.6: f = rnd.choice(['pir', 'swiss', 'refseqp', 'gff2', 'bed', 'dasgff', 'draw'])
                r = rnd.random()
                if r < 0.07: names[k - 4] = '%s::%s' % (rnd.choice(FORMATS), names[k - 4])
                elif r < 0.1: rq += ['-ofopenfile%d' % k, 'f%d.out' % k]
                elif r < 0.13: rq += ['-ofname%d' % k, 'nm']
                elif r < 0.16: rq += ['-ofextension%d' % k, 'ext']
                elif r < 0.19: rq += ['-ofsingle%d' % k]
                if rnd.random() < 0.93: rq += ['-offormat%d' % k, f]
            args = ['a.fa', 'b.fa', '-wordsize', str(ws)] + (['-globaldifferences'] if glob else []) + rq + ['-outfile', 'o.rep', '-aoutfeat', names[0], '-boutfeat', names[1], '-auto']
            for nme in os.listdir(d):
                if nme.endswith(('.out', '.rep', '.fa', '.diffgff')) or nme.startswith(('nm', 'ext')): os.remove(os.path.join(d, nme))
            open(os.path.join(d, 'a.fa'), 'w', newline='\n').write(f1); open(os.path.join(d, 'b.fa'), 'w', newline='\n').write(f2)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            try:
                p = subprocess.run(['diffseq'] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=20)
            except subprocess.TimeoutExpired:
                continue
            def read(nme):
                q = os.path.join(d, nme)
                return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) else None
            if p.returncode != 0 and not read('o.rep'): continue
            # the files the program wrote (the name on the command line, after the `format::` prefix, or the -ofopenfile name)
            written = sorted(x for x in os.listdir(d) if x not in ('a.fa', 'b.fa', 'o.rep', 'faketime.so'))
            entries = ['[-asequence] a.fa', '[-bsequence] b.fa', '-wordsize %d' % ws] + (['-globaldifferences'] if glob else [])
            i = 0
            while i < len(rq):
                if i + 1 < len(rq) and not rq[i + 1].startswith('-') or (i + 1 < len(rq) and rq[i + 1].lstrip('-').isdigit()): entries.append(rq[i] + ' ' + rq[i + 1]); i += 2
                else: entries.append(rq[i]); i += 1
            entries += ['-outfile o.rep', '-aoutfeat ' + names[0], '-boutfeat ' + names[1], '-auto']
            rows.append({'id': len(rows), 'nuc': nuc, 'fasta1': f1, 'fasta2': f2, 'wordsize': ws, 'glob': glob, 'rq': ' '.join(rq), 'names': names, 'cmd': entries, 'epoch': epoch,
                         'report': read('o.rep'), 'files': {w: read(w) for w in written}, 'error': p.stderr.strip() or None, 'rc': p.returncode})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if r['error']), 'with a message')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261226)
