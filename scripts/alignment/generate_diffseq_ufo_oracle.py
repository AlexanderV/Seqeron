#!/usr/bin/env python3
"""A5-8c3e (B09 F241) oracle: EMBOSS 6.6.0 `diffseq` with the features of its sequences read from a SEPARATE file (-ufo1 / -ufo2, -fformat, -fopenfile; the UFO is `[format::]file`): the first_feature / second_feature notes
of the report and both GFF3 files, for feature files of every input format (generate_featread_oracle.py), with and without the format, a file that does not exist, a format that is none, sequences that are FASTA, GenBank or
EMBL (features of the sequence file are not read when a UFO is given, and are when none is), -sbegin / -send / -sreverse (the table follows the sequence).
JSONL.GZ rows as generate_diffseq_features_oracle.py plus ffiles (the feature files by name).
Usage: python3 scripts/alignment/generate_diffseq_ufo_oracle.py OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_diffseq_oracle as D  # noqa: E402
import generate_diffseq_features_oracle as D2  # noqa: E402
import generate_featread_oracle as FR  # noqa: E402


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        while len(rows) < n:
            a = ''
            while len(a) < 60 or not all(c in 'ACGT' for c in a):
                a = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(80, 300)))
            b = D.mutate(rnd, a, 'ACGT', rnd.randint(1, 6))
            n1 = rnd.choice(['SEQA', 'AB123456', 'gb_a']); n2 = rnd.choice(['SEQB', 'chr7', 'Y99999'])
            fmt1 = rnd.choice(['fasta', 'fasta', 'genbank', 'embl']); fmt2 = rnd.choice(['fasta', 'fasta', 'genbank', 'embl'])
            def seqfile(fmt, name, seq):
                if fmt == 'fasta': return D.fasta(name, 'test', seq)
                return D2.entry(rnd, fmt, name, 'test', seq, D2.feature_lines(rnd, fmt, len(seq)))
            f1 = seqfile(fmt1, n1, a); f2 = seqfile(fmt2, n2, b)
            ext = {'genbank': 'gb', 'embl': 'embl', 'fasta': 'fa'}
            fa, fb = 'a.' + ext[fmt1], 'b.' + ext[fmt2]
            ffiles = {}
            rq = []
            for k, name in ((1, n1), (2, n2)):
                r = rnd.random()
                kind = rnd.choice(FR.KINDS)
                text = FR.make(rnd, kind)
                if rnd.random() < 0.1: text = FR.damaged(rnd, text)
                fname = 'f%d.dat' % k
                if r < 0.45:
                    ffiles[fname] = text; rq += ['-ufo%d' % k, fname]
                elif r < 0.7:
                    ffiles[fname] = text; rq += ['-ufo%d' % k, '%s::%s' % (rnd.choice([FR.FMTNAME[kind], FR.FMTNAME[kind], 'gff3', 'embl', 'bogus', 'GFF3', 'swiss']), fname)]
                elif r < 0.8:
                    ffiles[fname] = text; rq += ['-ufo%d' % k, fname, '-fformat%d' % k, rnd.choice([FR.FMTNAME[kind], 'gff3'])]
                elif r < 0.88:
                    ffiles['g%d.dat' % k] = text; rq += ['-fformat%d' % k, FR.FMTNAME[kind], '-fopenfile%d' % k, 'g%d.dat' % k]
                elif r < 0.93:
                    ffiles['%s.%s' % (name, FR.FMTNAME[kind])] = text; rq += ['-fformat%d' % k, FR.FMTNAME[kind]]
                elif r < 0.96:
                    rq += ['-ufo%d' % k, 'missing.dat']
            if rnd.random() < 0.3:
                for k, length in ((1, len(a)), (2, len(b))):
                    if rnd.random() < 0.6: rq += ['-sbegin%d' % k, str(rnd.randint(1, max(1, length // 3)))]
                    if rnd.random() < 0.6: rq += ['-send%d' % k, str(rnd.randint(max(1, length // 2), length))]
                    if rnd.random() < 0.3: rq.append('-sreverse%d' % k)
            rq += ['-rformat', rnd.choice(['diffseq'] * 6 + ['table', 'gff', 'tagseq', 'simple'])]
            ws = rnd.choice([3, 4, 5, 6, 8]); glob = rnd.random() < 0.3
            args = [fa, fb, '-wordsize', str(ws)] + (['-globaldifferences'] if glob else []) + rq + ['-outfile', 'o.rep', '-aoutfeat', 'a.gff', '-boutfeat', 'b.gff', '-auto']
            for nme in os.listdir(work): os.remove(os.path.join(work, nme))
            open(os.path.join(work, fa), 'w', newline='\n').write(f1); open(os.path.join(work, fb), 'w', newline='\n').write(f2)
            for nme, text in ffiles.items(): open(os.path.join(work, nme), 'w', newline='\n', encoding='latin-1').write(text)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            try:
                p = subprocess.run(['diffseq'] + args, cwd=work, capture_output=True, errors='replace', env=env, timeout=8)
            except subprocess.TimeoutExpired:
                subprocess.run(['pkill', '-9', 'diffseq']); continue
            def read(nme):
                q = os.path.join(work, nme)
                return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) else ''
            entries = ['[-asequence] ' + fa, '[-bsequence] ' + fb, '-wordsize %d' % ws] + (['-globaldifferences'] if glob else [])
            i = 0
            while i < len(rq):
                if i + 1 < len(rq) and (not rq[i + 1].startswith('-') or rq[i + 1].lstrip('-').isdigit()): entries.append(rq[i] + ' ' + rq[i + 1]); i += 2
                else: entries.append(rq[i]); i += 1
            entries += ['-outfile o.rep', '-aoutfeat a.gff', '-boutfeat b.gff', '-auto']
            rows.append({'id': len(rows), 'file1': f1, 'file2': f2, 'name1': fa, 'name2': fb, 'wordsize': ws, 'glob': glob, 'rq': ' '.join(rq), 'cmd': entries, 'epoch': epoch, 'ffiles': ffiles,
                         'report': read('o.rep'), 'feat1': read('a.gff'), 'feat2': read('b.gff'), 'error': p.stderr.strip() or None, 'rc': p.returncode})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows', sum(1 for r in rows if 'first_feature' in r['report'] or 'second_feature' in r['report']), 'with feature notes')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 300, int(sys.argv[3]) if len(sys.argv) > 3 else 20261502)
