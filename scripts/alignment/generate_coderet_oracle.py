#!/usr/bin/env python3
"""A5-8c3h (B09 F259, WP58) oracle: EMBOSS 6.6.0 `coderet` (emboss/coderet.c: the coding sequences, messenger RNAs, translations and non-coding regions of sequences with feature tables)
    coderet -auto -seqall in.dat [-sformat1 F] -outfile o.cod [-cdsoutseq [F::]f1 | ''] [-mrnaoutseq ...] [-translationoutseq ...] [-restoutseq ...] [-osformat3 F ... -osformat6 F]
Inputs: EMBL entries made here (CDS / mRNA / exon / gene features with joins, complements, orders, partial ends, /translation tags and notes, several entries per file) and the entries with feature tables of the reader oracle
(generate_seqret_oracle.py: embl, genbank, gff3, gff2, nbrf, swiss, refseqp, FASTA without features, damaged files).
JSONL.GZ rows: id, fmt, how, text, args, epoch, rc, log, files (the sequence files the program wrote), err, hang, memerr
Usage: python3 scripts/alignment/generate_coderet_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0, valgrind, gcc)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402

FMTS = ['embl', 'embl', 'embl', 'embl', 'embl', 'genbank', 'gff3', 'gff2', 'fasta', 'nbrf', 'swiss', 'refseqp', 'embl-own']
OSFORMATS = ['fasta', 'fasta', 'fasta', 'fasta', 'gcg', 'gcg8', 'pir', 'nbrf', 'raw', 'text', 'plain', 'staden', 'ig', 'iguspto', 'codata', 'fitch', 'mase', 'strider', 'gde', 'acedb', 'ncbi', 'gifasta', 'dasdna', 'das', 'asn1', 'nib', 'nibble', 'unknown', 'pearson', 'embl', 'em', 'emblnew', 'emold', 'emblold', 'genbank', 'gb', 'ddbj', 'refseq', 'genpept', 'refseqp', 'swiss', 'sw', 'swissprot', 'uniprot', 'swissnew', 'swold', 'fastq', 'fastq-sanger', 'fastq-illumina', 'fastq-solexa', 'sam', 'gff', 'gff2', 'gff3', 'scf', 'msf', 'clustal', 'aln', 'selex', 'phylip', 'phylipnon', 'phylip3', 'hennig86', 'mega', 'meganon', 'nexus', 'paup', 'nexusnon', 'paupnon', 'jackknifer', 'jackknifernon', 'treecon']   # the writers ported so far (every writer but experiment, debug and bam: B09 WP59)


def wrap_ft(tag, value, quote=True):
    text = '/%s="%s"' % (tag, value) if quote else '/%s=%s' % (tag, value)
    out, first = [], True
    while text:
        out.append('FT                   ' + text[:58]); text = text[58:]
    return out


def loc(rnd, n):
    kind = rnd.random()
    def seg():
        a = rnd.randint(1, n); b = min(n, a + rnd.randint(0, 40)); return a, b
    if kind < 0.35:
        a, b = seg(); s = '%d..%d' % (a, b)
        if rnd.random() < 0.1: s = '<' + s
        if rnd.random() < 0.1: s = s.replace('..', '..>')
        return s
    if kind < 0.45:
        a, b = seg(); return 'complement(%d..%d)' % (a, b)
    if kind < 0.75:
        segs = sorted(seg() for _ in range(rnd.randint(2, 4)))
        s = ','.join('%d..%d' % x for x in segs)
        return rnd.choice(['join(%s)', 'join(%s)', 'order(%s)', 'complement(join(%s))']) % s
    if kind < 0.85:
        segs = sorted(seg() for _ in range(rnd.randint(2, 3)))
        return 'join(%s)' % ','.join('complement(%d..%d)' % x for x in segs)
    if kind < 0.9:
        a, b = seg(); return '%d^%d' % (a, a + 1) if a < n else '%d' % a
    a, _ = seg(); return '%d' % a


def embl_entry(rnd, name):
    n = rnd.randint(20, 400)
    seq = ''.join(rnd.choice('acgt') if rnd.random() > 0.03 else rnd.choice('nrykmswbdhv') for _ in range(n))
    lines = ['ID   %s; SV 1; linear; genomic DNA; STD; UNC; %d BP.' % (name, n), 'XX', 'AC   %s;' % name, 'XX', 'DE   coding test', 'XX', 'FH   Key             Location/Qualifiers', 'FH']
    for _ in range(rnd.randint(0, 7)):
        t = rnd.choice(['CDS', 'CDS', 'CDS', 'mRNA', 'mRNA', 'exon', 'exon', 'gene', 'misc_feature', 'repeat_region', 'intron', 'tRNA', 'CDS'])
        lines.append('FT   %-16s%s' % (t, loc(rnd, n)))
        if t == 'CDS' and rnd.random() < 0.7:
            aa = ''.join(rnd.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(rnd.randint(1, 80)))
            lines += wrap_ft('translation', aa)
        if rnd.random() < 0.3: lines += wrap_ft('gene', rnd.choice(['g1', 'abc', 'x y']))
        if rnd.random() < 0.15: lines += wrap_ft('note', '*translation: MKV' + rnd.choice(['', 'LL']))
        if rnd.random() < 0.15: lines += wrap_ft('codon_start', str(rnd.randint(1, 3)), False)
        if rnd.random() < 0.1: lines += wrap_ft('transl_table', '11', False)
    lines.append('XX')
    lines.append('SQ   Sequence %d BP; %d A; %d C; %d G; %d T; %d other;' % (n, seq.count('a'), seq.count('c'), seq.count('g'), seq.count('t'), n - sum(seq.count(c) for c in 'acgt')))
    for i in range(0, n, 60):
        chunk = seq[i:i + 60]
        groups = ' '.join(chunk[j:j + 10] for j in range(0, len(chunk), 10))
        lines.append('     %-66s%9d' % (groups, min(n, i + 60)))
    lines.append('//')
    return '\n'.join(lines) + '\n'


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        while len(rows) < n:
            fmt = FMTS[len(rows) % len(FMTS)]
            if fmt == 'embl-own':
                text, how = ''.join(embl_entry(rnd, rnd.choice(['X%d' % rnd.randint(1, 99), 'AB%06d' % rnd.randint(1, 999999), 'Seq_%d' % rnd.randint(1, 9)])) for _ in range(rnd.randint(1, 3))), 'own'
                fmt = 'embl'
            else:
                text, how = S.make_input(rnd, gen, fmt)
            if rnd.random() < 0.08: text, how = S.damaged(rnd, text), how + '+damaged'
            args = ['-sformat1', fmt] if rnd.random() < 0.4 else []
            names = ['a.out', 'b.out', 'c.out', 'd.out']
            for opt, nm in zip(['-cdsoutseq', '-mrnaoutseq', '-translationoutseq', '-restoutseq'], names):
                r = rnd.random()
                if r < 0.45: continue
                if r < 0.6: args += [opt, '']
                elif r < 0.85: args += [opt, nm]
                else: args += [opt, rnd.choice(OSFORMATS) + '::' + nm]
            for i in range(3, 7):
                if rnd.random() < 0.3: args += ['-osformat%d' % i, rnd.choice(OSFORMATS)]
            if rnd.random() < 0.04: args += ['-sbegin1', str(rnd.randint(1, 20))]
            if rnd.random() < 0.03: args += ['-sreverse1']
            for nme in os.listdir(work): os.remove(os.path.join(work, nme))
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            cmd = ['-auto', '-seqall', 'in.dat'] + args + ['-outfile', 'o.cod']
            hang = False
            try:
                import resource
                def limit():
                    resource.setrlimit(resource.RLIMIT_FSIZE, (4 << 20, 4 << 20))
                p = subprocess.run(['coderet'] + cmd, cwd=work, capture_output=True, timeout=5, preexec_fn=limit, env=env)
                rc, err = p.returncode, p.stderr.decode('latin-1')
                if rc == -25: hang = True
            except subprocess.TimeoutExpired:
                rc, err, hang = -9, '', True
                subprocess.run(['pkill', '-9', 'coderet'])
            files = {}
            if not hang:
                for nme in sorted(os.listdir(work)):
                    if nme == 'in.dat': continue
                    files[nme] = open(os.path.join(work, nme), encoding='latin-1', newline='').read()
            log = files.pop('o.cod', '')
            memerr = False
            if not hang and shutil.which('valgrind') and (rc != 0 or len(rows) % 4 == 0):
                for nme in os.listdir(work):
                    if nme != 'in.dat': os.remove(os.path.join(work, nme))
                try:
                    vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'coderet'] + cmd, cwd=work, capture_output=True, timeout=120)
                    memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                except subprocess.TimeoutExpired:
                    memerr = True
            if err.startswith('Extract CDS, mRNA and translations from feature tables\n'): err = err[len('Extract CDS, mRNA and translations from feature tables\n'):]
            rows.append({'id': len(rows), 'fmt': fmt, 'how': how, 'text': text, 'args': args, 'epoch': epoch, 'rc': rc, 'log': log, 'files': files, 'err': err, 'hang': hang, 'memerr': memerr})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'died' if r['rc'] else 'files' if r['files'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20265805)
