#!/usr/bin/env python3
"""A5-8c3h (B09 F255, WP58) oracle: the USA forms of the sequence input of EMBOSS 6.6.0 (ajseqread.c `seqinUsaProcess` / `seqinListProcess`, ajquery.c `ajQuerystrParse*`, ajacd.c `acdSetSeqall`) through
    seqret -auto [-feature] -sequence SPEC [-sformat1 F] [-sbegin1 N] [-send1 N] [-sreverse1] [-sid1 ID] [-sdbname1 DB] [-squick1] [-slower1] [-supper1] [-firstonly]
          -outseq o.fa -osformat2 fasta [-offormat2 gff3 -ofname2 o.gff]
SPEC is a file name with a `[begin:end:r]` range, the `:begin:end` / `:begin..end:r` forms of the range, a `format::` prefix, a list file (`@list`, `list::list`, `LIST:list`, nested, with comments, bad records, missing
files, per-entry ranges and formats) or a query on a file (`file:id`, `file:acc:value`, `file:des:value`, wildcards, `file{a|b}`). The files are made from the inputs of the reader oracle
(generate_seqret_oracle.py: writer outputs of seqret, entries with feature tables) so that the range and the reversal reach the sequences AND their features.
JSONL.GZ rows: id, memerr (the reference run has a valgrind error: no defined output), how, files ({name: text}), spec, args, feature, epoch, rc, fa, gff, err, hang.
Usage: python3 scripts/alignment/generate_usa_oracle.py OUT.jsonl.gz COUNT [seed]      (needs EMBOSS 6.6.0)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as G  # noqa: E402

SRC = ['fasta', 'fasta', 'fasta', 'embl', 'genbank', 'genbank', 'swiss', 'nbrf', 'gff3', 'gcg', 'codata', 'ig', 'raw', 'fastq', 'staden', 'phylip', 'clustal', 'msf', 'text']


def multi_fasta(rnd, nuc=None):
    out = ''
    nuc = rnd.random() < 0.8 if nuc is None else nuc
    names = ['a', 'b', 'seq3', 'AB123456', 'x.y', 'Seq5', 'zz']
    rnd.shuffle(names)
    for i in range(rnd.choice([1, 2, 2, 3, 4, 5])):
        n = rnd.choice([3, 8, 20, 33, 61, 100])
        s = ''.join(rnd.choice('ACGT' if nuc else 'ACDEFGHIKLMNPQRSTVWY') for _ in range(n))
        acc = rnd.choice(['', '', 'XX%06d' % rnd.randint(1, 999999), 'P12345'])
        desc = rnd.choice(['', 'desc', 'some kinase', 'alpha beta'])
        out += '>%s%s%s\n%s\n' % (names[i], (' ' + acc) if acc else '', (' ' + desc) if desc else '', ''.join(s[j:j + 60] + '\n' for j in range(0, len(s), 60)).rstrip('\n'))
    return out


def one_input(rnd, gen):
    src = rnd.choice(SRC)
    if src == 'fasta' and rnd.random() < 0.6:
        return multi_fasta(rnd), src, 'multi'
    text, how = G.make_input(rnd, gen, src)
    if rnd.random() < 0.1:
        text, how = G.damaged(rnd, text), how + '+damaged'
    return text, src, how


def rng(rnd):
    """a [begin:end:r] range text"""
    b = rnd.choice(['', '', '1', '3', '5', '10', '-4', '-1', '0', '25', '70', '150'])
    e = rnd.choice(['', '', '8', '20', '2', '-3', '-1', '0', '60', '120', '400'])
    rev = rnd.choice(['', '', ':r', ':R', ':', ':x', ':rr'])
    if b == '' and e == '' and rnd.random() < 0.5:
        rev = rnd.choice([':r', ':R'])
    return '[%s:%s%s]' % (b, e, rev)


def usa(rnd, fname, fmt, allow_range=True, queries=True):
    """one USA for the file"""
    r = rnd.random()
    s = fname
    if r < 0.40 and allow_range:
        s += rng(rnd)
    elif r < 0.46 and allow_range:
        s += rnd.choice([':3:10', ':3..10', ':3..10:r', ':5:2', ':1:1:r', ':0:5', ':-3:5', ':4..', ':..9', ':3:10:R'])
    elif r < 0.62 and queries:
        s += rnd.choice([':a', ':b', ':a*', ':*', ':?', ':seq3', ':AB123456', ':id:a', ':ID:b', ':acc:XX*', ':acc:P12345', ':des:desc', ':des:*kinase*', ':sv:x', ':org:Homo*', ':key:a', ':z', ':a|b', ':a,b', ':a b',
                         ':{a|b}', '{a}', '{a|b|zz}', '{id:a}', ':id:{a|b}', ':id:a|acc:P12345', ':id:a&des:desc', ':id:a!des:desc', ':nosuch:a', ':a*:', '{a', 'a}', ':%d' % rnd.randint(0, 40), '%0', '%10', '%3'])
        if rnd.random() < 0.3 and allow_range and '[' not in s:
            s += rng(rnd)
    if rnd.random() < 0.18:
        s = rnd.choice([fmt, fmt.upper(), '', 'fasta', 'embl', 'unknownfmt']) + '::' + s
    return s


def main(out, count, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        made = 0
        while made < count:
            files = {}
            text, fmt, how = one_input(rnd, gen)
            files['in.dat'] = text
            text2, fmt2, how2 = one_input(rnd, gen) if rnd.random() < 0.6 else (multi_fasta(rnd), 'fasta', 'multi')
            files['two.dat'] = text2
            fmts = {'in.dat': fmt, 'two.dat': fmt2}
            kind = rnd.random()
            args = []
            spec = ''
            if kind < 0.30:
                spec = usa(rnd, 'in.dat', fmt)
            elif kind < 0.70:
                # list files
                lines = []
                for _ in range(rnd.choice([1, 2, 2, 3, 4])):
                    t = rnd.random()
                    if t < 0.7:
                        f = rnd.choice(['in.dat', 'two.dat'])
                        lines.append(usa(rnd, f, fmts[f], True, rnd.random() < 0.5))
                    elif t < 0.78:
                        lines.append(rnd.choice(['missing.dat', 'in.dat extra', 'two.dat # comment', '# only comment', '', '   ', '\t', 'in.dat#x', 'l2.lst', '@l2.lst', 'list::l2.lst', '@l1.lst', '@missing.lst', 'asis::ACGT', 'fasta::two.dat[2:4]']))
                    else:
                        lines.append(rnd.choice(['in.dat', 'two.dat']) + rnd.choice(['', '  # c', '\t#c']))
                files['l1.lst'] = '\n'.join(lines) + rnd.choice(['\n', '', '\n\n'])
                files['l2.lst'] = '\n'.join(usa(rnd, rnd.choice(['in.dat', 'two.dat']), 'fasta', True, False) for _ in range(rnd.choice([1, 2]))) + '\n'
                prefix = rnd.choice(['@', '@', 'list::', 'LIST:', 'List::', '@ ', 'fasta::@', '@@'])
                spec = prefix + 'l1.lst'
                if rnd.random() < 0.25:
                    spec += rng(rnd)
            elif kind < 0.78:
                spec = rnd.choice(['missing.dat', 'in.dat[', 'in.dat]', '[1:2]', 'in.dat[1:2][3:4]', '  in.dat', 'in.dat  ', 'in.dat[ 1:2]', 'in.dat[1: 2]', 'in.dat[1:2]x', '@', 'list::', 'asis::', 'in.dat:', 'in.dat::', '::in.dat', 
                                   'dat:in.dat', 'nosuchdb:a', 'nosuch:id:a', 'in.dat%5', 'in.dat%x'])
            else:
                spec = rnd.choice(['in.dat', 'two.dat'])
            if rnd.random() < 0.25:
                args += ['-sformat1', fmt]
            if rnd.random() < 0.30:
                args += ['-sbegin1', rnd.choice(['2', '5', '-6', '30', '1', '0'])]
            if rnd.random() < 0.30:
                args += ['-send1', rnd.choice(['8', '-3', '40', '3', '0'])]
            if rnd.random() < 0.22:
                args += ['-sreverse1']
            if rnd.random() < 0.10:
                args += ['-sid1', rnd.choice(['a', 'b', 'zz', 'idx', 'My_id'])]
            if rnd.random() < 0.05:
                args += ['-sdbname1', rnd.choice(['mydb', 'x'])]
            if rnd.random() < 0.05:
                args += ['-squick1']
            if rnd.random() < 0.05:
                args += [rnd.choice(['-slower1', '-supper1', '-sprotein1', '-snucleotide1', '-scircular1'])]
            if rnd.random() < 0.05:
                args += ['-iquery1', rnd.choice(['a', 'a*', 'id:a', '{a|b}'])]
            if rnd.random() < 0.04:
                args += ['-ioffset1', rnd.choice(['0', '5', '20', '60'])]
            if rnd.random() < 0.06:
                args += ['-firstonly']
            feature = rnd.random() < 0.7
            for nme in os.listdir(work):
                os.remove(os.path.join(work, nme))
            for nme, t in files.items():
                open(os.path.join(work, nme), 'w', newline='\n', encoding='latin-1').write(t)
            cmd = ['-auto'] + (['-feature'] if feature else []) + ['-sequence', spec] + args + ['-outseq', 'o.fa', '-osformat2', 'fasta'] + (['-offormat2', 'gff3', '-ofname2', 'o.gff'] if feature else [])
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            rc, err, hang = G.seqret(work, cmd, env=env)
            def rd(nme):
                q = os.path.join(work, nme)
                return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
            if 'ajSysCreateInpipe' in err or 'ajSysCreate' in err:
                continue                # a file name that ends with '|' is a command the library runs: not part of the oracle
            if err.startswith('Read and write (return) sequences\n'):
                err = err[len('Read and write (return) sequences\n'):]
            memerr = False
            if not hang and shutil.which('valgrind') and (fmt in G.VG_FORMATS or fmt2 in G.VG_FORMATS):
                # a reference run that reads or writes memory it does not own has no defined output
                vgargs = [a.replace('o.fa', 'ovg.fa').replace('o.gff', 'ovg.gff') for a in cmd]
                try:
                    vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'seqret'] + vgargs, cwd=work, capture_output=True, timeout=120, env=dict(os.environ))
                    memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                except subprocess.TimeoutExpired:
                    memerr = True
            rows.append({'id': len(rows), 'memerr': memerr, 'how': '%s/%s/%s' % (fmt, how, fmt2), 'files': files, 'spec': spec, 'args': args, 'feature': feature, 'epoch': epoch, 'rc': rc, 'fa': rd('o.fa'), 'gff': rd('o.gff'), 'err': err, 'hang': hang})
            made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    c = Counter('hang' if r['hang'] else 'died' if r['rc'] else 'fa' if r['fa'] else 'none' for r in rows)
    print(len(rows), 'rows;', sorted(c.items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20265801)
