#!/usr/bin/env python3
"""A5-8c3h (B09 F252) oracle: the sequence WRITERS of EMBOSS 6.6.0 (ajseqwrite.c, `seqOutFormat`) and the output options of seqret
    seqret -auto -sequence in.dat [-sformat1 F] -outseq OUT -osformat2 FMT [-ossingle2 [-osextension2 E] [-osdirectory2 sub]] [-osname2 N] [-firstonly]
(the output files are collected: the sequence file and, with -ossingle, one file per sequence). The input sequences are the inputs of the reader oracle (generate_seqret_oracle.py: writer outputs of seqret, entries with feature
tables, damaged texts) so that names, descriptions, accessions, versions, nucleotide / protein types, gaps and lower case reach the writer in every way the readers give them.
JSONL.GZ rows: id, fmt (the output format), text, spec, args, osformat, outseq, single, osext, osdir, osname, firstonly, epoch, rc, fa (the sequence file), files ({path: text} of -ossingle), err.
Usage: python3 scripts/alignment/generate_seqout_oracle.py OUT.jsonl.gz PER_FORMAT [seed] [format ...]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as G  # noqa: E402

FORMATS = ['unknown', 'fasta', 'pearson', 'gcg', 'gcg8', 'nbrf', 'pir', 'codata', 'strider', 'staden', 'text', 'plain', 'raw', 'ig', 'iguspto', 'gde', 'fitch', 'acedb', 'mase']
# the inputs of the readers whose reference runs are defined (the nexus and mega readers of the library read freed or foreign memory on some inputs, see generate_seqret_oracle.py: memerr)
INPUTS = ['fasta', 'fasta', 'fasta', 'gcg', 'embl', 'genbank', 'nbrf', 'swiss', 'ig', 'codata', 'strider', 'staden', 'gde', 'fitch', 'treecon', 'pir', 'clustal', 'msf', 'phylip', 'raw', 'textonly', 'asis', 'fastq', 'sam']


def collect(work):
    files = {}
    for root, _, names in os.walk(work):
        for n in names:
            p = os.path.join(root, n)
            rel = os.path.relpath(p, work)
            if rel == 'in.dat':
                continue
            files[rel] = open(p, encoding='latin-1', newline='').read()
    return files


def main(out, per, seed, formats):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        for fmt in formats:
            made = 0
            while made < per:
                src = rnd.choice(INPUTS)
                text, how = G.make_input(rnd, gen, src)
                if rnd.random() < 0.15: text, how = G.damaged(rnd, text), how + '+damaged'
                if rnd.random() < 0.05: text = text.rstrip('\n')
                spec, args = 'in.dat', []
                if rnd.random() < 0.7: args = ['-sformat1', src]
                if rnd.random() < 0.08: args += ['-firstonly']
                if rnd.random() < 0.08: args += [rnd.choice(['-slower1', '-supper1', '-sprotein1', '-snucleotide1'])]
                osformat, outseq = fmt, rnd.choice(['o.out', 'out', 'o.txt'])
                if rnd.random() < 0.3:
                    osformat, outseq = None, '%s::%s' % (rnd.choice([fmt, fmt.upper()]), outseq)
                single = rnd.random() < 0.25
                osext = rnd.choice(['abc', 'X1', 'Txt']) if single and rnd.random() < 0.4 else None
                osdir = rnd.choice(['sub', 'sub/']) if single and rnd.random() < 0.4 else None
                osname = rnd.choice(['nm', 'N1']) if rnd.random() < 0.2 else None
                if rnd.random() < 0.03: osformat, outseq = 'bogus', 'o.out'
                for nme in os.listdir(work):
                    p = os.path.join(work, nme)
                    if os.path.isdir(p):
                        for f in os.listdir(p): os.remove(os.path.join(p, f))
                        os.rmdir(p)
                    else:
                        os.remove(p)
                os.mkdir(os.path.join(work, 'sub'))
                open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
                epoch = rnd.randint(978307200, 2051222400)
                env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
                cmd = ['-auto', '-sequence', spec] + args + ['-outseq', outseq]
                if osformat: cmd += ['-osformat2', osformat]
                if single: cmd += ['-ossingle2']
                if osext: cmd += ['-osextension2', osext]
                if osdir: cmd += ['-osdirectory2', osdir]
                if osname: cmd += ['-osname2', osname]
                rc, err, hang = G.seqret(work, cmd, env=env)
                if err.startswith('Read and write (return) sequences\n'): err = err[len('Read and write (return) sequences\n'):]
                files = collect(work) if not hang else {}
                rows.append({'id': len(rows), 'fmt': fmt, 'how': how + '/' + src, 'text': text, 'spec': spec, 'args': args, 'osformat': osformat, 'outseq': outseq, 'single': single, 'osext': osext, 'osdir': osdir, 'osname': osname,
                             'epoch': epoch, 'rc': rc, 'files': files, 'err': err, 'hang': hang})
                made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    c = Counter((r['fmt'], 'files' if r['files'] else 'none') for r in rows)
    print(len(rows), 'rows;', sorted(c.items()))


if __name__ == '__main__':
    fm = sys.argv[4:] or FORMATS
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20265711, fm)
