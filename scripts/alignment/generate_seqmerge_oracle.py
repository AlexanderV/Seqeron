#!/usr/bin/env python3
"""A5-8c3h (B09 F253, F254) oracle: EMBOSS 6.6.0 `union` (emboss/union.c: the sequences of a file joined into one, with the overlaps of consecutive sequences found and counted) and `seqretsplit` (emboss/seqretsplit.c:
seqret that writes one file per sequence and renames duplicate names)
    union -auto -sequence in.dat [-sformat1 F] -outseq OUT [-osformat2 FMT] [-findoverlap] [-overlapfile ov.txt] [-firstonly ...]
    seqretsplit -auto -sequence in.dat [-sformat1 F] -outseq OUT [-osformat2 FMT] [-ossingle2 | -noossingle2] [-osextension2 E] [-osdirectory2 sub] [-firstonly]
The inputs are fragments of one sequence that overlap by random lengths (union) or sequences with repeated names in any case (seqretsplit), written in the formats of the reader oracle, plus the reader oracle texts.
JSONL.GZ rows: id, prog, text, spec, args, osformat, outseq, findoverlap, overlapfile, single (null / true / false), osext, osdir, epoch, rc, files ({path: text}), err, hang.
Usage: python3 scripts/alignment/generate_seqmerge_oracle.py OUT.jsonl.gz PER_PROGRAM [seed]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as G  # noqa: E402
import generate_seqout_oracle as O  # noqa: E402

FMTS = ['fasta', 'gcg', 'nbrf', 'codata', 'strider', 'staden', 'text', 'ig', 'gde', 'fitch', 'acedb', 'pir', 'unknown', 'pearson']


def fragments(rnd):
    nuc = rnd.random() < 0.8
    alpha = G.NUC if nuc else G.PROT
    base = ''.join(rnd.choice(alpha) for _ in range(rnd.choice([30, 60, 120, 200])))
    out, pos = [], 0
    names = ['frag1', 'frag2', 'Frag3', 'ctg.4', 'seq_5', 'x6']
    for i in range(rnd.choice([1, 2, 2, 3, 4, 5])):
        ln = rnd.randint(5, 60)
        seq = base[pos:pos + ln] or base[-ln:]
        out.append((names[i % len(names)] if rnd.random() < 0.9 else '', seq))
        pos += max(1, ln - rnd.choice([0, 0, 1, 3, 5, 10, ln]))
        if pos >= len(base): pos = max(0, len(base) - 10)
    if rnd.random() < 0.1:
        out.append(('dup', out[0][1]))
    fa = ''.join('>%s %s\n%s\n' % (n, rnd.choice(['', 'part', 'a b c']), s) for n, s in out)
    return fa


def dupnames(rnd):
    names = ['a', 'A', 'b', 'a', 'B', 'a.002', 'a', 'c', 'C', 'a.003', 'a']
    k = rnd.randint(1, 9)
    return ''.join('>%s\n%s\n' % (rnd.choice(names) if rnd.random() < 0.8 else ('' if rnd.random() < 0.5 else 'a'), ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(3, 40)))) for _ in range(k))


def main(out, per, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        for prog in ('union', 'seqretsplit'):
            made = 0
            while made < per:
                r = rnd.random()
                if r < 0.55:
                    text, how = (fragments(rnd) if prog == 'union' else dupnames(rnd)), 'fragments'
                else:
                    src = rnd.choice(O.INPUTS)
                    text, how = G.make_input(rnd, gen, src); how += '/' + src
                if rnd.random() < 0.15: text, how = G.damaged(rnd, text), how + '+damaged'
                args = []
                if rnd.random() < 0.5 and how != 'fragments': args = ['-sformat1', how.split('/')[-1].split('+')[0]]
                if rnd.random() < 0.06: args += [rnd.choice(['-slower1', '-supper1', '-sprotein1', '-snucleotide1'])]
                osformat = rnd.choice(FMTS) if rnd.random() < 0.5 else None
                outseq = rnd.choice(['o.out', 'out'])
                findoverlap = prog == 'union' and rnd.random() < 0.6
                overlapfile = 'ov.txt' if findoverlap and rnd.random() < 0.7 else (('ov.txt') if prog == 'union' and rnd.random() < 0.1 else None)
                single = None
                osext = osdir = None
                firstonly = prog == 'seqretsplit' and rnd.random() < 0.08
                if prog == 'seqretsplit':
                    single = rnd.choice([None, None, None, True, False])
                    if rnd.random() < 0.3: osext = rnd.choice(['abc', 'X1'])
                    if rnd.random() < 0.3: osdir = 'sub'
                for nme in os.listdir(work):
                    p = os.path.join(work, nme)
                    if os.path.isdir(p):
                        for f in os.listdir(p): os.remove(os.path.join(p, f))
                        os.rmdir(p)
                    else: os.remove(p)
                os.mkdir(os.path.join(work, 'sub'))
                open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
                epoch = rnd.randint(978307200, 2051222400)
                env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
                cmd = ['-auto', '-sequence', 'in.dat'] + args + ['-outseq', outseq]
                if osformat: cmd += ['-osformat2', osformat]
                if findoverlap: cmd += ['-findoverlap']
                if overlapfile: cmd += ['-overlapfile', overlapfile]
                if firstonly: cmd += ['-firstonly']
                if single is True: cmd += ['-ossingle2']
                elif single is False: cmd += ['-noossingle2']
                if osext: cmd += ['-osextension2', osext]
                if osdir: cmd += ['-osdirectory2', osdir]
                exe = [prog] + cmd
                hang = False
                proc = subprocess.Popen(exe, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
                try:
                    _, errb = proc.communicate(timeout=5)
                    rc, err = proc.returncode, errb.decode('latin-1')
                except subprocess.TimeoutExpired:
                    proc.kill(); proc.communicate(); rc, err, hang = -9, '', True
                if err.startswith('Concatenate multiple sequences into a single sequence\n'): err = err[len('Concatenate multiple sequences into a single sequence\n'):]
                if err.startswith('Read sequences and write them to individual files\n'): err = err[len('Read sequences and write them to individual files\n'):]
                files = O.collect(work) if not hang else {}
                rows.append({'id': len(rows), 'prog': prog, 'how': how, 'text': text, 'spec': 'in.dat', 'args': args, 'osformat': osformat, 'outseq': outseq, 'findoverlap': findoverlap, 'overlapfile': overlapfile, 'single': single,
                             'osext': osext, 'osdir': osdir, 'firstonly': firstonly, 'epoch': epoch, 'rc': rc, 'files': files, 'err': err, 'hang': hang})
                made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['prog'], 'files' if r['files'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20265712)
