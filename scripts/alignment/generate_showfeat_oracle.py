#!/usr/bin/env python3
"""A5-8c3h (B09 F260, WP58) oracle: EMBOSS 6.6.0 `showfeat` (emboss/showfeat.c: the features of the sequences as lines of ASCII graphics)
    showfeat -auto -sequence in.dat [-sformat1 F] -outfile o.sf [-sourcematch S] [-typematch T] [-tagmatch T] [-valuematch V] [-sort source|start|type|nosort] [-joinfeatures] [-html] [-noid] [-nodescription] [-noscale]
             [-width N] [-collapse] [-noforward] [-noreverse] [-nounknown] [-strand] [-sourceshow] [-position] [-notypeshow] [-tagshow] [-novalueshow] [-stricttags] [-annotation RANGES] [-sbegin1 N -send1 N -sreverse1]
Inputs: entries with feature tables in the formats whose readers carry the features (embl, genbank, swiss, nbrf, gff3, refseqp, gff2) made as in generate_seqret_oracle.py, entries with CDS / mRNA / exon features made as in
generate_coderet_oracle.py, FASTA files without features, damaged files.
JSONL.GZ rows: id, fmt, how, text, spec, args, rc, out, err, hang, memerr
Usage: python3 scripts/alignment/generate_showfeat_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0, valgrind)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402
import generate_coderet_oracle as C  # noqa: E402
import generate_extractfeat_oracle as E  # noqa: E402

FMTS = ['embl', 'genbank', 'swiss', 'nbrf', 'gff3', 'refseqp', 'gff2', 'embl', 'genbank', 'gff3', 'swiss', 'embl-own', 'embl-own', 'embl-own', 'genbank', 'fasta']
WIDTHS = [60, 60, 60, 10, 30, 80, 120, 4, 3, 1, 0, 2, 200, 59, 61, -5, 7]


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        while len(rows) < n:
            fmt = FMTS[len(rows) % len(FMTS)]
            if fmt == 'embl-own':
                text, how = ''.join(C.embl_entry(rnd, rnd.choice(['X%d' % rnd.randint(1, 99), 'AB%06d' % rnd.randint(1, 999999), 'Seq_%d' % rnd.randint(1, 9)])) for _ in range(rnd.randint(1, 3))), 'own'
                fmt = 'embl'
            else:
                text, how = S.make_input(rnd, gen, fmt)
            if rnd.random() < 0.1: text, how = S.damaged(rnd, text), how + '+damaged'
            spec = 'in.dat'
            args = ['-sformat1', fmt] if rnd.random() < 0.4 else []
            if rnd.random() < 0.3: args += ['-typematch', rnd.choice(E.TYPES)]
            if rnd.random() < 0.08: args += ['-sourcematch', rnd.choice(E.SOURCES)]
            if rnd.random() < 0.15: args += ['-tagmatch', rnd.choice(E.TAGS)]
            if rnd.random() < 0.1: args += ['-valuematch', rnd.choice(E.VALUES)]
            if rnd.random() < 0.5: args += ['-sort', rnd.choice(['source', 'start', 'type', 'nosort'])]
            if rnd.random() < 0.25: args.append('-joinfeatures')
            if rnd.random() < 0.08: args.append('-html')
            if rnd.random() < 0.1: args.append('-noid')
            if rnd.random() < 0.1: args.append('-nodescription')
            if rnd.random() < 0.1: args.append('-noscale')
            if rnd.random() < 0.4: args += ['-width', str(rnd.choice(WIDTHS))]
            if rnd.random() < 0.1: args.append('-collapse')
            if rnd.random() < 0.08: args.append('-noforward')
            if rnd.random() < 0.08: args.append('-noreverse')
            if rnd.random() < 0.05: args.append('-nounknown')
            if rnd.random() < 0.3: args.append('-strand')
            if rnd.random() < 0.3: args.append('-sourceshow')
            if rnd.random() < 0.3: args.append('-position')
            if rnd.random() < 0.1: args.append('-notypeshow')
            if rnd.random() < 0.35: args.append('-tagshow')
            if rnd.random() < 0.1: args.append('-novalueshow')
            if rnd.random() < 0.1: args.append('-stricttags')
            if rnd.random() < 0.15: args += ['-annotation', rnd.choice(['5-10 foo', '3-9 bar,20-30 baz', '1-100', '10-20 x y', '50-60,70-80', '5 10 a', '1-5', '2-3 a b c', '300-400 far', '8-12 a,8-12 a'])]
            if rnd.random() < 0.08: args += ['-sbegin1', str(rnd.randint(1, 30))]
            if rnd.random() < 0.08: args += ['-send1', str(rnd.randint(20, 120))]
            if rnd.random() < 0.04: args += ['-sreverse1']
            for nme in os.listdir(work): os.remove(os.path.join(work, nme))
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            cmd = ['-auto', '-sequence', spec] + args + ['-outfile', 'o.sf']
            hang = False
            try:
                import resource
                def limit():
                    resource.setrlimit(resource.RLIMIT_FSIZE, (4 << 20, 4 << 20))
                p = subprocess.run(['showfeat'] + cmd, cwd=work, capture_output=True, timeout=5, preexec_fn=limit)
                rc, err = p.returncode, p.stderr.decode('latin-1')
                if rc == -25: hang = True
            except subprocess.TimeoutExpired:
                rc, err, hang = -9, '', True
                subprocess.run(['pkill', '-9', 'showfeat'])
            q = os.path.join(work, 'o.sf')
            outtext = open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
            if err.startswith('Display features of a sequence in pretty format\n'): err = err[len('Display features of a sequence in pretty format\n'):]
            memerr = False
            if not hang and shutil.which('valgrind') and (rc != 0 or len(rows) % 5 == 0):
                vgargs = [a.replace('o.sf', 'ovg.sf') for a in cmd]
                try:
                    vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'showfeat'] + vgargs, cwd=work, capture_output=True, timeout=120)
                    memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                except subprocess.TimeoutExpired:
                    memerr = True
            rows.append({'id': len(rows), 'fmt': fmt, 'how': how, 'text': text, 'spec': spec, 'args': args, 'rc': rc, 'out': outtext, 'err': err, 'hang': hang, 'memerr': memerr})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'died' if r['rc'] else 'out' if r['out'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20265806)
