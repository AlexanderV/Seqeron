#!/usr/bin/env python3
"""A5-8c3h (B09 F262, WP58) oracle: EMBOSS 6.6.0 `featreport` (emboss/featreport.c: the first feature table of a feature file reported for a sequence) and `twofeat` (see generate_twofeat_oracle.py)
    featreport -auto -sequence in.dat -features f.dat -outfile o.rep [-sformat1 F] [-sbegin1 N] [-send1 N] [-sreverse1] [-fbegin2 N] [-fend2 N] [-rformat3 FORMAT] [-raccshow3] [-rdesshow3] [-rscoreshow3|-norscoreshow3]
              [-rstrandshow3|-norstrandshow3] [-rusashow3] [-rmaxall3 N] [-rmaxseq3 N]
Inputs: feature files made as in generate_seqret_oracle.py (embl, genbank, swiss, nbrf, gff3, refseqp, gff2 with feature tables) and EMBL entries with CDS / mRNA / exon features (generate_coderet_oracle.py), as the sequence file itself or beside
FASTA sequences, damaged files; every report format of the library.
JSONL.GZ rows: id, fmt, how, text (the sequence file), ftext (the feature file), args, cmd (the entries of the command line the header echoes), epoch, rc, out, err, hang, memerr
Usage: python3 scripts/alignment/generate_featreport_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0, valgrind, gcc)
"""
import gzip, json, os, random, re, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402
import generate_coderet_oracle as C  # noqa: E402

RFILES = float(os.environ.get('REPORT_FILES_SHARE', '0'))      # the share of the runs with the file name qualifiers -rname3 -rextension3 -rdirectory3 and no -outfile (B09 WP60)
FQUAL = float(os.environ.get('FEATURE_QUAL_SHARE', '0'))       # the share of the runs with -fformat2 -fopenfile2 -iquery2 -ioffset2 -fask2 -freverse2 -fcircular2 (B09 WP60)
FMTS = ['own', 'own', 'own', 'embl', 'genbank', 'gff3', 'gff2', 'swiss', 'nbrf', 'refseqp', 'own', 'embl', 'genbank']
RFORMATS = ['simple', 'embl', 'genbank', 'gff', 'gff3', 'pir', 'swiss', 'dasgff', 'debug', 'listfile', 'dbmotif', 'diffseq', 'draw', 'restrict', 'excel', 'feattable', 'motif', 'nametable', 'regions', 'seqtable', 'srs', 'summary', 'table', 'tagseq']
RFORMATS_W = ['gff', 'gff', 'gff3', 'table', 'table', 'simple', 'excel', 'embl', 'genbank', 'dasgff', 'debug', 'srs', 'listfile', 'regions', 'feattable', 'seqtable', 'nametable', 'tagseq', 'motif', 'summary', 'draw', 'restrict', 'dbmotif', 'diffseq', 'pir', 'swiss']


def echo(name, value):
    if value is None: return name
    if value == '' or any(c in value for c in ' \t*?[]{}|!&^'): return name + ' "' + value + '"'
    return name + ' ' + value


def entries_of(args):
    entries, i = [], 0
    while i < len(args):
        a = args[i]
        if a.startswith('-') and i + 1 < len(args) and (not args[i + 1].startswith('-') or args[i + 1] == '-' or re.match(r'^-\d+(\.\d+)?$', args[i + 1])): entries.append(echo(a, args[i + 1])); i += 2
        else: entries.append(echo(a, None)); i += 1
    return entries


def run(prog, work, args, shim, epoch, timeout=5):
    env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
    try:
        import resource
        def limit():
            resource.setrlimit(resource.RLIMIT_FSIZE, (4 << 20, 4 << 20))
        p = subprocess.run([prog] + args, cwd=work, capture_output=True, timeout=timeout, preexec_fn=limit, env=env)
        return p.returncode, p.stderr.decode('latin-1'), rc_hang(p.returncode)
    except subprocess.TimeoutExpired:
        subprocess.run(['pkill', '-9', prog])
        return -9, '', True


def rc_hang(rc): return rc == -25


def valgrind(prog, work, args):
    if not shutil.which('valgrind'): return False
    vgargs = [a.replace('o.rep', 'ovg.rep') for a in args]
    try:
        vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', prog] + vgargs, cwd=work, capture_output=True, timeout=120)
        return vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
    except subprocess.TimeoutExpired:
        return True


def clean(work):
    for nme in os.listdir(work):
        p = os.path.join(work, nme)
        if os.path.isdir(p):
            for f in os.listdir(p): os.remove(os.path.join(p, f))
            os.rmdir(p)
        else:
            os.remove(p)


def report_file_args(rnd, args, num):
    """the report file qualifiers (-rnameN -rextensionN -rdirectoryN) and a default file name (no -outfile)"""
    if rnd.random() >= RFILES: return args
    args = list(args)
    if rnd.random() < 0.55:
        i = args.index('-outfile'); del args[i:i + 2]
    elif rnd.random() < 0.3:
        i = args.index('-outfile'); args[i + 1] = rnd.choice(['sub/o.rep', 'other/o.rep', 'o.rep', ''])
    if rnd.random() < 0.5: args += ['-rname%d' % num, rnd.choice(['nm', 'N1', 'a b', 'x:y', 'dir/nm', '  lead'])]
    if rnd.random() < 0.4: args += ['-rextension%d' % num, rnd.choice(['abc', 'X1', '.dot', ''])]
    if rnd.random() < 0.5: args += ['-rdirectory%d' % num, rnd.choice(['sub', 'sub', 'sub/', 'nosuch', ''])]
    return args


def report_output(work, hang):
    """(path of the report file the program wrote, its text): the one file of the work directory that is not an input"""
    if hang: return '', ''
    found = []
    for root, _, names in os.walk(work):
        for n in sorted(names):
            rel = os.path.relpath(os.path.join(root, n), work)
            if rel not in ('in.dat', 'f.dat'): found.append(rel)
    if not found: return '', ''
    rel = found[0] if len(found) == 1 else sorted(found)[0]
    return rel, open(os.path.join(work, rel), encoding='latin-1', newline='').read()


def make_text(rnd, gen, fmt):
    if fmt == 'own':
        return ''.join(C.embl_entry(rnd, rnd.choice(['X%d' % rnd.randint(1, 99), 'AB%06d' % rnd.randint(1, 999999), 'Seq_%d' % rnd.randint(1, 9)])) for _ in range(rnd.randint(1, 2))), 'own', 'embl'
    text, how = S.make_input(rnd, gen, fmt)
    return text, how, fmt


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
            ftext, how, ffmt = make_text(rnd, gen, fmt)
            if rnd.random() < 0.08: ftext, how = S.damaged(rnd, ftext), how + '+damaged'
            if rnd.random() < 0.6: text = ftext
            else: text = S.rand_fasta(rnd, 'N', 1)
            args = ['-auto', '-sequence', 'in.dat', '-features', 'f.dat', '-outfile', 'o.rep']
            if rnd.random() < 0.15: args += ['-sformat1', 'embl' if text == ftext and ffmt == 'embl' else 'fasta']
            if rnd.random() < 0.08: args += ['-sbegin1', str(rnd.randint(1, 20))]
            if rnd.random() < 0.08: args += ['-send1', str(rnd.randint(25, 80))]
            if rnd.random() < 0.05: args += ['-sreverse1']
            if rnd.random() < 0.08: args += ['-fbegin2', str(rnd.randint(1, 20))]
            if rnd.random() < 0.08: args += ['-fend2', str(rnd.randint(25, 80))]
            if rnd.random() < 0.85: args += ['-rformat3', rnd.choice(RFORMATS_W)]
            if rnd.random() < 0.15: args.append('-raccshow3')
            if rnd.random() < 0.15: args.append('-rdesshow3')
            if rnd.random() < 0.08: args += ['-rscoreshow3', 'N']
            if rnd.random() < 0.08: args += ['-rstrandshow3', 'N']
            if rnd.random() < 0.12: args.append('-rusashow3')
            if rnd.random() < 0.1: args += ['-rmaxall3', str(rnd.randint(1, 4))]
            if rnd.random() < 0.1: args += ['-rmaxseq3', str(rnd.randint(1, 4))]
            args = report_file_args(rnd, args, 3)
            if rnd.random() < FQUAL:
                if rnd.random() < 0.4: args += ['-fformat2', rnd.choice(['embl', 'genbank', 'gff3', 'junk', 'swiss', 'EMBL'])]
                if rnd.random() < 0.3: args.append('-freverse2')
                if rnd.random() < 0.3: args.append('-fcircular2')
                if rnd.random() < 0.3: args += ['-fopenfile2', rnd.choice(['f.dat', 'nosuch', 'x y'])]
                if rnd.random() < 0.2: args.append('-fask2')
                if rnd.random() < 0.3: args += ['-iquery2', rnd.choice(['X1', 'id:X2', 'seq*', '', 'acc:A'])]
                if rnd.random() < 0.3: args += ['-ioffset2', rnd.choice(['0', '10', '1000', '-5'])]
            clean(work)
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            open(os.path.join(work, 'f.dat'), 'w', newline='\n', encoding='latin-1').write(ftext)
            os.mkdir(os.path.join(work, 'sub'))
            epoch = rnd.randint(978307200, 2051222400)
            rc, err, hang = run('featreport', work, args, shim, epoch)
            outfile, outtext = report_output(work, hang)
            if err.startswith('Read and write a feature table\n'): err = err[len('Read and write a feature table\n'):]
            memerr = (not hang) and (rc != 0 or len(rows) % 5 == 0) and valgrind('featreport', work, args)
            rows.append({'id': len(rows), 'fmt': ffmt, 'how': how, 'text': text, 'ftext': ftext, 'args': args, 'cmd': entries_of(args), 'epoch': epoch, 'rc': rc, 'out': outtext, 'outfile': outfile, 'err': err, 'hang': hang, 'memerr': memerr})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', sorted(Counter((r['fmt'], 'hang' if r['hang'] else 'died' if r['rc'] else 'out' if r['out'] else 'none') for r in rows).items()))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20265808)
