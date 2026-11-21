#!/usr/bin/env python3
"""A3-9e2 oracle: `combineMUMs` and the script `run-mummer3` of MUMmer 3.23 (Debian 3.23+dfsg-8).
Rows (JSON lines):
  {"kind":"cmums","ref","qry","gaps","args":[...],"files":{name:text},"out","err","rc"}
        combineMUMs <args> ref.fa qry.fa gaps.txt run in a scratch directory; args are the options (-D -e x -n -N k -q tag -r tag -S -t -v k -x -W name); "files" = every file the program wrote (the -W file, the .cover files);
        "out" / "err" = standard output / error (latin-1, file names untouched); "rc" = the exit status (134 = abort of an assert, 139 = segmentation fault, 255 = exit (EXIT_FAILURE)).
  {"kind":"runmummer3","ref","qry","prefix","files":{ext:text},"out","err","rc"}
        run-mummer3 ref.fa qry.fa <prefix>: the files .out .gaps .align .errorsgaps, the script's standard output and error (the mummer timing lines removed)
  {"kind":"ub", ...}  a crafted row whose result is undefined in the program (valgrind: invalid read; or the output changes under MALLOC_PERTURB_) - recorded, not compared.
The gaps files come from the real pipeline (mummer -mumreference -b | mgaps with random -l / -f / -s), from raw mummer matches (overlapping matches: negative gap lengths) and from hand-written lists.
Usage: python3 scripts/alignment/generate_combinemums_oracle.py [seed=20261204] | gzip > combinemums_oracle.jsonl.gz"""
import json, os, random, shutil, subprocess, sys

WORK = '/tmp/b09cm_%d' % os.getpid()
COMP = str.maketrans('ACGTacgt', 'TGCAtgca')


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261204
    n_cm = int(sys.argv[2]) if len(sys.argv) > 2 else 330
    n_rm3 = int(sys.argv[3]) if len(sys.argv) > 3 else 230
    rnd = random.Random(seed)

    def rseq(k, alpha='acgt'): return ''.join(rnd.choice(alpha) for _ in range(k))
    def rc(s): return s.translate(COMP)[::-1]
    def mutate(s, rate):
        o = []
        for c in s:
            r = rnd.random()
            if r < rate * .6: o.append(rnd.choice('acgt'))
            elif r < rate * .75: pass
            elif r < rate: o.append(c + rnd.choice('acgt'))
            else: o.append(c)
        return ''.join(o)
    def fasta(recs, width=60, upper=False):
        out = []
        for nm, sq in recs:
            out.append('>' + nm + '\n')
            if upper: sq = sq.upper()
            for i in range(0, len(sq), width): out.append(sq[i:i + width] + '\n')
        return ''.join(out)
    def run(cmd, **kw): return subprocess.run(cmd, cwd=WORK, capture_output=True, **kw)
    def reset():
        shutil.rmtree(WORK, ignore_errors=True); os.makedirs(WORK)
    def lat(b): return b.decode('latin-1')
    def noise(s, frac):
        """N and other IUPAC characters inside the sequence"""
        o = list(s)
        for _ in range(int(len(s) * frac) + (1 if rnd.random() < .3 else 0)):
            if not o: break
            p = rnd.randrange(len(o)); o[p] = rnd.choice('nnnnrykm')
            if rnd.random() < .3:
                for q in range(p, min(len(o), p + rnd.randint(2, 12))): o[q] = 'n'
        return ''.join(o)

    def make_pair(big=False):
        """a reference and a (multi-record) query that share diverged segments, some of them reverse complemented"""
        L = rnd.randint(900, 2800) if not big else rnd.randint(10600, 11600)
        base = rseq(L)
        if rnd.random() < .25:           # a repeat inside the reference
            r0 = rnd.randint(0, L // 2); seg = base[r0:r0 + rnd.randint(40, 300)]
            base = base[:L // 2] + mutate(seg, .02) + base[L // 2:]
        ref = base
        recs = []
        for i in range(rnd.choice([1, 1, 2, 2, 3])):
            a = rnd.randint(0, len(base) // 3); b = rnd.randint(len(base) // 2, len(base))
            piece = mutate(base[a:b], rnd.choice([0, .01, .02, .04, .07, .12, .2]))
            if rnd.random() < .4:        # a second, independent divergent block
                piece = piece[:len(piece) // 2] + rseq(rnd.randint(1, 120)) + piece[len(piece) // 2:]
            if rnd.random() < .45: piece = rc(piece)
            if rnd.random() < .5: piece = rseq(rnd.randint(0, 150)) + piece + rseq(rnd.randint(0, 150))
            recs.append(('q%d desc %d' % (i, i), piece))
        return ref, recs

    rows = []

    def runmummer3_row(ref_text, qry_text):
        reset()
        open(WORK + '/ref.fa', 'w').write(ref_text); open(WORK + '/qry.fa', 'w').write(qry_text)
        p = run(['run-mummer3', 'ref.fa', 'qry.fa', 'pf'])
        files = {}
        for ext in ('out', 'gaps', 'align', 'errorsgaps', 'cover'):
            fp = '%s/pf.%s' % (WORK, ext)
            files[ext] = open(fp, encoding='latin-1', newline='').read() if os.path.exists(fp) else None
        extra = sorted(f for f in os.listdir(WORK) if f not in ('ref.fa', 'qry.fa') and not f.startswith('pf.'))
        err = '\n'.join(l for l in lat(p.stderr).split('\n') if not l.startswith(('# CONSTRUCTIONTIME', '# COMPLETETIME', '# SPACE')))
        row = {'kind': 'runmummer3', 'ref': ref_text, 'qry': qry_text, 'prefix': 'pf', 'files': files, 'out': lat(p.stdout), 'err': err, 'rc': p.returncode, 'extra': extra}
        if p.returncode == 0:
            q = run(['run-mummer3', 'ref.fa', 'qry.fa', 'pg'], env=dict(os.environ, MALLOC_PERTURB_='165'))
            row['ub'] = any((open('%s/pg.%s' % (WORK, e), encoding='latin-1', newline='').read() if os.path.exists('%s/pg.%s' % (WORK, e)) else None) != files[e] for e in ('out', 'gaps', 'align', 'errorsgaps'))
        return row

    # ---------- run-mummer3
    for k in range(n_rm3):
        ref, recs = make_pair(big=(k % 100 == 99))
        if k % 7 == 3: ref = noise(ref, .004)
        if k % 9 == 4: recs = [(nm, noise(sq, .003)) for nm, sq in recs]
        ref_text = '>ref sequence\n' + '\n'.join(ref[i:i + 60] for i in range(0, len(ref), 60)) + '\n'
        qry_text = fasta(recs, rnd.choice([50, 60, 70, 80]), upper=(k % 11 == 5))
        rows.append(runmummer3_row(ref_text, qry_text))
    # the failing runs of the script
    for ref_text, qry_text in [('>r\nacgtacgtacgtacgtacgtacgtacgtacgt\n', '>q\nttttttttttttttttttttttt\n'),          # no match: mgaps prints the labels only
                               ('>r\n\n', '>q\nacgtacgtacgtacgtacgtacgtacgt\n'),                                        # empty reference
                               ('', '>q\nacgtacgtacgtacgtacgtacgtacgt\n'),
                               ('>r\nacgtacgtacgtacgtacgtacgtacgt\n', ''),
                               ('>r\nacgtacgtacgtacgtacgtacgtacgt\n', '>q\n')]:
        rows.append(runmummer3_row(ref_text, qry_text))

    # ---------- combineMUMs
    OPTS = [['-D'], ['-S'], ['-n'], ['-x'], ['-t'], ['-e', '.10'], ['-e', '.05'], ['-e', '0.2'], ['-e', '0'], ['-e', '1'], ['-N', '3'], ['-N', '5'], ['-N', '1'], ['-q', 'QQ'], ['-r', 'RR'],
            ['-v', '1'], ['-W', 'my.errors'], ['-W', 'sub/e.out']]
    def pipeline(ref_text, qry_text, l, mg_l, mg_f, mg_s, mode):
        open(WORK + '/ref.fa', 'w').write(ref_text); open(WORK + '/qry.fa', 'w').write(qry_text)
        mm = run(['mummer', mode, '-b', '-l', str(l), 'ref.fa', 'qry.fa'])
        if mm.returncode != 0: return None
        raw = mm.stdout
        mg = subprocess.run(['mgaps', '-l', str(mg_l), '-f', str(mg_f), '-s', str(mg_s)], cwd=WORK, input=raw, capture_output=True)
        return lat(raw), lat(mg.stdout)

    def cm_row(ref_text, qry_text, gaps, flat):
        reset()
        os.makedirs(WORK + '/sub', exist_ok=True)
        open(WORK + '/ref.fa', 'w').write(ref_text); open(WORK + '/qry.fa', 'w').write(qry_text); open(WORK + '/gaps.txt', 'w').write(gaps)
        p = run(['combineMUMs'] + flat + ['ref.fa', 'qry.fa', 'gaps.txt'])
        files = {}
        for dp, _, fs in os.walk(WORK):
            for f in fs:
                rel = os.path.relpath(os.path.join(dp, f), WORK)
                if rel not in ('ref.fa', 'qry.fa', 'gaps.txt'): files[rel] = open(os.path.join(dp, f), encoding='latin-1', newline='').read()
        row = {'kind': 'cmums', 'ref': ref_text, 'qry': qry_text, 'gaps': gaps, 'args': flat, 'files': files, 'out': lat(p.stdout), 'err': lat(p.stderr), 'rc': p.returncode}
        if p.returncode == 0:
            # the genuine non-deterministic class: a run whose output changes when the allocator is tuned differently has read memory it never wrote (the reference buffer is realloc'd to Ref_Len + 1 bytes; the
            # coverage report reads its terminator one byte past the block: when the length is 7 modulo 16 that byte is the size field of the next heap chunk)
            ub = False
            for env in (dict(os.environ, MALLOC_PERTURB_='165'), dict(os.environ, GLIBC_TUNABLES='glibc.malloc.tcache_count=0'), dict(os.environ, MALLOC_MMAP_THRESHOLD_='50000')):
                reset_files = [f for f in files]
                q = run(['combineMUMs'] + flat + ['ref.fa', 'qry.fa', 'gaps.txt'], env=env)
                same = q.stdout == p.stdout and q.stderr == p.stderr
                for f, t in files.items():
                    fp = os.path.join(WORK, f)
                    same = same and os.path.exists(fp) and open(fp, encoding='latin-1', newline='').read() == t
                ub = ub or not same
            ref_len = sum(len(l.strip()) for l in ref_text.split('\n')[1:])
            if ref_len % 16 == 7 and '-x' not in flat:
                row['ref7'] = True
                ub = True
                vg = run(['valgrind', '-q', '--error-exitcode=9'] + ['combineMUMs'] + flat + ['ref.fa', 'qry.fa', 'gaps.txt'])
                row['valgrind'] = vg.returncode == 9 and 'Invalid read of size 1' in lat(vg.stderr) and 'after a block of size %d alloc' % (ref_len + 1) in lat(vg.stderr).replace(',', '')
            row['ub'] = ub
        return row

    for k in range(n_cm):
        ref, recs = make_pair(big=(k % 125 == 124))
        if k % 31 == 7: ref = ref[:len(ref) - (len(ref) - 7) % 16]       # a length of 7 modulo 16
        if k % 5 == 1: ref = noise(ref, .006)
        if k % 6 == 2: recs = [(nm, noise(sq, .006)) for nm, sq in recs]
        ref_text = '>ref sequence\n' + '\n'.join(ref[i:i + 60] for i in range(0, len(ref), 60)) + '\n'
        qry_text = fasta(recs, rnd.choice([50, 60, 70]), upper=(k % 13 == 6))
        reset()
        l = rnd.choice([10, 12, 15, 20, 20, 25])
        got = pipeline(ref_text, qry_text, l, rnd.choice([20, 40, 65, 100]), rnd.choice([.05, .12, .2]), rnd.choice([90, 200, 600, 1000]), rnd.choice(['-mumreference', '-mumreference', '-mum', '-maxmatch']))
        if got is None: continue
        raw, gaps = got
        variant = rnd.random()
        if variant < .12:        # raw mummer output: overlapping and unordered matches, the header lines of mummer
            gaps = raw
        elif variant < .2:       # the first lines only
            gaps = ''.join(gaps.splitlines(True)[:rnd.randint(1, 14)])
        elif variant < .26:      # comment lines
            ls = gaps.splitlines(True)
            ls.insert(rnd.randint(0, len(ls)), '# comment\n')
            gaps = ''.join(ls)
        if not any(not l.startswith(('>', '#')) for l in gaps.splitlines()):
            continue
        args = []
        for _ in range(rnd.choice([0, 1, 1, 2, 3])):
            o = rnd.choice(OPTS)
            if o not in args and not (o[0] == '-W' and any(a[0] == '-W' for a in args)): args.append(o)
        flat = [x for o in args for x in o]
        if '-v' in flat and len(gaps) > 6000: flat = [x for x in flat if x != '-v']; flat = [x for i, x in enumerate(flat)]; args = [o for o in args if o[0] != '-v']
        flat = [x for o in args for x in o]
        rows.append(cm_row(ref_text, qry_text, gaps, flat))

    # ---------- crafted combineMUMs inputs: the failures of the program and its special cases
    R1 = '>r1 first\n' + 'acgtgctagctagcgatcgatcgtagctagctagcatcgatcgactagctagcatgcatcgatcgatgcatgcatgcatcgatgc' * 3 + '\n>r2 second\nacgtacgt\n'
    Q1 = '>q1 one\n' + 'gctagcgatcgatcgtagctagctagcatcgatcgactagctagcatgcatcgatcgatgcatgcatgca' * 3 + '\n>q2 two\n' + 'tcgatcgatgcatgcatgcatcgatgcacgtgctagctagc' * 2 + '\n'
    G1 = ('> q1\n       1        1     40    none      -      -\n      60       60     30    none     19     19\n> q1 Reverse\n'
          '       5        7     35    none      -      -\n')
    crafted = [
        ('', Q1, G1, []),                                              # no reference
        (R1, '', G1, []),                                              # no query record
        (R1, Q1, '', []),                                              # an empty gaps file
        (R1, Q1, '# only a comment\n', []),
        (R1, Q1, '> q1\n> q1 Reverse\n', []),                          # headers only
        (R1, Q1, '> q1\n> q1 Reverse\n> q2\n> q2 Reverse\n> q3\n', []),          # more forward blocks than query records
        (R1, Q1, '       1        1     40    none      -      -\n', []),   # a match before the first header
        (R1, Q1, '> q1\nnot numbers\n', []),
        (R1, Q1, '> q1\n\n', []),
        (R1, Q1, '> q1\n       1        1\n', []),
        (R1, Q1, G1, []), (R1, Q1, G1, ['-t']), (R1, Q1, G1, ['-t', '-D']), (R1, Q1, G1, ['-S']), (R1, Q1, G1, ['-v', '1']), (R1, Q1, G1, ['-v', '2', '-S']),
        (R1, Q1, G1.replace('\n', '\r\n'), []),
        (R1, Q1, '>\n' + G1, ['-t']),                                   # a header line without a tag
        (R1, Q1, G1 + '> q2\n  1 1 30 none - -\n', []),
        (R1, Q1, '> ' + 'x' * 1100 + '\n' + G1, []),                    # a header longer than the line buffer
        (R1, Q1, '> q1\n' + ' ' * 1100 + '1 1 40\n', []),
        ('>r1\nacgtNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNacgtgatcgatcgatcgatgctagctagctagcatgcatgc\n',
         '>q1\nacgtNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNacgtgatcgatcgatcgatgctagctagctagcatgcatgc\n',
         '> q1\n       1        1      4    none      -      -\n      70       70     30    none     65     65\n', ['-N', '5']),
        ('>r1\nacgtNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNacgtgatcgatcgatcgatgctagctagctagcatgcatgc\n',
         '>q1\nacgtNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNacgtgatcgatcgatcgatgctagctagctagcatgcatgc\n',
         '> q1\n       1        1      4    none      -      -\n      70       70     30    none     65     65\n', ['-n']),
        ('>r1\nacgt*1acgt\n', '>q1\nacgt-9acgt\n', '> q1\n  1 1 4 none - -\n', []),         # characters that are not letters
        ('>r1\nAcGtAcGtAcGtAcGtAcGtAcGt\n', '>q1\nacgtacgtacgtacgtacgtacgt\n', '> q1\n  1 1 12 none - -\n  13 13 12 none 1 1\n', []),
        (R1, Q1, '> q1\n       1        1     40    none      -      -\n    9000     9000     30    none     19     19\n', []),      # coordinates beyond the sequences
        (R1, Q1, '> q1\n     200        1     40    none      -      -\n', []),
        (R1, Q1, '> q1\n       1      200     40    none      -      -\n', []),
        (R1, Q1, '> q1\n      50       50     40    none      -      -\n      60       60     40    none     -30    -30\n', []),   # overlapping matches: negative gap lengths
        (R1, Q1, '> q1\n      50       50     40    none      -      -\n      30       30     40    none     -30    -30\n', []),   # a match before the previous one
    ]
    for ref_text, qry_text, gaps, flat in crafted:
        reset()
        os.makedirs(WORK + '/sub', exist_ok=True)
        open(WORK + '/ref.fa', 'w').write(ref_text); open(WORK + '/qry.fa', 'w').write(qry_text); open(WORK + '/gaps.txt', 'w').write(gaps)
        rows.append(cm_row(ref_text, qry_text, gaps, flat))

    # long extensions: a single match in the middle of two 24 kb sequences that differ by 1 - 2 percent is extended in steps of 10 000 bases (MAX_EXTENSION), both ways, on both strands
    for k, (rate, rev, flat) in enumerate([(.01, False, []), (.02, True, []), (.01, False, ['-v', '1']), (.015, True, ['-v', '1', '-S']), (.02, False, ['-D'])]):
        reset()
        big = rseq(24000)
        mid = rnd.randint(11000, 13000)
        q = mutate(big, rate)
        if rev: q = rc(q)
        ref_text = '>big\n' + '\n'.join(big[i:i + 70] for i in range(0, len(big), 70)) + '\n'
        qry_text = '>bq\n' + '\n'.join(q[i:i + 70] for i in range(0, len(q), 70)) + '\n'
        got = pipeline(ref_text, qry_text, 20, 20, .12, 600, '-mumreference')
        if got is None: continue
        ls = got[1].splitlines(True)
        hdr = [i for i, l in enumerate(ls) if l.startswith('>')]
        if len(hdr) < 2: continue
        fwd = [l for l in ls[hdr[0] + 1:hdr[1]] if not l.startswith('#')]
        rv = [l for l in ls[hdr[1] + 1:] if not l.startswith(('>', '#'))]
        pick = (rv if rev else fwd)
        if not pick: continue
        m = pick[len(pick) // 2]
        gaps = ls[hdr[0]] + (m if not rev else '') + ls[hdr[1]] + (m if rev else '')
        rows.append(cm_row(ref_text, qry_text, gaps, flat))

    shutil.rmtree(WORK, ignore_errors=True)
    for r in rows: print(json.dumps(r))
    kinds = {}
    for r in rows:
        key = r['kind'] + (':rc%d' % r['rc'] if r['rc'] else '') + (':ub' if r.get('ub') else '')
        kinds[key] = kinds.get(key, 0) + 1
    sys.stderr.write('rows=%d %s\n' % (len(rows), kinds))


main()
