#!/usr/bin/env python3
"""B09 WP61 (F280): the common part of the oracles of the EMBOSS 6.6.0 sequence programs (getorf, transeq, revseq, sixpack, ...): the installed program is run as
    <program> -auto <args> ...      in a work directory that holds the input files
and the standard output, the error stream, the exit status and every file it wrote are kept (JSONL.GZ rows: id, prog, args, files, rc, out, err, outfiles, hang). The C# test
EmbossPrograms_Oracle_Tests runs EmbossProgramRunner.Run(prog, args, files) on every row and compares the four results byte for byte.
A program-specific script defines gen(rnd, i) -> (args, files) and calls main(program, gen, default_n)."""
import gzip, json, os, random, re, resource, subprocess, sys, tempfile, shutil

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_seqret_oracle as S  # noqa: E402

IUPAC_N = 'RYSWKMBDHVN'


def dna(rnd, n, orfs=True, iupac=0.0, case='upper', rna=False, x=False):
    """a random DNA sequence of n bases with planted ORFs (ATG ... stop) and, optionally, ambiguity codes, lower case, U"""
    s = []
    while len(s) < n:
        if orfs and rnd.random() < 0.25:
            body = ''.join(rnd.choice(['GCT', 'GAA', 'AAA', 'CTG', 'GGT', 'AGC', 'CAG', 'TTC', 'CCA', 'GTT', 'ATG', 'TAC']) for _ in range(rnd.randint(2, 40)))
            s.extend('ATG' + body + rnd.choice(['TAA', 'TAG', 'TGA']))
        else:
            s.extend(rnd.choice('ACGT') for _ in range(rnd.randint(1, 30)))
    s = s[:n]
    if iupac:
        s = [rnd.choice(IUPAC_N) if rnd.random() < iupac else c for c in s]
    if x and s and rnd.random() < 0.5:
        s[rnd.randrange(len(s))] = rnd.choice('Xx')
    t = ''.join(s)
    if case == 'lower':
        t = t.lower()
    elif case == 'mixed':
        t = ''.join(c.lower() if rnd.random() < 0.3 else c for c in t)
    if rna:
        t = t.replace('T', 'U').replace('t', 'u')
    return t


def protein(rnd, n, extra=0.05):
    return ''.join(rnd.choice('ACDEFGHIKLMNPQRSTVWY') if rnd.random() > extra else rnd.choice('BZX*') for _ in range(n))


def fasta(rnd, seqs, width=None):
    """seqs: list of (name, desc, residues)"""
    out = ''
    for name, desc, s in seqs:
        w = width or rnd.choice([60, 60, 60, 70, 80, 50, 1000, 10])
        out += '>' + name + (' ' + desc if desc else '') + '\n' + ''.join(s[j:j + w] + '\n' for j in range(0, len(s), w))
    return out


def rand_seqs(rnd, kind='n', nseq=None, lens=(0, 300), **kw):
    """a list of (name, desc, residues)"""
    res = []
    for i in range(nseq or rnd.choice([1, 1, 1, 2, 2, 3, 4])):
        n = rnd.choice([rnd.randint(lens[0], lens[1]), rnd.randint(lens[0], max(lens[0], 60)), rnd.randint(1, 12)])
        s = dna(rnd, n, **kw) if kind == 'n' else protein(rnd, n)
        res.append((S.rand_name(rnd) + ('' if i == 0 else str(i)), S.rand_desc(rnd), s))
    return res


def convert(text, fmt, d):
    """the text in another sequence format, by the installed seqret"""
    if fmt == 'fasta':
        return text
    p = os.path.join(d, 'conv.fa')
    open(p, 'w', newline='\n', encoding='latin-1').write(text)
    r = subprocess.run(['seqret', '-auto', '-sequence', 'conv.fa', '-osformat', fmt, '-outseq', 'stdout'], cwd=d, capture_output=True, timeout=10)
    return r.stdout.decode('latin-1') if r.returncode == 0 and r.stdout else text


def clock_shim(d):
    """compiles scripts/alignment/faketime.c (time() = $FAKE_EPOCH, gettimeofday() = $FAKE_EPOCH / $FAKE_USEC) into d and returns its path"""
    shim = os.path.join(d, 'faketime.so')
    subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c'), '-ldl'], check=True)
    return shim


def run(prog, args, files, work, timeout=8, env=None, timetoday=None):
    """runs the program; returns (rc, stdout, stderr, outfiles, hang). timetoday: the value of the EMBOSS variable of that name (HOME = the work directory with a .embossrc that holds `SET timetoday "VALUE"`): the date of the run and seed 0"""
    for nme in os.listdir(work):
        path = os.path.join(work, nme)
        shutil.rmtree(path) if os.path.isdir(path) else os.remove(path)
    if timetoday is not None:
        rc_text = 'SET timetoday "%s"\n' % timetoday
        open(os.path.join(work, '.embossrc'), 'w').write(rc_text)
        env = dict(env or os.environ, HOME=work)
        files = dict(files)
        files['.embossrc'] = rc_text
    for k, v in files.items():
        if '/' in k:
            os.makedirs(os.path.join(work, os.path.dirname(k)), exist_ok=True)
        open(os.path.join(work, k), 'w', newline='\n', encoding='latin-1').write(v)

    def limit():
        resource.setrlimit(resource.RLIMIT_FSIZE, (8 << 20, 8 << 20))
    hang = False
    proc = subprocess.Popen([prog, '-auto'] + args, cwd=work, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.DEVNULL, preexec_fn=limit, env=env)
    try:
        outb, errb = proc.communicate(timeout=timeout)
        rc = proc.returncode
        if rc in (-25, -11, -6):
            hang = rc == -25
    except subprocess.TimeoutExpired:
        proc.kill()
        outb, errb = proc.communicate()
        rc, hang = -9, True
        outb = errb = b''
    outfiles = {}
    for root, _, names in os.walk(work):
        for nme in names:
            rel = os.path.relpath(os.path.join(root, nme), work)
            if rel not in files:
                outfiles[rel] = open(os.path.join(root, nme), 'rb').read().decode('latin-1')
            elif open(os.path.join(root, nme), 'rb').read().decode('latin-1') != files[rel]:
                outfiles[rel] = open(os.path.join(root, nme), 'rb').read().decode('latin-1')
    return rc, outb.decode('latin-1'), errb.decode('latin-1'), outfiles, hang


def valgrind_bad(prog, args, files, work, env=None, tolerate_underflow=False):
    """True when the run reads or writes invalid memory (valgrind -q reports an error): such a run has no defined result and is refused by the port"""
    if not shutil.which('valgrind'):
        return False
    for nme in os.listdir(work):
        path = os.path.join(work, nme)
        shutil.rmtree(path) if os.path.isdir(path) else os.remove(path)
    for k, v in files.items():
        if '/' in k:
            os.makedirs(os.path.join(work, os.path.dirname(k)), exist_ok=True)
        open(os.path.join(work, k), 'w', newline='\n', encoding='latin-1').write(v)
    try:
        vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', prog, '-auto'] + args, cwd=work, capture_output=True, timeout=120, env=env, stdin=subprocess.DEVNULL)
    except subprocess.TimeoutExpired:
        raise TimeoutError('valgrind did not finish: the row is dropped (a slow machine must not make a refusal)')
    if tolerate_underflow:
        # a read of one byte before a heap block is the high byte of the size field of the chunk header (always 0): the loop `for(...; p[pos]; ...)` of the programs stops there as at a NUL (benign, deterministic)
        records = re.split(rb'\n==\d+== \n', vg.stderr)
        records = [r for r in records if not (re.search(rb'Invalid read of size 1\b', r) and re.search(rb'is 1 bytes before a block of size', r))]
        return vg.returncode == 99 and any(re.search(rb'(Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', r) for r in records) or any(re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', r, re.M) for r in records)
    return vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None


def main(prog, gen, default_n=300, default_seed=20261701):
    """gen(rnd, i, dir) -> (args, files) or (args, files, extra): extra = {'usec': N} (the microseconds gettimeofday returns when the random number generator is seeded) or {'timetoday': 'yyyy-mm-dd hh:mm:ss'} (the value of the variable).
    Every row also holds 'epoch': the date of the run (time() of the run; the writers of genbank / embl / ... print it), drawn from a generator of its own so that the rows do not depend on it."""
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else default_n
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else default_seed
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        work = os.path.join(d, 'w')
        os.mkdir(work)
        gdir = os.path.join(d, 'g')
        os.mkdir(gdir)
        shim = clock_shim(d)
        while len(rows) < n:
            got = gen(rnd, len(rows), gdir)
            args, files = got[0], got[1]
            extra = got[2] if len(got) > 2 else {}
            epoch = random.Random(seed * 1000003 + len(rows)).randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch), TZ='UTC')
            if 'usec' in extra:
                env['FAKE_USEC'] = str(extra['usec'])
            rc, o, e, outfiles, hang = run(prog, args, files, work, env=env, timetoday=extra.get('timetoday'))
            rc2, o2, e2, outfiles2, hang2 = run(prog, args, files, work, env=env, timetoday=extra.get('timetoday'))
            if hang and rc == -9:
                print('row', len(rows), 'timed out: dropped', args, file=sys.stderr)       # slow rather than defined (e.g. the 28 s of an NCBI header with two barred tokens): no row
                continue
            if (rc, o, e, outfiles) != (rc2, o2, e2, outfiles2):
                print('row', len(rows), 'not deterministic: dropped', args, file=sys.stderr)
                continue
            try:
                vgbad = extra.get('vg') and rc >= 0 and not hang and valgrind_bad(prog, args, dict(files, **({'.embossrc': 'SET timetoday "%s"\n' % extra['timetoday']} if extra.get('timetoday') is not None else {})), work, dict(env, **({'HOME': work} if extra.get('timetoday') is not None else {})), tolerate_underflow=extra['vg'] == 'underflow')
            except TimeoutError as ex:
                print('row', len(rows), ex, args, file=sys.stderr)
                continue
            if vgbad:
                rc, hang, o, e, outfiles = -11, False, '', '', {}                  # invalid memory access in the reference: refused by the port
            if rc >= 0 and not hang and 'Memory allocation failed' in e:
                rc, hang, o, e, outfiles = -11, False, '', '', {}              # the allocation of a gigantic array (an unsigned index below 0): the result depends on the memory of the machine: refused by the port
            row = {'id': len(rows), 'prog': prog, 'args': args, 'files': {k: v for k, v in files.items()}, 'rc': rc, 'out': o, 'err': e, 'outfiles': outfiles, 'hang': hang, 'epoch': epoch}
            row.update({k: v for k, v in extra.items()})
            rows.append(row)
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows:
            fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', 'rc', sorted(Counter(r['rc'] for r in rows).items()), 'nonempty out', sum(1 for r in rows if r['out'] or r['outfiles']), 'err', sum(1 for r in rows if r['err']))
