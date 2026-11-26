#!/usr/bin/env python3
"""B09 WP75 (F346 - F349): converts the returned archive of the Docker package docs/Validation/review-2026-09/handoff/A2-19c-fastapair (returned/A2-19c-fastapair-out.zip: 180 sets of
`mafft --fastapair|--fastswpair ...` run WITHOUT --quiet with a logging wrapper around the real fasta34 (FASTA 34.26.5, sha256 7491909b...), MAFFT 7.505-1; per set the return code, rows, names, tree / hat2 file, standard error, every raw
`fasta34` call of `dndfast7` (arguments, query file, library file, raw output) and a direct run of `dndfast7 [-A]` in an empty directory (its hat2, hat3 and standard error)) into two test oracles in
tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment:
  fasta34_align_oracle.tsv.gz    id (set.call) <TAB> arguments <TAB> query file <TAB> library file <TAB> raw output (escaped: \\ \\n \\t \\r); the process-dependent parts are normalised: /tmp/query-N -> /tmp/query-P (and data), the path of the
                                 program -> fasta34, the clock lines -> "Scan time:  T", " start: T done: T"; the values that the program prints from uninitialised memory ("mp_KS: -0.0000 (N=33) at  86": n_chi_sq and max_i are never assigned when
                                 the library has at most 20 sequences; the same two numbers in all 724 such calls of this build) -> "(N=UNINIT) at UNINIT"
  mafft_fastapair_oracle.tsv.gz  id <TAB> options <TAB> seqs ('|') <TAB> rc <TAB> rows ('|') <TAB> compact JSON: names, tree, hat2 (the .hat2 file of the run), stderr (the process number of "Segmentation fault" -> P),
                                 direct = {rc, stderr, hat2, hat3 (the %p of the next block: (nil) -> NIL, an address -> PTR)}
The archive is untrusted data: it is parsed, never executed. Usage: python3 -I convert_fastapair_returned.py <returned.zip> [<more.zip> ...] [suffix: '2' for the second run, handoff/A2-19c-fastapair2]"""
import gzip, io, json, os, re, sys, tarfile, tempfile, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment')


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def unesc(s):
    return re.sub(r'\\(.)', lambda m: {'n': '\n', 't': '\t', 'r': '\r'}.get(m.group(1), m.group(1)), s, flags=re.S)


def normalise_raw(t):
    t = re.sub(r'/tmp/(query|data|fasta)-\d+', r'/tmp/\1-P', t)
    t = t.replace('/usr/local/bin/fasta34', 'fasta34')
    t = re.sub(r'^ start: .* done: .*$', ' start: T done: T', t, flags=re.M)
    t = re.sub(r'(Scan time: )\s*[0-9.]+', r'\1 T', t)
    t = re.sub(r'(Display time: )\s*[0-9.]+', r'\1 T', t)
    # n_chi_sq and max_i of prhist are never assigned without a histogram: the program prints what is on its stack
    t = re.sub(r'^(; mp_KS: -0\.0000) \(N=\d+\) at\s+\d+$', r'\1 (N=UNINIT) at UNINIT', t, flags=re.M)
    return t


def normalise_err(t):
    return re.sub(r'(/usr/bin/mafft: line \d+:)\s+\d+ (Segmentation fault)', r'\1 P \2', t)


def normalise_hat3(t):
    return re.sub(r' (\(nil\)|0x[0-9a-f]+)$', lambda m: ' NIL' if m.group(1) == '(nil)' else ' PTR', t, flags=re.M)


def read_zip(archive, td, tag):
    with zipfile.ZipFile(archive) as z:
        for n in z.namelist():
            assert not n.startswith('/') and '..' not in n, n
        d = os.path.join(td, tag)
        z.extractall(d)
    return os.path.join(d, 'out')


def main(archives, suffix=''):
    """archives: one or more returned zips; a set takes every part it finds in any of them (a run that was stopped before it wrote results.tsv has the raw calls, the direct runs and the
    standard error of the sets it finished; a later run supplies the complete `mafft` runs). Sets with no results row get rc 'NA' (no complete run to compare)."""
    inputs = os.path.join(HERE, '..', '..', 'docs', 'Validation', 'review-2026-09', 'handoff', 'A2-19c-fastapair' + suffix, 'inputs', 'sets.tsv')
    in_sets = {}
    if suffix != '':
        for line in open(inputs, encoding='utf-8', newline='\n'):
            f = line.rstrip('\n').split('\t')
            in_sets[f[0]] = (f[1], f[2], json.loads(f[3]) if len(f) > 3 else {})
    with tempfile.TemporaryDirectory() as td:
        outs = [read_zip(a, td, 'z%d' % i) for i, a in enumerate(archives)]
        results = {}
        for out in outs:
            path = os.path.join(out, 'results.tsv')
            if os.path.exists(path):
                for line in open(path, encoding='latin-1', newline='\n'):
                    f = line.rstrip('\n').split('\t')
                    assert len(f) in (9, 10), len(f)          # run 4 of the package: a 10th column holds the error of the set (the script catches errors per set)
                    results[f[0]] = (f + [''] * (10 - len(f)), out)
        ids = sorted({n.split('.')[0] for out in outs for n in os.listdir(os.path.join(out, 'sets'))}, key=int)
        calls, sets = [], []
        for k in ids:
            def find(ext):
                for out in outs:
                    p = os.path.join(out, 'sets', k + ext)
                    if os.path.exists(p):
                        return p
                return None
            p_dir = find('.direct.json.gz')
            p_log = find('.fastalog.tar.gz')
            if p_dir is None or p_log is None:
                continue
            direct = json.load(gzip.open(p_dir, 'rt', encoding='utf-8'))
            more = {}
            if k in results:
                f, out = results[k]
                _, opt, seqs, rc, rws, names, sha, tree, hat2, error = f
                # the standard error of the complete run belongs to the run that has its results row
                stderr = gzip.open(os.path.join(out, 'sets', k + '.stderr.gz'), 'rb').read().decode('latin-1')
                full = error == ''
                if not full:
                    rc = 'NA'
                for ext, key in (('.stdout.gz', 'stdout'), ('.files.json.gz', 'outfiles')):      # run 4: the raw standard output (CLUSTAL, PHYLIP, ...) and the files the run left in its directory
                    pth = os.path.join(out, 'sets', k + ext)
                    if os.path.exists(pth):
                        raw = gzip.open(pth, 'rb').read()
                        more[key] = raw.decode('latin-1') if key == 'stdout' else json.loads(raw.decode('utf-8'))
                if error:
                    more['error'] = unesc(error)
            else:
                opt, seqs, _ = in_sets[k]
                rc, rws, names, tree, hat2 = 'NA', '', '', '', ''
                stderr = gzip.open(find('.stderr.gz'), 'rb').read().decode('latin-1')
                full = False
            sets.append('\t'.join([k, opt, seqs, rc, rws, json.dumps({
                'names': names.split('|'), 'tree': unesc(tree), 'hat2': unesc(hat2), 'stderr': normalise_err(stderr), 'full': full, 'files': in_sets[k][2] if k in in_sets else {}, **more,
                'direct': {'rc': direct['rc'], 'stderr': normalise_err(direct['stderr']), 'hat2': direct['hat2'], 'hat3': normalise_hat3(direct['hat3'])}}, separators=(',', ':'))]))
            with tarfile.open(fileobj=io.BytesIO(open(p_log, 'rb').read()), mode='r:gz') as t:
                members = {m.name: t.extractfile(m).read().decode('latin-1') for m in t.getmembers() if m.isfile()}
            for i in sorted({n.split('.')[0] for n in members}):
                args = members[i + '.args'].strip()
                args = re.sub(r'/tmp/(query|data|fasta)-\d+', r'/tmp/\1-P', args)
                calls.append('\t'.join(['%s.%s' % (k, i), args, esc(members[i + '.query']), esc(members[i + '.data']), esc(normalise_raw(members[i + '.out']))]))
        assert suffix != '' or len(sets) == 180, len(sets)
        for name, rows in (('fasta34_align%s_oracle.tsv.gz' % suffix, calls), ('mafft_fastapair%s_oracle.tsv.gz' % suffix, sets)):
            with open(os.path.join(DATA, name), 'wb') as raw_f, gzip.GzipFile(fileobj=raw_f, mode='wb', compresslevel=9, mtime=0, filename='') as f:
                f.write(('\n'.join(rows) + '\n').encode('utf-8'))
        print(len(calls), 'fasta34 calls,', len(sets), 'sets,', sum(1 for r in sets if '"full":true' in r), 'with a complete run')


if __name__ == '__main__':
    args = sys.argv[1:]
    suffix = ''
    if args and not args[-1].endswith('.zip'):
        suffix = args.pop()
    main(args, suffix)
