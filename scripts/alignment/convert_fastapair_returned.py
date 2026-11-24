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
The archive is untrusted data: it is parsed, never executed. Usage: python3 -I convert_fastapair_returned.py <returned.zip> [suffix: '2' for the second run, handoff/A2-19c-fastapair2]"""
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


def main(archive, suffix=''):
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(archive) as z:
            for n in z.namelist():
                assert not n.startswith('/') and '..' not in n, n
            z.extractall(td)
        out = os.path.join(td, 'out')
        calls, sets = [], []
        for line in open(os.path.join(out, 'results.tsv'), encoding='latin-1', newline='\n'):
            f = line.rstrip('\n').split('\t')
            assert len(f) == 9, len(f)
            k, opt, seqs, rc, rws, names, sha, tree, hat2 = f
            stderr = gzip.open(os.path.join(out, 'sets', k + '.stderr.gz'), 'rb').read().decode('latin-1')
            direct = json.load(gzip.open(os.path.join(out, 'sets', k + '.direct.json.gz'), 'rt', encoding='utf-8'))
            sets.append('\t'.join([k, opt, seqs, rc, rws, json.dumps({
                'names': names.split('|'), 'tree': unesc(tree), 'hat2': unesc(hat2), 'stderr': normalise_err(stderr),
                'direct': {'rc': direct['rc'], 'stderr': normalise_err(direct['stderr']), 'hat2': direct['hat2'], 'hat3': normalise_hat3(direct['hat3'])}}, separators=(',', ':'))]))
            with tarfile.open(fileobj=io.BytesIO(open(os.path.join(out, 'sets', k + '.fastalog.tar.gz'), 'rb').read()), mode='r:gz') as t:
                members = {m.name: t.extractfile(m).read().decode('latin-1') for m in t.getmembers() if m.isfile()}
            ids = sorted({n.split('.')[0] for n in members})
            for i in ids:
                args = members[i + '.args'].strip()
                args = re.sub(r'/tmp/(query|data|fasta)-\d+', r'/tmp/\1-P', args)
                calls.append('\t'.join(['%s.%s' % (k, i), args, esc(members[i + '.query']), esc(members[i + '.data']), esc(normalise_raw(members[i + '.out']))]))
        assert suffix != '' or len(sets) == 180, len(sets)
        for name, rows in (('fasta34_align%s_oracle.tsv.gz' % suffix, calls), ('mafft_fastapair%s_oracle.tsv.gz' % suffix, sets)):
            with open(os.path.join(DATA, name), 'wb') as raw_f, gzip.GzipFile(fileobj=raw_f, mode='wb', compresslevel=9, mtime=0, filename='') as f:
                f.write(('\n'.join(rows) + '\n').encode('utf-8'))
        print(len(calls), 'fasta34 calls,', len(sets), 'sets')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else '')
