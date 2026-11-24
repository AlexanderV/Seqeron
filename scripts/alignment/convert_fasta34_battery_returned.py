#!/usr/bin/env python3
"""B09 WP71 (F323 / F324): converts the returned archive of the Docker package docs/Validation/review-2026-09/handoff/A2-19c-fasta34-battery (returned/A2-19c-fasta34-battery-out.zip: the raw `fasta34 -z3 -m10` output of 340 calls
(20 probes of getfastascoremtx + 320 random queries / libraries as splittbfast.c makes them) and 200 sets of `mafft --fastaparttree` run WITHOUT --quiet; MAFFT 7.505-1, fasta34 34.26.5 sha256 7491909b...) into two test oracles in
tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment:
  fasta34_battery_oracle.tsv.gz            id <TAB> kind (p | n) <TAB> query <TAB> library ('|' separated) <TAB> rc <TAB> raw output (escaped: \\ \\n \\t \\r; the temporary directory name -> /tmp/T, the clock lines -> "Scan time:  T", " start: T done: T")
  mafft_fastaparttree_full_oracle.tsv.gz   id <TAB> options <TAB> seqs ('|') <TAB> rc <TAB> rows ('|') <TAB> extra (compact JSON: names, tree (the .tree file, '' when there is none), stderr WITHOUT the process-number line `rm -f /tmp/data-<pid>*  /tmp/query-<pid> /tmp/fasta-<pid>`)
The archive is untrusted data: it is parsed, never executed. Usage: python3 -I convert_fasta34_battery_returned.py <returned.zip>"""
import gzip, json, os, re, sys, tempfile, zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment')
INPUTS = os.path.join(HERE, '..', '..', 'docs', 'Validation', 'review-2026-09', 'handoff', 'A2-19c-fasta34-battery', 'inputs')


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def normalise_raw(t):
    t = re.sub(r'/tmp/tmp[A-Za-z0-9_]+/', '/tmp/T/', t)
    t = re.sub(r'^ start: .* done: .*$', ' start: T done: T', t, flags=re.M)
    t = re.sub(r'(Scan time: )\s*[0-9.]+', r'\1 T', t)
    t = re.sub(r'(Display time: )\s*[0-9.]+', r'\1 T', t)
    return t


def main(archive):
    with tempfile.TemporaryDirectory() as td:
        with zipfile.ZipFile(archive) as z:
            for n in z.namelist():
                assert not n.startswith('/') and '..' not in n, n
            z.extractall(td)
        out = os.path.join(td, 'out')
        # inputs of the battery: the probes are written by the container script, the 340 calls are inputs/battery.tsv
        battery = {}
        lib = ''.join(c * 6 + 'X' * 6 for c in 'ACDEFGHIKLMNPQRSTVWY')
        for c in 'ACDEFGHIKLMNPQRSTVWY':
            battery['probe-' + c] = ('p', c * 6, lib)
        for line in open(os.path.join(INPUTS, 'battery.tsv')):
            k, kind, q, l = line.rstrip('\n').split('\t')
            battery['b' + k] = (kind, q, l)
        rows = []
        for line in open(os.path.join(out, 'battery.tsv'), newline='\n'):
            i, kind, rc = line.rstrip('\n').split('\t')[:3]
            assert battery[i][0] == kind, i
            raw = gzip.open(os.path.join(out, 'battery', i + '.out.gz'), 'rb').read().decode('latin-1')      # bytes: the programs print '\r' (no newline translation)
            rows.append('\t'.join([i, kind, battery[i][1], battery[i][2], rc, esc(normalise_raw(raw))]))
        assert len(rows) == 340, len(rows)   # 20 probes + 320 calls of inputs/battery.tsv
        with open(os.path.join(DATA, 'fasta34_battery_oracle.tsv.gz'), 'wb') as raw_f, gzip.GzipFile(fileobj=raw_f, mode='wb', compresslevel=9, mtime=0, filename='') as f:
            f.write(('\n'.join(rows) + '\n').encode('utf-8'))
        sets = []
        for line in open(os.path.join(out, 'results.tsv'), encoding='latin-1', newline='\n'):
            f = line.rstrip('\n').split('\t')
            assert len(f) == 8, len(f)
            k, opt, seqs, rc, rws, names, sha, tree = f
            stderr = gzip.open(os.path.join(out, 'sets', k + '.stderr.gz'), 'rb').read().decode('latin-1')
            stderr, n = re.subn(r'^rm -f /tmp/data-(\d+)\*  /tmp/query-\1 /tmp/fasta-\1\n', '', stderr, flags=re.M)
            assert n == (0 if rc == '1' else 1), (k, n, rc)       # the one set that stops (rc 1: "Seq 31, too short.") never reaches the cleaning line
            tree = tree.replace('\\n', '\n').replace('\\t', '\t').replace('\\\\', '\\')
            sets.append('\t'.join([k, opt, seqs, rc, rws, json.dumps({'names': names.split('|'), 'tree': tree, 'stderr': stderr}, separators=(',', ':'))]))
        assert len(sets) == 200
        with open(os.path.join(DATA, 'mafft_fastaparttree_full_oracle.tsv.gz'), 'wb') as raw_f, gzip.GzipFile(fileobj=raw_f, mode='wb', compresslevel=9, mtime=0, filename='') as f:
            f.write(('\n'.join(sets) + '\n').encode('utf-8'))
        print(len(rows), 'battery rows,', len(sets), 'sets')


if __name__ == '__main__':
    main(sys.argv[1])
