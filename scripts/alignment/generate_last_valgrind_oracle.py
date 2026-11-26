#!/usr/bin/env python3
"""B09 WP85 (A2-19c-last, F376): the valgrind proof of the heap overflow of `pairlocalalign -R` for EVERY set of the returned archive of the Docker package handoff/A2-19c-last (run 2).
The installed Debian mafft 7.505-1 runs each set with LAST replaced by a shim that answers every call of lastdb / lastal with the recorded output of the package (so the pair stage sees exactly the alignments of LAST 1256),
`pairlocalalign` runs under valgrind (memcheck). Output: tests/Seqeron/Seqeron.Genomics.Last.Tests/TestData/Alignment/last_valgrind_oracle.tsv: set <TAB> invalid writes counted by valgrind <TAB> 1 if the oracle's run crashed.
Needs: /usr/bin/mafft 7.505, /usr/bin/valgrind, the archive unzipped in an EMPTY directory (argument 1: its `out` directory). Nothing from the archive is executed; the shim is written here.
Usage: python3 -I scripts/alignment/generate_last_valgrind_oracle.py <out directory of the unzipped archive>"""
import glob, gzip, json, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Last.Tests', 'TestData', 'Alignment')
SHIM = r'''#!/usr/bin/python3
import sys, os, json, gzip
rec = json.load(gzip.open(os.environ['REPLAY_SET'], 'rt'))
args = sys.argv[1:]
q = open(args[-1]).read() if os.path.exists(args[-1]) else None
want = '\n'.join(args) + '\n'
for c in rec['run']['calls']:
    it = c['items']
    if c['prog'] == 'lastal' and it.get('args') == want and it.get('in.' + args[-1]) == q:
        sys.stdout.write(it.get('stdout', '')); sys.stderr.write(it.get('stderr', '')); sys.exit(int(it.get('rc', '0').strip() or 0))
sys.exit(99)
'''


def main(out):
    work = tempfile.mkdtemp(prefix='lastvg')
    bindir = os.path.join(work, 'bin')
    os.makedirs(bindir)
    for f in glob.glob('/usr/lib/mafft/lib/mafft/*'):
        os.symlink(f, os.path.join(bindir, os.path.basename(f)))
    os.unlink(os.path.join(bindir, 'pairlocalalign'))
    open(os.path.join(bindir, 'pairlocalalign'), 'w').write('#!/bin/sh\nexec valgrind --log-file="$VG_LOG.%p" -q /usr/lib/mafft/lib/mafft/pairlocalalign "$@"\n')
    open(os.path.join(bindir, 'lastal'), 'w').write(SHIM)
    open(os.path.join(bindir, 'lastdb'), 'w').write('#!/bin/sh\nexit 0\n')
    for n in ('pairlocalalign', 'lastal', 'lastdb'):
        os.chmod(os.path.join(bindir, n), 0o755)
    rows = []
    for line in open(os.path.join(out, 'index.tsv'), encoding='utf-8'):
        k = line.split('\t')[0]
        setfile = os.path.join(out, 'sets', k + '.json.gz')
        rec = json.load(gzip.open(setfile, 'rt', encoding='utf-8'))
        d = os.path.join(work, 'run_' + k)
        os.makedirs(d)
        open(os.path.join(d, 'in.fa'), 'w').write(rec['in.fa'])
        if rec['new.fa']:
            open(os.path.join(d, 'new.fa'), 'w').write(rec['new.fa'])
        log = os.path.join(work, 'vg_' + k)
        env = dict(os.environ, MAFFT_BINARIES=bindir, REPLAY_SET=setfile, VG_LOG=log)
        subprocess.run(rec['argv'], cwd=d, env=env, capture_output=True, timeout=600)
        invalid = sum(len(re.findall('Invalid write', open(f, errors='replace').read())) for f in glob.glob(log + '.*'))
        err = rec['run']['stderr']
        rows.append('%s\t%d\t%d' % (k, invalid, 1 if ('Segmentation fault' in err or 'Aborted' in err) else 0))
    with open(os.path.join(DATA, 'last_valgrind_oracle.tsv'), 'w', newline='\n') as f:
        f.write('\n'.join(rows) + '\n')
    print(len(rows), 'sets', sum(1 for r in rows if int(r.split('\t')[1]) > 0), 'with invalid writes')


if __name__ == '__main__':
    main(sys.argv[1])
