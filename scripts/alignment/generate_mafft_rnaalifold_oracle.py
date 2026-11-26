#!/usr/bin/env python3
"""B09 WP83 (item A2-19c-rnaalifold): the oracle of `mafft --rnaalifold ...` produced ON THIS MACHINE: the installed Debian mafft 7.505-1 (script and `tbfast` of the same source as the extended build of the Docker package) and `RNAalifold` of the
apt package vienna-rna 2.5.1 (multiverse) - the SAME versions as in the archive of run 2 of handoff/A2-19c-rna (the 25 `--rnaalifold` sets of that archive were run here too: return code, standard output and standard error equal byte for byte, WP83),
so the sandbox oracle is acceptable. 100 new RNA families with stems (deterministic, seed 83083) x 10 option lines (`--rnaalifold` with --globalpair / --localpair / --genafpair / --maxiterate / --retree / --clustalout / --thread / --op / --ep),
each run twice (both runs must be equal), written in the format of every mafft oracle (see convert_rna_returned.py) to tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/mafft_rnaalifold_oracle.tsv.gz (group `rnaalifold`).
The version check refuses to run with other versions. No network access, nothing downloaded.
Usage: python3 scripts/alignment/generate_mafft_rnaalifold_oracle.py"""
import gzip, hashlib, json, os, random, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment')
PAIR = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
OPTS = ['--rnaalifold --globalpair', '--rnaalifold --localpair', '--rnaalifold --genafpair', '--rnaalifold --globalpair --maxiterate 1000', '--rnaalifold --localpair --maxiterate 2',
        '--rnaalifold --globalpair --retree 1', '--rnaalifold --localpair --clustalout', '--rnaalifold --globalpair --thread 1', '--rnaalifold --localpair --op 2.0 --ep 0.2', '--rnaalifold --genafpair --maxiterate 2']


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def family(rng, n, stems):
    base = ''
    for _ in range(stems):
        stem = ''.join(rng.choice('ACGU') for _ in range(rng.randint(4, 7)))
        loop = ''.join(rng.choice('ACGU') for _ in range(rng.randint(3, 6)))
        base += ''.join(rng.choice('ACGU') for _ in range(rng.randint(2, 5))) + stem + loop + ''.join(PAIR[c] for c in reversed(stem))
    out = []
    for _ in range(n):
        s = []
        for c in base:
            x = rng.random()
            if x < 0.03:
                continue
            if x < 0.06:
                s.append(c + rng.choice('ACGU'))
                continue
            s.append(rng.choice('ACGU') if rng.random() < 0.08 else c)
        out.append(''.join(s))
    return out


def run(args, text):
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, 'in.fa'), 'w').write(text)
        p = subprocess.run(['mafft'] + args + ['in.fa'], capture_output=True, cwd=td, env=dict(os.environ, TMPDIR=td), timeout=600)
        files = {f: open(os.path.join(td, f), 'rb').read().decode('latin-1') for f in sorted(os.listdir(td)) if f != 'in.fa'}
        return p.returncode, p.stdout.decode('latin-1'), p.stderr.decode('latin-1'), files


def main():
    v = subprocess.run(['mafft', '--version'], capture_output=True, text=True)
    r = subprocess.run(['RNAalifold', '--version'], capture_output=True, text=True)
    assert 'v7.505' in v.stderr + v.stdout and '2.5.1' in r.stdout + r.stderr, (v.stderr, r.stdout)
    rng = random.Random(83083)
    rows, k = [], 0
    for opts in OPTS:
        for _ in range(10):
            seqs = family(rng, rng.choice([2, 3, 4, 5, 6, 8]), rng.choice([1, 2, 3]))
            text = ''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs))
            args = opts.split()
            a, b = run(args, text), run(args, text)
            assert a == b, 'not deterministic: ' + opts
            rc, out, err, files = a
            cols = ['rnaalifold-s%03d' % k, 'rnaalifold', json.dumps(args), text, '{}', str(rc), out, hashlib.sha256(err.encode('utf-8')).hexdigest(), err, json.dumps(files)]
            rows.append('\t'.join(esc(c) for c in cols))
            k += 1
    with gzip.open(os.path.join(DATA, 'mafft_rnaalifold_oracle.tsv.gz'), 'wt', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(rows) + '\n')
    print('sets', k, file=sys.stderr)


if __name__ == '__main__':
    main()
