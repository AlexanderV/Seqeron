#!/usr/bin/env python3
"""B09 WP86 (item A2-19c-rna, F377): the sandbox oracle of the RNA step of `mafft --contrafold ...` and `mafft --mccaskill ...` (tbfast with -r -c -o: `foldrna`, `Lalignmm_hmout`, the importance matrix of the pair library) with a NON-ZERO structure weight, produced ON THIS MACHINE:
the installed Debian mafft 7.505-1 (the script and ALL its programs, among them the real `contrafoldwrap` and `tbfast`) run with MAFFT_BINARIES pointing at a private directory of symbolic links to the Debian programs plus ONE program of our own, `contrafold`, which stands for
CONTRAfold: it answers `contrafold predict IN --posteriors 0.01 OUT` with a table of base-pairing probabilities that this script made up (random stems, six significant digits, the BPSEQ-like text of CONTRAfold 2.02) and records in canned.json;
and (group `mccaskill-canned`) ONE program that stands for MXSCARNA, `mxscarnamod -m -writebpp IN` (the real `mccaskillwrap` of Debian runs it): it prints made-up pair lines "left right p" (17 significant digits, some below 0.01, some with left = right).
So the oracle checks everything of `mafft --contrafold` except CONTRAfold itself (the true CONTRAfold is checked by the 16 sets of the Docker package A2-19c-rna): contrafoldwrap (the reading of the posteriors, hat4), the progressive stage with the structure term, the options --weightr / --weightm / --rop.
Format of every mafft oracle (see convert_rna_returned.py), group `contrafold-canned` (and `contrafold-canned-iterate` for --maxiterate, which the port refuses); the `files` column holds {"canned": <canned.json text>}; the text of the private directory in standard error is replaced by /usr/lib/mafft/lib/mafft.
Deterministic (seed 86086), each run twice (equal). No network access, nothing downloaded or compiled.
Usage: python3 scripts/alignment/generate_mafft_contrafold_oracle.py"""
import gzip, hashlib, json, os, random, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment')
DEBIAN = '/usr/lib/mafft/lib/mafft'
PAIR = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G'}
OPTS = ['--contrafold --globalpair --weightr 8.0 --weightm 2.0', '--contrafold --localpair --weightr 8.0 --weightm 2.0', '--contrafold --genafpair --weightr 8.0 --weightm 2.0',
        '--contrafold --globalpair --weightr 4.0', '--contrafold --localpair --weightr 12.0 --weightm 1.5', '--contrafold --globalpair --weightr 8.0 --rop -2.0',
        '--contrafold --localpair --weightr 8.0 --op 2.0 --ep 0.2', '--contrafold --globalpair --weightr 8.0 --weightm 2.0 --clustalout', '--contrafold --globalpair --weightr 8.0 --retree 1',
        '--contrafold --localpair --weightr 8.0 --thread 1', '--contrafold --globalpair --weightr 0.0 --weightm 2.0', '--contrafold --localpair --weightr 8.0 --weightm 2.0 --reorder']
MCC = ['--mccaskill --globalpair --weightr 8.0 --weightm 2.0', '--mccaskill --localpair --weightr 8.0 --weightm 2.0', '--mccaskill --genafpair --weightr 4.0', '--mccaskill --globalpair --weightr 8.0 --clustalout',
       '--mccaskill --localpair --weightr 12.0 --thread 1', '--mccaskill --globalpair --weightr 8.0 --retree 1', '--mccaskill --globalpair']
MCC_ITER = ['--mccaskill --globalpair --maxiterate 2 --weightr 8.0 --weightm 2.0']
ITER = ['--contrafold --globalpair --maxiterate 2 --weightr 8.0 --weightm 2.0', '--contrafold --localpair --maxiterate 1000 --weightr 8.0 --weightm 2.0']
FAKE_MX = '''#!/usr/bin/python3
import json, os, sys
a = sys.argv
if a[1:3] != ['-m', '-writebpp'] or len(a) != 4:
    sys.stderr.write('fake mxscarnamod: unexpected arguments %r\\n' % (a[1:],)); sys.exit(2)
seq = ''.join(l.strip() for l in open(a[3]).read().split('\\n')[1:]).upper()
canned = json.load(open(os.environ['CF_CANNED']))
sys.stdout.write(canned[seq])
'''
FAKE = '''#!/usr/bin/python3
import json, os, sys
a = sys.argv
if len(a) != 6 or a[1] != 'predict' or a[3] != '--posteriors' or a[4] != '0.01':
    sys.stderr.write('fake contrafold: unexpected arguments %r\\n' % (a[1:],)); sys.exit(2)
seq = ''.join(l.strip() for l in open(a[2]).read().split('\\n')[1:])
canned = json.load(open(os.environ['CF_CANNED']))
open(a[5], 'w').write(canned[seq])
'''


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def family(rng, n, stems):
    base = ''
    for _ in range(stems):
        stem = ''.join(rng.choice('ACGU') for _ in range(rng.randint(4, 8)))
        loop = ''.join(rng.choice('ACGU') for _ in range(rng.randint(3, 7)))
        base += ''.join(rng.choice('ACGU') for _ in range(rng.randint(2, 6))) + stem + loop + ''.join(PAIR[c] for c in reversed(stem))
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


def posteriors(rng, seq):
    """A table like the output of `contrafold predict --posteriors 0.01`: line i = 'i letter j:p ...' (j > i, p >= 0.01, six significant digits)."""
    L = len(seq)
    pairs = {}
    for _ in range(rng.randint(1, 3)):
        k = rng.randint(3, 7)
        if L < 2 * k + 5:
            continue
        a = rng.randint(1, L - 2 * k - 4)
        b = rng.randint(a + k + 3, L - k + 1)
        p0 = rng.uniform(0.3, 0.99)
        for t in range(k):
            if a + t < b - t - 2:
                pairs[(a + t, b - t)] = max(0.011, p0 * rng.uniform(0.8, 1.0))
    for _ in range(rng.randint(2, 8)):
        i = rng.randint(1, max(1, L - 4))
        j = rng.randint(i + 3, L) if i + 3 <= L else None
        if j:
            pairs[(i, j)] = rng.uniform(0.011, 0.2)
    lines = []
    for i in range(1, L + 1):
        row = ' '.join('%d:%s' % (j, '%.6g' % p) for (a, j), p in sorted(pairs.items(), key=lambda kv: kv[0][1]) if a == i)
        lines.append('%d %s%s' % (i, seq[i - 1], ' ' + row if row else ''))
    return '\n'.join(lines) + '\n'


def bpp_text(rng, seq):
    """What `mxscarnamod -m -writebpp` could print: a header and lines 'left right p' (zero-based, 17 significant digits)."""
    L = len(seq)
    lines = ['>in']
    for _ in range(rng.randint(3, max(4, L // 3))):
        i = rng.randint(0, L - 4)
        j = rng.randint(i + 3, L - 1) if i + 3 <= L - 1 else None
        if j is not None:
            lines.append('%d %d %s' % (i, j, '%.17g' % (rng.uniform(0.3, 0.99) if rng.random() < 0.7 else rng.uniform(0.001, 0.2))))
    lines.append('%d %d %s' % (0, 0, '%.17g' % rng.uniform(0.1, 0.9)))
    return '\n'.join(lines) + '\n'


def run(libdir, canned_path, args, text):
    with tempfile.TemporaryDirectory() as td:
        open(os.path.join(td, 'in.fa'), 'w').write(text)
        env = dict(os.environ, TMPDIR=td, MAFFT_BINARIES=libdir, CF_CANNED=canned_path)
        p = subprocess.run(['mafft'] + args + ['in.fa'], capture_output=True, cwd=td, env=env, timeout=600)
        files = {f: open(os.path.join(td, f), 'rb').read().decode('latin-1') for f in sorted(os.listdir(td)) if f != 'in.fa'}
        return p.returncode, p.stdout.decode('latin-1'), p.stderr.decode('latin-1').replace(libdir, DEBIAN), files


def main():
    v = subprocess.run(['mafft', '--version'], capture_output=True, text=True)
    assert 'v7.505' in v.stderr + v.stdout, v.stderr
    libdir = tempfile.mkdtemp(prefix='mafftlib')
    for f in os.listdir(DEBIAN):
        os.symlink(os.path.join(DEBIAN, f), os.path.join(libdir, f))
    fake = os.path.join(libdir, 'contrafold')
    if os.path.lexists(fake):
        os.unlink(fake)
    open(fake, 'w').write(FAKE)
    os.chmod(fake, 0o755)
    fake_mx = os.path.join(libdir, 'mxscarnamod')
    if os.path.lexists(fake_mx):
        os.unlink(fake_mx)
    open(fake_mx, 'w').write(FAKE_MX)
    os.chmod(fake_mx, 0o755)
    rng = random.Random(86086)
    rows, k = [], 0
    for group, options, per in (('contrafold-canned', OPTS, 8), ('contrafold-canned-iterate', ITER, 3), ('mccaskill-canned', MCC, 6), ('mccaskill-canned-iterate', MCC_ITER, 3)):
        for opts in options:
            for _ in range(per):
                seqs = family(rng, rng.choice([2, 3, 4, 5, 6, 8]), rng.choice([1, 2, 3]))
                text = ''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs))
                canned = {s.upper(): (bpp_text(rng, s) if group.startswith('mccaskill') else posteriors(rng, s)) for s in seqs}
                cpath = os.path.join(libdir, 'canned.json')
                json.dump(canned, open(cpath, 'w'))
                args = opts.split()
                a, b = run(libdir, cpath, args, text), run(libdir, cpath, args, text)
                assert a == b, 'not deterministic: ' + opts
                rc, out, err, files = a
                cols = ['%s-%03d' % (group, k), group, json.dumps(args), text, json.dumps({'canned': json.dumps(canned, sort_keys=True)}), str(rc), out, hashlib.sha256(err.encode('utf-8')).hexdigest(), err, json.dumps(files)]
                rows.append('\t'.join(esc(c) for c in cols))
                k += 1
    with gzip.open(os.path.join(DATA, 'mafft_contrafold_oracle.tsv.gz'), 'wt', encoding='utf-8', newline='\n') as f:
        f.write('\n'.join(rows) + '\n')
    print('sets', k, file=sys.stderr)


if __name__ == '__main__':
    main()
