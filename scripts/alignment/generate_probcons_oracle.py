#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/probcons_oracle.tsv (B09 audit A2-23) with the installed ProbCons 1.12 (`probcons`, Debian 1.12-14).

Every row = one run `probcons ARGS` in an empty directory that holds the input files; the row records the arguments, the input files, the exit code,
standard output, standard error (the whole text: heading, parameter set, tree, progress lines, errors) and every file the program wrote (pairwise
alignments, -annot file, trained parameter file). The sets are generated here (python `random.Random`, fixed seeds): families of protein sequences
with substitutions and indels, lower-case runs, ambiguity codes (X, B, Z, J, U, O), '*', gap characters; every option of the program (-c, -ir, -pre,
-a, -clustalw, -annot, -pairs, -viterbi, -co, -go, -ge, -v, -t, -e, -p), MSF / Clustal / odd FASTA inputs and the program's error exits.

Escaping of the text columns: backslash -> \\\\, newline -> \\n, tab -> \\t, carriage return -> \\r. Files: name=content joined by the character U+001F.
Run: python3 scripts/alignment/generate_probcons_oracle.py
"""
import os, random, subprocess, tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'probcons_oracle.tsv')
binary = 'probcons'
SEP = '\x1f'
PROT = 'ACDEFGHIKLMNPQRSTVWY'


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def mutate(rng, base, sub, indel):
    out = []
    for ch in base:
        r = rng.random()
        if r < indel / 2:
            continue
        if r < indel:
            out.append(ch)
            out.append(rng.choice(PROT))
            continue
        out.append(rng.choice(PROT) if rng.random() < sub else ch)
    return ''.join(out) or base[:1]


def family(rng, n, length, sub, indel):
    base = ''.join(rng.choice(PROT) for _ in range(length))
    return [mutate(rng, base, sub, indel) for _ in range(n)]


def decorate(rng, seqs):
    res = []
    for s in seqs:
        s = list(s)
        if rng.random() < 0.25 and len(s) > 4:
            a = rng.randrange(len(s)); b = min(len(s), a + rng.randrange(1, 8))
            for i in range(a, b):
                s[i] = s[i].lower()
        if rng.random() < 0.25:
            for _ in range(rng.randrange(1, 3)):
                s[rng.randrange(len(s))] = rng.choice('XBZJUO*')
        if rng.random() < 0.15:
            for _ in range(rng.randrange(1, 3)):
                s.insert(rng.randrange(len(s) + 1), rng.choice('-.'))
        res.append(''.join(s))
    return res


def fasta(seqs, names=None, width=0):
    out = []
    for i, s in enumerate(seqs):
        nm = names[i] if names and i < len(names) else 's%d' % (i + 1)
        out.append('>' + nm)
        if width:
            out.extend(s[j:j + width] for j in range(0, len(s), width))
        else:
            out.append(s)
    return '\n'.join(out) + '\n'


cases = []  # (id, kind, args, {file: text})


def add(kind, args, files):
    cases.append(('P%04d' % (len(cases) + 1), kind, args, files))


rng = random.Random(20261120)


def fam(nmin=2, nmax=9, lmin=10, lmax=110):
    n = rng.randrange(nmin, nmax + 1)
    return decorate(rng, family(rng, n, rng.randrange(lmin, lmax + 1), rng.choice([0.05, 0.2, 0.4]), rng.choice([0.0, 0.03, 0.1])))


# A. defaults
for i in range(50):
    add('default', [], {'in.fa': fasta(fam())})
# B. -c / -ir
for i in range(28):
    add('c', ['-c', str(i % 6)], {'in.fa': fasta(fam(3, 8, 10, 80))})
for i in range(28):
    add('ir', ['-ir', str([0, 1, 2, 5, 10, 30, 200][i % 7])], {'in.fa': fasta(fam(3, 8, 10, 80))})
for i in range(10):
    add('ir', ['-ir', '1000', '-c', '1'], {'in.fa': fasta(fam(4, 7, 10, 40))})
# C. -clustalw, -a, -annot
for i in range(20):
    add('clustalw', ['-clustalw'], {'in.fa': fasta(fam(2, 8, 10, 130), names=None)})
for i in range(14):
    seqs = fam(3, 8, 10, 70)
    add('order', ['-a'] + (['-clustalw'] if i % 2 else []), {'in.fa': fasta(seqs)})
for i in range(20):
    add('annot', ['-annot', 'ann.txt'] + (['-a'] if i % 4 == 3 else []), {'in.fa': fasta(fam(2, 8, 10, 90))})
# D. all pairs
for i in range(10):
    add('pairs', ['-pairs'], {'in.fa': fasta(fam(2, 6, 10, 70))})
for i in range(10):
    add('viterbi', ['-viterbi'], {'in.fa': fasta(fam(2, 6, 10, 70))})
for i in range(8):
    add('pairs', ['-pairs', '-co', '0.3', '-c', str(i % 3)] + (['-clustalw'] if i % 2 else []), {'in.fa': fasta(fam(2, 5, 10, 60))})
for i in range(4):
    add('viterbi', ['-viterbi', '-clustalw'], {'in.fa': fasta(fam(2, 4, 10, 60))})
# E. cosmetic gap penalties and the cutoff. The program allocates `openingPenalty2` with swapped dimensions (Main.cc / ProbabilisticModel.h
# `VVF openingPenalty2 (numSeqs1+1, VF (numSeqs2+1))` filled for `i <= numSeqs2, j <= numSeqs1`): a heap overflow whenever the two profiles differ in size (rc -11 /
# -6 or garbage). -go / -ge are therefore only run where every merge has two profiles of equal size: two sequences (1 + 1) and four sequences with -ir 0 and a balanced tree.
import re as _re
GOGE = [('-1', '0'), ('0', '-1'), ('-2', '-1'), ('-0.5', '-0.5'), ('-5', '-2'), ('-10', '-10'), ('-1', '-0.5'), ('-3', '0'), ('-0.1', '-0.01'), ('-7', '-3')]
for go, ge in GOGE:
    for _ in range(3):
        add('penalties', ['-go', go, '-ge', ge], {'in.fa': fasta(decorate(rng, family(rng, 2, rng.randrange(10, 90), 0.2, 0.05)))})
    add('penalties', ['-go', go], {'in.fa': fasta(decorate(rng, family(rng, 2, rng.randrange(10, 90), 0.2, 0.05)))})
    add('penalties', ['-ge', ge, '-co', '0.1'], {'in.fa': fasta(decorate(rng, family(rng, 2, rng.randrange(10, 90), 0.2, 0.05)))})
balanced = 0
tries = 0
with tempfile.TemporaryDirectory() as dd:
    while balanced < 24 and tries < 400:
        tries += 1
        text = fasta(decorate(rng, family(rng, 4, rng.randrange(10, 70), 0.2, 0.05)))
        with open(os.path.join(dd, 'x.fa'), 'w', newline='') as f:
            f.write(text)
        pr = subprocess.run([binary, '-ir', '0', 'x.fa'], cwd=dd, capture_output=True, text=True)
        if _re.search(r'Alignment tree: \(\(\S+ \S+\) \(\S+ \S+\)\)\n', pr.stderr):
            go, ge = GOGE[balanced % len(GOGE)]
            add('penalties', ['-ir', '0', '-go', go, '-ge', ge], {'in.fa': text})
            balanced += 1
for co in ('0.1', '0.5', '1', '0', '0.05'):
    for _ in range(3):
        add('cutoff', ['-co', co], {'in.fa': fasta(fam(3, 8, 10, 80))})
for _ in range(4):
    add('cutoff', ['-co', '0.2', '-c', '0', '-ir', '5'], {'in.fa': fasta(fam(3, 8, 10, 80))})
# F. pre-training, training, parameter files
for i in range(14):
    add('pre', ['-pre', str(1 + i % 4)], {'in.fa': fasta(fam(3, 7, 10, 70))})
train_sets = [fam(4, 7, 30, 90) for _ in range(3)]
add('train', ['-t', 'params.out', 'a.fa'], {'a.fa': fasta(train_sets[0])})
add('train', ['-t', 'params.out', 'a.fa', 'b.fa'], {'a.fa': fasta(train_sets[0]), 'b.fa': fasta(train_sets[1])})
add('train', ['-t', 'params.out', '-e', 'a.fa'], {'a.fa': fasta(train_sets[0])})
add('train', ['-t', 'params.out', '-e', 'a.fa', 'b.fa', 'c.fa'], {'a.fa': fasta(train_sets[0]), 'b.fa': fasta(train_sets[1]), 'c.fa': fasta(train_sets[2])})
for i in range(6):
    add('train', ['-t', 'params.out'] + (['-e'] if i % 2 else []) + ['a.fa'], {'a.fa': fasta(fam(3, 6, 20, 80))})
add('train', ['-t', 'params.out', '-v', 'a.fa'], {'a.fa': fasta(train_sets[1])})

# G. verbose
for i in range(10):
    add('verbose', ['-v'] + (['-pre', '1'] if i % 5 == 4 else []) + (['-go', '-1', '-ge', '-1'] if i % 3 == 2 else []),
        {'in.fa': fasta(decorate(rng, family(rng, 2, 30, 0.2, 0.05)) if i % 3 == 2 else fam(2, 6, 10, 60))})
add('verbose', ['-v', '-pairs'], {'in.fa': fasta(fam(3, 4, 10, 50))})
add('verbose', ['-v', '-viterbi'], {'in.fa': fasta(fam(3, 4, 10, 50))})

# H. input formats
for i in range(8):
    seqs = fam(3, 6, 20, 90)
    add('format', [], {'in.fa': fasta(seqs, width=rng.choice([10, 60, 70]))})
add('format', [], {'in.fa': fasta(fam(3, 5, 20, 60)).replace('\n', '\r\n')})                     # CRLF
add('format', [], {'in.fa': fasta(fam(3, 5, 20, 60)).rstrip('\n')})                              # no trailing newline
add('format', [], {'in.fa': '\n\n' + fasta(fam(3, 5, 20, 60), names=['a b c', ' x  ', 'seq|3', 'd', 'e'], ) + '\n\n'})
add('format', [], {'in.fa': '>s1\nMKV LAAG\nIVG\n>s2\nMKV.LAAGIVG-\n>s3 \nMKVLAGIVG\n'})        # spaces, '.', '-' inside the sequence
add('format', [], {'in.fa': '>s1\nMKVLAAGIVG\n>empty\n>s3\nMKVLAGIVG\n'})                        # an empty sequence stops the reader
add('format', [], {'in.fa': '>empty\n>s2\nMKVLAAGIVG\n>s3\nMKVLAGIVG\n'})                        # empty first sequence: MSF branch
add('format', [], {'in.fa': '>s1\nMKVLAAGIVG\n'})
add('format', [], {'in.fa': '>s1\nMKVLAAGIVG\n>s2\nMKVLAAGIVG\n'})
add('format', [], {'in.fa': '>s1\nM\n>s2\nMK\n>s3\nK\n'})
add('format', [], {'in.fa': '>s1\nMKVLAAGIVG\n>s2\nMKV1LAGIVG\n'})                               # unknown character
for kind in ('msf', 'clustal', 'clustal2', 'noheader'):
    seqs = family(rng, 4, 70, 0.1, 0.05)
    # an alignment (gaps inserted)
    rows = []
    for s in seqs:
        rows.append(s)  # no gap characters: the program keeps them in MSF / Clustal input and then overflows a heap buffer (BuildPosterior)
    L = max(len(r) for r in rows)
    names = ['alpha', 'beta', 'gamma', 'delta']
    if kind == 'msf':
        text = '!!AA_MULTIPLE_ALIGNMENT 1.0\n\n  MSF: %d  Type: P  Check: 0 ..\n\n' % L
        for nm, r in zip(names, rows):
            text += ' Name: %s oo  Len: %d  Check: 0  Weight: 1.00\n' % (nm, L)
        text += '\n//\n\n'
        for j in range(0, L, 50):
            for nm, r in zip(names, rows):
                if j < len(r):
                    text += '%-12s%s\n' % (nm, ' '.join(r[j + k:j + k + 10] for k in range(0, 50, 10)))
            text += '\n'
    elif kind in ('clustal', 'clustal2'):
        text = 'CLUSTAL W (1.83) multiple sequence alignment\n\n' if kind == 'clustal' else 'PROBCONS version 1.12 multiple sequence alignment\n\n'
        for j in range(0, L, 60):
            for nm, r in zip(names, rows):
                if j < len(r):
                    text += '%-14s%s\n' % (nm, r[j:j + 60])
            text += '%-14s%s\n\n' % ('', '*' * min(60, L - j))
    else:
        text = 'some title\n//\n\n'
        for nm, r in zip(names, rows):
            text += '%s %s\n' % (nm, r)
        text += '\n'
    add('format', [], {'in.fa': text})
    add('format', ['-clustalw', '-ir', '3'], {'in.fa': text})
# K. several files
for i in range(6):
    add('files', ['a.fa', 'b.fa'], {'a.fa': fasta(fam(2, 4, 10, 60)), 'b.fa': fasta(fam(1, 3, 10, 60), names=['x1', 'x2', 'x3'])})
# L. longer sequences
for i in range(6):
    add('long', ['-ir', '3'], {'in.fa': fasta(decorate(rng, family(rng, rng.randrange(3, 5), rng.randrange(200, 400), 0.15, 0.05)))})
# M. identical / similar / tiny sets
add('tiny', [], {'in.fa': fasta(['ACDEFGHIK'] * 4)})
add('tiny', [], {'in.fa': fasta(['A', 'A', 'A'])})
add('tiny', [], {'in.fa': fasta(['AC', 'CA', 'AC', 'GG'])})
add('tiny', [], {'in.fa': fasta(['MKV', 'MKVLAAAAAAAAAAAAAA'])})
add('tiny', [], {'in.fa': fasta(['XXXXBBBZZZ', 'XXXXBBBZZZ', 'XXBBZZ'])})
add('tiny', [], {'in.fa': fasta(['mkvlaagivg', 'MKVLAAGIVG', 'MkVlAaGiVg'])})
add('tiny', [], {'in.fa': fasta(['MKV*LAAG', 'MKVLAAG*', '*MKVLAAG'])})
# N. error exits
for args in (['-q'], ['-c', '9'], ['-c', 'x'], ['-c'], ['-ir', '1001'], ['-ir', '-1'], ['-pre', '21'], ['-go', '1'], ['-ge', '0.5'], ['-co', '2'], ['-co', '-0.1'], ['-go', 'abc'],
             ['-e'], ['-t'], ['-p'], ['-annot'], ['-c', '0x2'], ['-c', '010'], ['-c', '2abc'], ['-p', 'nofile.txt']):
    add('error', args + ['in.fa'], {'in.fa': fasta(fam(3, 4, 10, 40))})
add('error', ['missing.fa'], {})
add('error', [], {})

rows = []
with tempfile.TemporaryDirectory() as d:
    # parameter files for -p: trained with the program itself (then perturbed copies of the default format)
    base = os.path.join(d, 'base'); os.makedirs(base)
    with open(os.path.join(base, 'a.fa'), 'w', newline='') as f:
        f.write(fasta(fam(5, 6, 40, 90)))
    subprocess.run([binary, '-t', 'p1.txt', '-e', 'a.fa'], cwd=base, capture_output=True)
    subprocess.run([binary, '-t', 'p2.txt', 'a.fa'], cwd=base, capture_output=True)
    pf = {n: open(os.path.join(base, n)).read() for n in ('p1.txt', 'p2.txt') if os.path.exists(os.path.join(base, n))}
    for n, txt in pf.items():
        for i in range(4):
            add('paramfile', ['-p', n] + (['-clustalw'] if i == 3 else []), {'in.fa': fasta(fam(3, 7, 10, 70)), n: txt})
    for ci, (cid, kind, args, files) in enumerate(cases):
        run = os.path.join(d, 'run%d' % ci); os.makedirs(run)
        for nm, txt in files.items():
            with open(os.path.join(run, nm), 'w', newline='') as f:
                f.write(txt)
        before = set(os.listdir(run))
        cmd = [binary] + args + (['in.fa'] if 'in.fa' in files and 'in.fa' not in args else [])
        p = subprocess.run(cmd, cwd=run, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=300, errors='surrogateescape')
        outs = {}
        for nm in sorted(os.listdir(run)):
            if nm not in before or open(os.path.join(run, nm), newline='', errors='surrogateescape').read() != files.get(nm):
                outs[nm] = open(os.path.join(run, nm), newline='', errors='surrogateescape').read()
        rows.append((cid, kind, cmd[1:], files, p.returncode, p.stdout, p.stderr, outs))

ver = subprocess.run(['dpkg-query', '-W', binary], capture_output=True, text=True).stdout.strip().replace('\t', ' ')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='\n') as fh:
    fh.write('# ' + ver + '; uname -m ' + os.uname().machine + '\n')
    fh.write('id\tkind\targs\tfiles\trc\tstdout\tstderr\toutfiles\n')
    for cid, kind, args, files, rc, so, se, outs in rows:
        fh.write('\t'.join([cid, kind, esc(SEP.join(args)) if args else '-',
                            esc(SEP.join(k + '=' + v for k, v in files.items())), str(rc), esc(so), esc(se),
                            esc(SEP.join(k + '=' + v for k, v in outs.items()))]) + '\n')
print(len(rows), 'rows;', sum(1 for r in rows if r[4] != 0), 'with a non-zero exit;', ver)
