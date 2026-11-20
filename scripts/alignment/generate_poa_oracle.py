#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/poa_oracle.tsv (B09 audit A2-26) with the installed POA 2.0 (`poa`, Debian 2.0+20060928-9).

Every row = one run `poa -read_fasta in.fa -clustal c.aln -pir p.pir -po p.po OPTIONS MATRIX` in an empty directory holding the input files; the row records the
arguments, the input files (FASTA, matrix, pair-score file), the exit code, standard output and the three output files. Standard error is not compared: the program prints
uninitialised `%ld` arguments there. Sets are generated here (python `random.Random`, fixed seeds): protein families (BLOSUM80, gap penalties 12 / 6 / 6) and
nucleotide families with a lower-case nucleotide matrix (affine truncation / decay gap model); every option of the port (-do_global, -do_progressive, -read_pairscores,
-fuse_all, -hb, -hbmin, -best, -preserve_seqorder, -tolower / -toupper, -v, -printmatrix).

Escaping of the text columns: backslash -> \\\\, newline -> \\n, tab -> \\t, carriage return -> \\r. Files: name=content joined by the character U+001F.
Run: python3 scripts/alignment/generate_poa_oracle.py
"""
import os, random, subprocess, tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'poa_oracle.tsv')
binary = 'poa'
SEP = '\x1f'
PROT = 'ACDEFGHIKLMNPQRSTVWY'
DNA = 'acgt'

MATFILES = [os.path.join(root, 'src', 'Seqeron', 'Algorithms', 'Seqeron.Genomics.Alignment', 'SequenceAligner.PoaMatrix.cs')]


def blosum80():
    # the matrix of the distribution (verbatim in PoaDefaultMatrix)
    src = open(MATFILES[-1]).read()
    i = src.index('"""\n') + 4
    j = src.index('\n""" +')
    return src[i:j] + '\n'


def blosum80_trunc(b):
    # blosum80_trunc.mat of the distribution: the same scores with a truncated / decaying gap model
    lines = b.split('\n')
    out = []
    for l in lines:
        if l.startswith('GAP-PENALTIES='):
            out += ['GAP-TRUNCATION-LENGTH=10', 'GAP-DECAY-LENGTH=5', 'GAP-PENALTIES=12 6 1']
        else:
            out.append(l)
    return '\n'.join(out)


def dna_matrix(match, mismatch, gaps, extra=''):
    syms = 'acgtn'
    lines = ['# test nucleotide matrix', gaps] + ([extra] if extra else []) + ['   ' + '  '.join(syms)]
    for a in syms:
        row = [a]
        for b in syms:
            if a == 'n' or b == 'n':
                row.append('0')
            else:
                row.append(str(match if a == b else mismatch))
        lines.append(' '.join(row))
    return '\n'.join(lines) + '\n'


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def mutate(rng, base, alpha, sub, indel):
    out = []
    for ch in base:
        r = rng.random()
        if r < indel / 2:
            continue
        if r < indel:
            out.append(ch)
            out.append(rng.choice(alpha))
            continue
        out.append(rng.choice(alpha) if rng.random() < sub else ch)
    return ''.join(out) or base[:1]


def family(rng, n, length, alpha, sub, indel):
    base = ''.join(rng.choice(alpha) for _ in range(length))
    return [mutate(rng, base, alpha, sub, indel) for _ in range(n)]


def fasta(seqs, names=None, titles=None, width=0):
    out = []
    for i, s in enumerate(seqs):
        nm = names[i] if names and i < len(names) else 's%d' % (i + 1)
        t = (' ' + titles[i]) if titles and i < len(titles) and titles[i] else ''
        out.append('>' + nm + t)
        if width:
            out.extend(s[j:j + width] for j in range(0, len(s), width))
        else:
            out.append(s)
    return '\n'.join(out) + '\n'


cases = []  # (id, kind, args, files)


def add(kind, args, files, matrix='blosum80.mat'):
    cases.append(('Q%04d' % (len(cases) + 1), kind, args, files, matrix))


rng = random.Random(20261121)
B80 = blosum80()
FILES_P = {'blosum80.mat': B80}


def prot(nmin=2, nmax=8, lmin=10, lmax=90):
    n = rng.randrange(nmin, nmax + 1)
    return family(rng, n, rng.randrange(lmin, lmax + 1), PROT, rng.choice([0.05, 0.2, 0.4]), rng.choice([0.0, 0.03, 0.1]))


def nuc(nmin=2, nmax=8, lmin=15, lmax=100):
    n = rng.randrange(nmin, nmax + 1)
    return family(rng, n, rng.randrange(lmin, lmax + 1), DNA, rng.choice([0.05, 0.15, 0.3]), rng.choice([0.0, 0.04, 0.1]))


def pf(seqs, **kw):
    return dict(FILES_P, **{'in.fa': fasta(seqs, **kw)})


# A. protein, every mode
for i in range(40):
    add('default', [], pf(prot()))
for i in range(30):
    add('global', ['-do_global'], pf(prot()))
for i in range(30):
    add('progressive', ['-do_progressive'], pf(prot()))
for i in range(20):
    add('progressive', ['-do_progressive', '-do_global'], pf(prot(3, 8)))
for i in range(16):
    add('fuse', ['-fuse_all'] + (['-do_global'] if i % 2 else []) + (['-do_progressive'] if i % 4 >= 2 else []), pf(prot(3, 8)))
for i in range(24):
    add('hb', ['-hb'] + (['-hbmin', ['0.5', '0.8', '0.95', '0.3'][i % 4]] if i % 3 else []) + (['-do_global'] if i % 2 else []), pf(prot(4, 10, 20, 100)))
for i in range(10):
    add('hb', ['-hb', '-best'], pf(prot(4, 9, 20, 90)))
for i in range(12):
    add('order', ['-do_progressive', '-preserve_seqorder'], pf(prot(3, 9)))
for i in range(6):
    add('order', ['-preserve_seqorder'], pf(prot(3, 6)))
for i in range(6):
    add('order', ['-do_progressive', '-preserve_seqorder', '-hb'], pf(prot(4, 8, 20, 80)))
# B. nucleotides, three matrices
mats = {
    'n1.mat': dna_matrix(5, -4, 'GAP-PENALTIES=8 6 6'),
    'n2.mat': dna_matrix(4, -3, 'GAP-PENALTIES=10 4 1', 'GAP-TRUNCATION-LENGTH=6'),
    'n3.mat': dna_matrix(3, -2, 'GAP-PENALTIES=7 3 0', 'GAP-DECAY-LENGTH=4'),
    'n4.mat': dna_matrix(5, -4, 'GAP-PENALTIES=12 6 6', 'GAP-PENALTIES-X=6 3 3'),
}
for mname, mtxt in mats.items():
    for i in range(14):
        a = []
        if i % 2: a.append('-do_global')
        if i % 4 >= 2: a.append('-do_progressive')
        if i % 5 == 4: a.append('-fuse_all')
        if i % 7 == 6: a.append('-hb')
        add('nuc', a, {mname: mtxt, 'in.fa': fasta(nuc())}, matrix=mname)
trunc = blosum80_trunc(B80)
for i in range(16):
    add('trunc', (['-do_global'] if i % 2 else []) + (['-do_progressive'] if i % 4 >= 2 else []), {'trunc.mat': trunc, 'in.fa': fasta(prot(3, 8))}, matrix='trunc.mat')
# C. case switching, unknown letters, odd FASTA
for i in range(8):
    seqs = prot(3, 6)
    seqs = [s if j % 2 else s.lower() for j, s in enumerate(seqs)]
    add('case', [], pf(seqs))
    add('case', ['-toupper'], pf(seqs))
    add('case', ['-tolower'], {'n1.mat': mats['n1.mat'], 'in.fa': fasta([s.upper() for s in nuc(3, 5)])}, matrix='n1.mat')
for i in range(10):
    seqs = [list(s) for s in prot(3, 6)]
    for s in seqs:
        for _ in range(rng.randrange(0, 3)):
            s.insert(rng.randrange(len(s) + 1), rng.choice('XBZ*-.@1'))
    add('unknown', [], pf([''.join(s) for s in seqs]))
for i in range(6):
    seqs = prot(3, 5)
    add('format', [], pf(seqs, width=rng.choice([10, 60])))
add('format', [], {'blosum80.mat': B80, 'in.fa': fasta(prot(3, 4)).replace('\n', '\r\n')})
add('format', [], {'blosum80.mat': B80, 'in.fa': fasta(prot(3, 4)).rstrip('\n')})
add('format', [], {'blosum80.mat': B80, 'in.fa': '# a comment line\n' + fasta(prot(3, 4), titles=['first title here', 'second', ''])})
add('format', [], {'blosum80.mat': B80, 'in.fa': '\n\n' + fasta(prot(3, 4), names=['a_very_long_sequence_name_exceeding_thirty_two_characters', 'b', 'c'], titles=['t1', 't2 with words', '']) + '\n'})
add('format', [], {'blosum80.mat': B80, 'in.fa': '>s1 x\nMKV LAAG\n*ignored\nIVG\n>s2\nMKVLAAGIVG\n>s3 y z\n\nMKVLAGIVG\n'})
add('format', [], {'blosum80.mat': B80, 'in.fa': '>s1\nMKVLAAGIVG\n'})
add('format', [], {'blosum80.mat': B80, 'in.fa': '>s1\nMKVLAAGIVG\n>s2\nMKVLAAGIVG\n'})
add('format', [], {'blosum80.mat': B80, 'in.fa': '>s1\nM\n>s2\nMK\n>s3\nK\n'})
add('format', [], {'blosum80.mat': B80, 'in.fa': '>s1\nMKVLAAGIVG\n>empty\n>s3\nMKVLAGIVG\n'})
add('format', [], {'blosum80.mat': B80, 'in.fa': '>s1\nMKVLAAGIVG\n#end of cluster\n>s2\nMKVLAAGIVG\n'})
# D. pair score files
for i in range(14):
    seqs = prot(3, 7)
    names = ['q%d' % (k + 1) for k in range(len(seqs))]
    lines = []
    for a in range(len(seqs)):
        for b in range(a):
            lines.append('%s\t%s\t%.2f' % (names[a], names[b], rng.uniform(-20, 200)))
    rng.shuffle(lines)
    add('pairscores', ['-read_pairscores', 'pairs.txt'] + (['-do_global'] if i % 2 else []), dict(pf(seqs, names=names), **{'pairs.txt': '\n'.join(lines) + '\n'}))
for i in range(4):
    seqs = prot(4, 6)
    names = ['q%d' % (k + 1) for k in range(len(seqs))]
    lines = ['%s %s %d' % (names[2], names[0], 50), '%s %s %d' % (names[3], names[1], 40)]
    add('pairscores', ['-read_pairscores', 'pairs.txt', '-preserve_seqorder'], dict(pf(seqs, names=names), **{'pairs.txt': '\n'.join(lines) + '\n'}))
add('pairscores', ['-read_pairscores', 'pairs.txt'], dict(pf(prot(3, 3)), **{'pairs.txt': 's1 nosuch 5\n'}))
# E. verbose, printmatrix
for i in range(6):
    add('verbose', ['-v'], pf(prot(2, 4)))
add('verbose', ['-v'], {'trunc.mat': trunc, 'in.fa': fasta(prot(2, 3))}, matrix='trunc.mat')
add('verbose', ['-v'], {'n3.mat': mats['n3.mat'], 'in.fa': fasta(nuc(2, 3))}, matrix='n3.mat')
add('verbose', ['-v'], {'n4.mat': mats['n4.mat'], 'in.fa': fasta(nuc(2, 3))}, matrix='n4.mat')
add('matrix', ['-printmatrix', 'ARND'], pf(prot(2, 3)))
add('matrix', ['-printmatrix', 'WYV?'], pf(prot(2, 3)))
add('matrix', ['-printmatrix', 'acgt'], {'n1.mat': mats['n1.mat'], 'in.fa': fasta(nuc(2, 3))}, matrix='n1.mat')
# F. long sets
for i in range(6):
    add('long', [] if i % 2 else ['-do_global'], pf(family(rng, rng.randrange(3, 6), rng.randrange(150, 300), PROT, 0.15, 0.04)))
# G. errors
add('error', [], {'blosum80.mat': B80})
add('error', ['-subset', 'x.txt'], pf(prot(2, 3)))
add('error', ['-hb'], {'blosum80.mat': B80, 'in.fa': '>s1\n'})
add('error', [], {'in.fa': fasta(prot(2, 3))}, matrix='nomatrix.mat')
add('error', [], {'bad.mat': 'A C\nA 1 2\nB 3 4\n', 'in.fa': fasta(prot(2, 3))}, matrix='bad.mat')
add('error', [], {'bad.mat': 'A C\nA 1\nC 3 4\n', 'in.fa': fasta(prot(2, 3))}, matrix='bad.mat')

rows = []
with tempfile.TemporaryDirectory() as d:
    for ci, (cid, kind, args, files, matrix) in enumerate(cases):
        run = os.path.join(d, 'run%d' % ci)
        os.makedirs(run)
        for nm, txt in files.items():
            with open(os.path.join(run, nm), 'w', newline='') as f:
                f.write(txt)
        before = set(os.listdir(run))
        if kind == 'error' and not args and 'in.fa' not in files:
            cmd = [binary]
            argl = []
        else:
            argl = (['-read_fasta', 'in.fa'] if '-read_fasta' not in args else []) + ['-clustal', 'c.aln', '-pir', 'p.pir', '-po', 'p.po'] + args + [matrix]
            cmd = [binary] + argl
        p = subprocess.run(cmd, cwd=run, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=300, errors='surrogateescape')
        outs = {}
        for nm in sorted(os.listdir(run)):
            if nm not in before:
                outs[nm] = open(os.path.join(run, nm), newline='', errors='surrogateescape').read()
        rows.append((cid, kind, argl, files, p.returncode, p.stdout, outs))

ver = subprocess.run(['dpkg-query', '-W', binary], capture_output=True, text=True).stdout.strip().replace('\t', ' ')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='\n') as fh:
    fh.write('# ' + ver + '; uname -m ' + os.uname().machine + '\n')
    fh.write('id\tkind\targs\tfiles\trc\tstdout\toutfiles\n')
    for cid, kind, args, files, rc, so, outs in rows:
        fh.write('\t'.join([cid, kind, esc(SEP.join(args)) if args else '-', esc(SEP.join(k + '=' + v for k, v in files.items())), str(rc), esc(so),
                            esc(SEP.join(k + '=' + v for k, v in outs.items()))]) + '\n')
print(len(rows), 'rows;', sum(1 for r in rows if r[4] != 0), 'with a non-zero exit;', ver)
