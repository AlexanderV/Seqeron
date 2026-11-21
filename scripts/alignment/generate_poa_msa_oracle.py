#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/poa_msa_oracle.tsv (B09 audit A2-26b) with the installed POA 2.0 (`poa`, Debian 2.0+20060928-9): MSA / PO file input.

Every row = one run `poa READ-OPTIONS -clustal c.aln -pir p.pir -po p.po OPTIONS MATRIX` in an empty directory. The input files are made by the program itself from random families
(its PO / PIR / CLUSTAL outputs, with -hb / -best variants) and by hand (CLUSTAL with conservation lines and without header, blocks, long rows (fgets chunks of 510 characters),
PIR with comments, blank lines (the program re-uses the previous row), digits, CRLF, subset / remove lists, a list of files). The row records the arguments, the files, the exit code,
standard output, the lines of standard error that start with `...` (Read / Wrote) and the written files (same escaping as poa_oracle.tsv).
Run: python3 scripts/alignment/generate_poa_msa_oracle.py
"""
import os, random, subprocess, tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'poa_msa_oracle.tsv')
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



rng = random.Random(20261302)
B80 = blosum80()
N1 = dna_matrix(5, -4, 'GAP-PENALTIES=8 6 6')
TRUNC = blosum80_trunc(B80)


def run_poa(args, files, d, stdin=None):
    run = tempfile.mkdtemp(dir=d)
    for nm, txt in files.items():
        with open(os.path.join(run, nm), 'w', newline='', encoding='utf-8') as f:
            f.write(txt)
    before = set(os.listdir(run))
    try:
        p = subprocess.run([binary] + args, cwd=run, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=10, errors='surrogateescape')
    except subprocess.TimeoutExpired:
        # the program loops forever (a PO file / alignment it cannot read): recorded as exit status -999
        p = subprocess.CompletedProcess([binary] + args, -999, '', '')
    outs = {}
    for nm in sorted(os.listdir(run)):
        if nm not in before:
            outs[nm] = open(os.path.join(run, nm), newline='', errors='surrogateescape').read()
    return p, outs


work = tempfile.mkdtemp()
# 1. generated inputs: protein and nucleotide families aligned by the program
gen = []  # (kind, matrix name, matrix text, {po, pir, aln}, names)
for k in range(90):
    nuc = k % 3 == 2
    n = rng.randrange(2, 8)
    if nuc:
        seqs = family(rng, n, rng.randrange(20, 70), DNA, rng.choice([0.05, 0.2]), rng.choice([0.0, 0.05, 0.1]))
        mname, mtxt = 'n1.mat', N1
    else:
        seqs = family(rng, n, rng.randrange(15, 70), PROT, rng.choice([0.05, 0.2, 0.4]), rng.choice([0.0, 0.05, 0.1]))
        mname, mtxt = ('blosum80.mat', B80) if k % 3 else ('trunc.mat', TRUNC)
    names = ['s%d' % (i + 1) for i in range(n)]
    titles = [rng.choice(['', 'first title', 'a b c']) for _ in range(n)]
    extra = rng.choice([[], [], ['-do_global'], ['-hb'], ['-hb', '-best'], ['-do_progressive']])
    p, outs = run_poa(['-read_fasta', 'in.fa', '-clustal', 'c.aln', '-pir', 'p.pir', '-po', 'p.po'] + extra + [mname], {mname: mtxt, 'in.fa': fasta(seqs, names=names, titles=titles)}, work)
    if p.returncode == 0 and 'p.po' in outs:
        gen.append((mname, mtxt, outs, names))

cases = []


def add(kind, args, files, matrix):
    cases.append((kind, args, files, matrix))


OUT = ['-clustal', 'c.aln', '-pir', 'p.pir', '-po', 'p.po']
for i, (mname, mtxt, outs, names) in enumerate(gen):
    fmt = ['p.po', 'p.pir', 'c.aln'][i % 3]
    mats = {mname: mtxt}
    opt = rng.choice([[], [], ['-do_global'], ['-fuse_all'], ['-hb'], ['-preserve_seqorder'], ['-do_progressive']])
    if i % 5 == 0:
        opt += ['-toupper'] if mname != 'n1.mat' else ['-tolower']
    add('single', ['-read_msa', 'm1'] + OUT + opt + [mname], dict(mats, m1=outs[fmt]), mname)
    # two MSAs (different files), or an MSA and FASTA
    o2 = gen[(i * 7 + 3) % len(gen)]
    if o2[0] == mname:
        fmt2 = ['p.pir', 'c.aln', 'p.po'][i % 3]
        add('two', ['-read_msa', 'm1', '-read_msa2', 'm2'] + OUT + opt + [mname], dict(mats, m1=outs[fmt], m2=o2[2][fmt2]), mname)
    seqs = family(rng, 2, 30, DNA if mname == 'n1.mat' else PROT, 0.1, 0.05)
    add('msa+fasta', ['-read_msa', 'm1', '-read_fasta', 'in.fa'] + OUT + opt + [mname], dict(mats, m1=outs[fmt], **{'in.fa': fasta(seqs, names=['x1', 'x2'])}), mname)
    if i % 4 == 0:
        sel = [nm for nm in names if rng.random() < 0.5] or names[:1]
        add('subset', ['-read_msa', 'm1', '-subset', 'sub.txt'] + OUT + opt + [mname], dict(mats, m1=outs[fmt], **{'sub.txt': ''.join('SOURCENAME=%s\n' % s for s in sel)}), mname)
        add('remove', ['-read_msa', 'm1', '-remove', 'sub.txt'] + OUT + opt + [mname], dict(mats, m1=outs[fmt], **{'sub.txt': ''.join('SOURCENAME=%s\n' % s for s in sel)}), mname)
    if i % 6 == 1:
        add('list', ['-read_msa_list', 'list.txt'] + OUT + opt + [mname], dict(mats, m1=outs['p.po'], m2=outs['p.pir'], m3=outs['c.aln'], **{'list.txt': 'm1\nm2\nm3\n'}), mname)
    if i % 6 == 3:
        o3 = [g for g in gen if g[0] == mname][(i // 6) % max(1, len([g for g in gen if g[0] == mname]))]
        add('subset2', ['-read_msa', 'm1', '-read_msa2', 'm2', '-subset', 's1.txt', '-subset2', 's2.txt'] + OUT + [mname],
            dict(mats, m1=outs['p.po'], m2=o3[2]['p.po'], **{'s1.txt': 'SOURCENAME=s1\nSOURCENAME=s2\n', 's2.txt': 'SOURCENAME=s1\n'}), mname)
        add('remove2', ['-read_msa', 'm1', '-read_msa2', 'm2', '-remove', 's1.txt', '-remove2', 's2.txt'] + OUT + [mname],
            dict(mats, m1=outs['p.pir'], m2=o3[2]['c.aln'], **{'s1.txt': 'SOURCENAME=s1\n', 's2.txt': 'SOURCENAME=s1\r\nSOURCENAME=s2\r\n'}), mname)

# 2. hand-made files
P = lambda s: s
clu = 'CLUSTAL W (1.83) multiple sequence alignment\n\n\nsa           MKVLAAGIVG-LA\nsb           MKVL-AGIVGALA\nsc           MKVLAAG--GALA\n             ****  **  ***\n\nsa           KKK\nsb           KK-\nsc           K--\n\n'
hm = [
    ('clustal', clu),
    ('clustal-nohdr', clu.split('\n', 3)[3]),
    ('clustal-comment', '# a comment\n' + clu),
    ('clustal-digits', clu.replace('MKVLAAGIVG-LA', 'MKVL AAGI VG-LA 13')),
    ('clustal-crlf', clu.replace('\n', '\r\n')),
    ('clustal-badname', clu.replace('sb           KK-', 'zz           KK-')),
    ('clustal-onelineperseq', 'CLUSTAL\nsa MKVLAAGIVG\nsb MKVLAGIVG\n'),
    ('clustal-short', 'CLUSTAL\nsa\n'),
    ('clustal-long', 'CLUSTAL\n\nsa ' + 'MKVLAAGIVGLA' * 50 + '\nsb ' + 'MKVLAGIVGALA' * 50 + '\n'),
    ('clustal-long2', 'CLUSTAL\n\nsa ' + 'MKVLAAGIVGLA' * 42 + '\nsb ' + 'MKVLAGIVGALA' * 42 + '\n'),
    ('clustal-allgap', 'CLUSTAL\n\nsa MKVLAAG\nsb -------\nsc MKVLAAG\n'),
    ('clustal-unknown-header', 'sa MKVLAAG\nsb MKVLA-G\n'),
    ('clustal-lower', 'CLUSTAL\n\nsa mkvlaag\nsb mkvla-g\nsc MKVLAAG\n'),
    ('clustal-3seq-unequal', 'CLUSTAL\n\nsa MKVLAAGIVG\nsb MKVLA\nsc MK-LAAGIVGAA\n'),
    ('pir', '>sa first one\nMKVLAAGIVG-LA\n>sb\nMKVL-AGIVGALA\n>sc some title\nMKVLAAG--GALA\n'),
    ('pir-wrapped', '>sa\nMKVLAA\nGIVG-L\nA\n>sb\nMKVL-A\nGIVGAL\nA\n'),
    ('pir-blank', '>sa\nMKVLAA\n\nGIVG\n>sb\nMKVL-A\nGIVGAL\nA\n'),
    ('pir-comment', '# c\n>sa\n#x\nMKVLAAGIVG\n*y\n>sb\nMKVLAGIVG\n'),
    ('pir-digits', '>sa\nMKVL AAGI 10\nVG\n>sb\nMKVLAGIVG\n'),
    ('pir-crlf', '>sa t\r\nMKVLAAGIVG\r\n>sb\r\nMKVLAGIVG\r\n'),
    ('pir-long', '>sa\n' + 'MKVLAAGIVGLA' * 50 + '\n>sb\n' + 'MKVLAGIVGALA' * 50 + '\n'),
    ('pir-one', '>sa\nMKVLAAGIVG\n'),
    ('pir-noname', '>\nMKVLAAG\n>sb\nMKVLAAG\n'),
    ('pir-emptyseq', '>sa\nMKVLAAG\n>sb\n>sc\nMKVLAAG\n'),
    ('pir-lower', '>sa\nmkvlaag\n>sb\nMKVLAAG\n'),
    ('pir-twice', '>sa\nMKVLAAG\n>sa\nMKVLAAG\n>sb\nMKV\n'),
    ('pir-brackets', '>sa\nMK[VL]AAG\n>sb\nMK?LAAG\n'),
    ('pir-longname', '>' + 'n' * 33 + '\nMKVLAAG\n>sb\nMKVLAAG\n'),
    ('pir-name35', '>' + 'n' * 35 + '\nMKVLAAG\n>sb\nMKVLAAG\n'),
    ('empty', ''),
    ('blank', '\n\n'),
    ('comment-only', '# nothing\n* nothing\n'),
    ('garbage', 'hello world\nthis is not an alignment\n'),
    ('po-minimal', 'VERSION=LPO.0.1\nNAME=x\nTITLE=t\nLENGTH=3\nSOURCECOUNT=1\nSOURCENAME=x\nSOURCEINFO=3 0 1 -1 t\nM:S0\nK:L0S0\nV:L1S0\n'),
    ('po-notitle', 'VERSION=LPO.0.1\nNAME=x\nTITLE=\nLENGTH=3\nSOURCECOUNT=1\nSOURCENAME=x\nSOURCEINFO=3 0 1 -1\nM:S0\nK:L0S0\nV:L1S0\n'),
    ('po-two', 'VERSION=LPO.0.1\nNAME=x\nTITLE=t\nLENGTH=4\nSOURCECOUNT=2\nSOURCENAME=x\nSOURCEINFO=3 0 1 -1 t\nSOURCENAME=y\nSOURCEINFO=3 0 1 -1 u v\nM:S0S1\nK:L0S0S1\nA:L1S0\nV:L1S1A2\n'),
    ('po-crlf', 'VERSION=LPO.0.1\r\nNAME=x\r\nTITLE=t\r\nLENGTH=3\r\nSOURCECOUNT=1\r\nSOURCENAME=x\r\nSOURCEINFO=3 0 1 -1 t\r\nM:S0\r\nK:L0S0\r\nV:L1S0\r\n'),
    ('po-truncated', 'VERSION=LPO.0.1\nNAME=x\nTITLE=t\nLENGTH=3\nSOURCECOUNT=1\nSOURCENAME=x\nSOURCEINFO=3 0 1 -1 t\nM:S0\n'),
    ('po-noheader', 'NAME=x\nTITLE=t\nLENGTH=3\nSOURCECOUNT=1\nSOURCENAME=x\nSOURCEINFO=3 0 1 -1 t\nM:S0\nK:L0S0\nV:L1S0\n'),
]
for kind, text in hm:
    for opt in ([], ['-do_global'], ['-hb'], ['-toupper'], ['-tolower']):
        add('hand-' + kind, ['-read_msa', 'm1'] + OUT + opt + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': text}, 'blosum80.mat')
    add('hand2-' + kind, ['-read_msa', 'm1', '-read_msa2', 'm2'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': text, 'm2': clu}, 'blosum80.mat')
    add('handsub-' + kind, ['-read_msa', 'm1', '-subset', 'sub.txt'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': text, 'sub.txt': 'SOURCENAME=sa\nSOURCENAME=sc\n'}, 'blosum80.mat')
    add('handrm-' + kind, ['-read_msa', 'm1', '-remove', 'sub.txt'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': text, 'sub.txt': 'SOURCENAME=sa\r\n'}, 'blosum80.mat')
# errors of the command line
add('err', ['-read_msa', 'nosuch'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-subset', 'nosuch'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu}, 'blosum80.mat')
add('err', ['-read_msa_list', 'nosuch'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80}, 'blosum80.mat')
add('err', ['-read_msa_list', 'l', '-read_msa', 'm1'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'l': 'm1\n', 'm1': clu}, 'blosum80.mat')
add('err', ['-read_msa_list', 'l'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'l': ''}, 'blosum80.mat')
add('err', ['-read_msa_list', 'l'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'l': 'm1 nosuch\n', 'm1': clu}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-subset', 's', '-remove', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu, 's': 'SOURCENAME=sa\n'}, 'blosum80.mat')
add('err', ['-subset', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 's': 'SOURCENAME=sa\n'}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-subset', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu, 's': 'SOURCENAME=zz\n'}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-subset', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu, 's': ''}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-remove', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu, 's': 'SOURCENAME=sa\nSOURCENAME=sb\nSOURCENAME=sc\n'}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-subset', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu, 's': 'SOURCENAME=sa'}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-subset', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu, 's': 'xx\nSOURCENAME=sa\n'}, 'blosum80.mat')
add('err', ['-read_msa', 'm1', '-subset', 's'] + OUT + ['blosum80.mat'], {'blosum80.mat': B80, 'm1': clu, 's': 'SOURCENAME=sa\n\n\nSOURCENAME=sb\n'}, 'blosum80.mat')

rows = []
for kind, args, files, matrix in cases:
    p, outs = run_poa(args, files, work)
    err = ''.join(l + '\n' for l in p.stderr.split('\n') if l.startswith('...'))
    rows.append((kind, args, files, p.returncode, p.stdout, err, outs))

ver = subprocess.run(['dpkg-query', '-W', binary], capture_output=True, text=True).stdout.strip().replace('\t', ' ')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='\n', encoding='utf-8') as fh:
    fh.write('# ' + ver + '; uname -m ' + os.uname().machine + '\n')
    fh.write('id\tkind\targs\tfiles\trc\tstdout\tstderr\toutfiles\n')
    for i, (kind, args, files, rc, so, err, outs) in enumerate(rows):
        fh.write('\t'.join(['M%04d' % (i + 1), kind, esc(SEP.join(args)), esc(SEP.join(k + '=' + v for k, v in files.items())), str(rc), esc(so), esc(err),
                            esc(SEP.join(k + '=' + v for k, v in outs.items()))]) + '\n')
print(len(rows), 'rows;', sum(1 for r in rows if r[3] != 0), 'with a non-zero exit;', sorted(set(r[3] for r in rows)), ver)
