#!/usr/bin/env python3
"""B09 WP61 (F280) oracle of the sequence manipulation programs of EMBOSS 6.6.0 that read sequences and write sequences:
revseq transeq extractseq cutseq degapseq trimseq skipseq nthseq notseq pasteseq (and the programs added to PROGS below).
    <program> -auto -sequence in.dat <program options> [-sbegin1 N -send1 N -sreverse1 -slower1 -supper1 -sformat1 F] -outseq stdout [-osformat F]
Inputs: one to four sequences (nucleic or protein; gaps, ambiguity codes, lower / mixed case, U, empty), damaged files, formats fasta embl genbank gcg nbrf raw ig ncbi staden codata fastq.
JSONL.GZ rows: id, prog, args, files, rc, out, err, outfiles, hang (see emboss_prog_oracle_lib.py).
Usage: python3 scripts/alignment/generate_emboss_prog_seq_oracle.py PROGRAM OUT.jsonl.gz [n=400] [seed]      (needs EMBOSS 6.6.0)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402
import generate_seqret_oracle as S  # noqa: E402

IN_FORMATS = ['fasta'] * 8 + ['embl', 'genbank', 'gcg', 'nbrf', 'raw', 'ig', 'ncbi', 'staden', 'codata', 'fastq']
OUT_FORMATS = ['fasta'] * 6 + ['embl', 'genbank', 'gcg', 'nbrf', 'raw', 'msf', 'clustal', 'phylip', 'ig', 'swiss', 'fastq', 'ncbi', 'text', 'staden', 'pir', 'codata', 'strider', 'mega', 'nexus', 'gff3', 'sam', 'fitch', 'selex', 'hennig86']
LAST_DAMAGED = False
TABLES = ['0', '1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23']


def gapped(rnd, s):
    """gap characters in a sequence: inner gaps, a leading / trailing run"""
    s = list(s)
    for _ in range(rnd.choice([0, 0, 1, 3, 8])):
        if s: s.insert(rnd.randrange(len(s) + 1), rnd.choice('-.~'))
    if rnd.random() < 0.3: s = [rnd.choice('-.~') for _ in range(rnd.randint(1, 5))] + s
    if rnd.random() < 0.3: s = s + [rnd.choice('-.~') for _ in range(rnd.randint(1, 5))]
    if not any(c.isalpha() for c in s): s.append('A')           # a sequence of gaps only makes the writers of some set formats overwrite memory
    return ''.join(s)


def ends(rnd, s):
    """ambiguity / unknown runs at the ends (trimseq)"""
    pad = lambda: ''.join(rnd.choice('NnXxRYKM-.') for _ in range(rnd.randint(0, 12)))
    return pad() + s + pad() if rnd.random() < 0.7 else s


def seqs_for(rnd, prog, kind='n', nseq=None, lens=(0, 200), gaps=False, endpad=False):
    opts = dict(orfs=True, iupac=rnd.choice([0, 0, 0, 0.02, 0.1]), case=rnd.choice(['upper', 'upper', 'upper', 'lower', 'mixed']), rna=rnd.random() < 0.08, x=rnd.random() < 0.1)
    items = L.rand_seqs(rnd, kind, nseq=nseq, lens=lens, **opts) if kind == 'n' else L.rand_seqs(rnd, 'p', nseq=nseq, lens=lens)
    out = []
    for name, desc, s in items:
        if gaps: s = gapped(rnd, s)
        if endpad: s = ends(rnd, s)
        out.append((name, desc, s))
    return out


def seq_file(rnd, d, items, fmt=None):
    text = L.fasta(rnd, items)
    fmt = fmt or rnd.choice(IN_FORMATS)
    text = L.convert(text, fmt, d)
    global LAST_DAMAGED
    LAST_DAMAGED = rnd.random() < 0.05
    if LAST_DAMAGED: text = S.damaged(rnd, text)
    return text, fmt


def seq_quals(rnd, fmt, suffix='1', p=0.15):
    a = []
    if rnd.random() < p: a += ['-sbegin' + suffix, str(rnd.randint(-15, 40))]
    if rnd.random() < p: a += ['-send' + suffix, str(rnd.randint(-40, 120))]
    if rnd.random() < p / 2: a += ['-sreverse' + suffix]
    if rnd.random() < p / 3: a += [rnd.choice(['-slower', '-supper']) + suffix]
    if rnd.random() < p and fmt != 'fasta': a += ['-sformat' + suffix, fmt]
    return a


def out_quals(rnd, p=0.25):
    f = 'fasta' if rnd.random() > p else rnd.choice(OUT_FORMATS)
    return ['-osformat', f] if f != 'fasta' else []


def common(rnd, d, prog, kind='n', nseq=None, lens=(0, 200), gaps=False, endpad=False, rev_ok=True):
    items = seqs_for(rnd, prog, kind, nseq, lens, gaps, endpad)
    text, fmt = seq_file(rnd, d, items)
    return text, fmt, ['-sequence', 'in.dat'] + seq_quals(rnd, fmt)


def gen_revseq(rnd, i, d):
    text, fmt, args = common(rnd, d, 'revseq', 'n', gaps=True)
    if rnd.random() < 0.05:
        text, fmt = seq_file(rnd, d, seqs_for(rnd, 'revseq', 'p', 1)); args = ['-sequence', 'in.dat']
    for o in ['-noreverse', '-nocomplement', '-notag']:
        if rnd.random() < 0.25: args.append(o)
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_transeq(rnd, i, d):
    text, fmt, args = common(rnd, d, 'transeq', 'n', lens=(0, 300))
    if rnd.random() < 0.7: args += ['-frame', rnd.choice(['1', '2', '3', 'F', '-1', '-2', '-3', 'R', '6', '1,2', '1,-1', '2,3,-2', 'f', 'r', '1,F', '-1,-3'])]
    if rnd.random() < 0.5: args += ['-table', rnd.choice(TABLES)]
    if rnd.random() < 0.3:
        n = rnd.randint(1, 3)
        args += ['-regions', ','.join('%d-%d' % (a, a + rnd.randint(0, 60)) for a in [rnd.randint(1, 150) for _ in range(n)])]
    for o in ['-trim', '-clean', '-alternative']:
        if rnd.random() < 0.25: args.append(o)
    if rnd.random() < 0.25: args.append('-nomethionine')
    if rnd.random() < 0.03: args += rnd.choice([['-table', '7'], ['-frame', '7'], ['-regions', '5-3'], ['-regions', '0-4'], ['-regions', '3']])
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_extractseq(rnd, i, d):
    kind = 'n' if rnd.random() < 0.7 else 'p'
    text, fmt = seq_file(rnd, d, seqs_for(rnd, 'extractseq', kind, 1, lens=(1, 200)))
    args = ['-sequence', 'in.dat'] + seq_quals(rnd, fmt)
    if rnd.random() < 0.8:
        n = rnd.choice([1, 1, 2, 3])
        args += ['-regions', rnd.choice([',', ' ', ';']).join('%d%s%d' % (a, rnd.choice(['-', ':', '..', '=']), a + rnd.randint(0, 50)) for a in sorted(rnd.randint(1, 180) for _ in range(n)))]
    if rnd.random() < 0.4: args.append('-separate')
    if rnd.random() < 0.03: args += rnd.choice([['-regions', '0-3'], ['-regions', '5-3'], ['-regions', '3']])
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_cutseq(rnd, i, d):
    kind = 'n' if rnd.random() < 0.6 else 'p'
    text, fmt = seq_file(rnd, d, seqs_for(rnd, 'cutseq', kind, 1, lens=(1, 120), gaps=rnd.random() < 0.2))
    args = ['-sequence', 'in.dat'] + seq_quals(rnd, fmt, p=0.2)
    if rnd.random() < 0.85: args += ['-from', str(rnd.randint(-3, 130))]
    if rnd.random() < 0.85: args += ['-to', str(rnd.randint(-3, 130))]
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_degapseq(rnd, i, d):
    kind = 'n' if rnd.random() < 0.6 else 'p'
    text, fmt, args = common(rnd, d, 'degapseq', kind, gaps=True)
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_trimseq(rnd, i, d):
    kind = 'n' if rnd.random() < 0.6 else 'p'
    items = seqs_for(rnd, 'trimseq', kind, None, (5, 150), rnd.random() < 0.4, True)
    items = [(n, dsc, x if len(x) >= 3 else x + 'ACG') for n, dsc, x in items]      # a window that is not shorter than the sequence makes the program read outside it
    text, fmt = seq_file(rnd, d, items)
    args = ['-sequence', 'in.dat'] + seq_quals(rnd, fmt)
    minlen = min(len(x[2]) for x in items)
    if rnd.random() < 0.5 and minlen >= 2 and fmt == 'fasta' and not LAST_DAMAGED and not any(c in args for c in ('-sbegin1', '-send1')):
        args += ['-window', str(rnd.choice([w for w in [1, 2, 3, 5, 8, 10, 20] if w < minlen]))]
    if rnd.random() < 0.4: args += ['-percent', str(rnd.choice([0, 10, 25, 33.3, 50, 75, 100, 101]))]
    for o in ['-strict', '-star', '-noleft', '-noright']:
        if rnd.random() < 0.25: args.append(o)
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_skipseq(rnd, i, d):
    kind = rnd.choice(['n', 'p'])
    text, fmt, args = common(rnd, d, 'skipseq', kind, nseq=rnd.randint(1, 6), lens=(1, 80))
    if rnd.random() < 0.85: args += ['-skip', str(rnd.randint(-1, 7))]
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_nthseq(rnd, i, d):
    kind = rnd.choice(['n', 'p'])
    text, fmt, args = common(rnd, d, 'nthseq', kind, nseq=rnd.randint(1, 6), lens=(1, 80))
    if rnd.random() < 0.9: args += ['-number', str(rnd.choice([1, 2, 3, 4, 5, 6, 7, 9]))]
    elif rnd.random() < 0.5: args += ['-number', '0']
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_notseq(rnd, i, d):
    kind = rnd.choice(['n', 'p'])
    items = seqs_for(rnd, 'notseq', kind, rnd.randint(1, 6), lens=(1, 60))
    text, fmt = seq_file(rnd, d, items)
    args = ['-sequence', 'in.dat'] + seq_quals(rnd, fmt, p=0.05)
    names = [n for n, _, _ in items]
    files = {'in.dat': text}
    pats = [names[rnd.randrange(len(names))], names[rnd.randrange(len(names))][:3] + '*', '*' + names[0][-2:], 'zzz', names[0].lower(), ','.join(names[:2]), ';'.join(names[-2:]), names[0][:2] + '?' + names[0][3:], '*']
    ex = rnd.choice(pats)
    if rnd.random() < 0.25:
        files['names.txt'] = '\n'.join(['# comment', '! other'] + [rnd.choice(pats) for _ in range(rnd.randint(1, 3))] + ['']) 
        ex = '@names.txt'
    args += ['-exclude', ex]
    out = 'stdout'
    if rnd.random() < 0.3:
        j = rnd.choice(['junk.fa', 'stdout'])
        args += ['-junkoutseq', j]
        if j == 'stdout': out = 'main.fa'
    return args + out_quals(rnd, 0.1) + ['-outseq', out], files


def gen_pasteseq(rnd, i, d):
    kind = 'n' if rnd.random() < 0.7 else 'p'
    ta, fa = seq_file(rnd, d, seqs_for(rnd, 'pasteseq', kind, 1, lens=(1, 80)))
    tb, fb = seq_file(rnd, d, seqs_for(rnd, 'pasteseq', kind if rnd.random() < 0.9 else rnd.choice(['n', 'p']), 1, lens=(1, 50)))
    args = ['-asequence', 'a.dat', '-bsequence', 'b.dat'] + seq_quals(rnd, fa, '1', 0.1) + seq_quals(rnd, fb, '2', 0.1)
    if rnd.random() < 0.85: args += ['-pos', str(rnd.randint(-2, 90))]
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'a.dat': ta, 'b.dat': tb}


def tailed(rnd, s):
    """a sequence with a poly-A tail and / or a poly-T head (trimest)"""
    s = s.replace('A', 'C') if rnd.random() < 0.3 else s
    if rnd.random() < 0.6: s = s + ''.join(rnd.choice('AAAAAAAAN' if rnd.random() < 0.5 else 'AAAAAAAAAAAC') for _ in range(rnd.randint(2, 25)))
    if rnd.random() < 0.4: s = ''.join(rnd.choice('TTTTTTTTN' if rnd.random() < 0.5 else 'TTTTTTTTTTTG') for _ in range(rnd.randint(2, 25))) + s
    return s


def gen_trimest(rnd, i, d):
    items = [(n, dsc, tailed(rnd, x)) for n, dsc, x in seqs_for(rnd, 'trimest', 'n', None, (5, 120))]
    text, fmt = seq_file(rnd, d, items)
    args = ['-sequence', 'in.dat'] + seq_quals(rnd, fmt, p=0.05)
    if rnd.random() < 0.5: args += ['-minlength', str(rnd.choice([1, 2, 3, 4, 6, 10]))]
    if rnd.random() < 0.4: args += ['-mismatches', str(rnd.choice([0, 1, 2, 3]))]
    for o in ['-noreverse', '-tolower', '-nofiveprime']:
        if rnd.random() < 0.25: args.append(o)
    if rnd.random() < 0.03: args += ['-minlength', '0']
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_newseq(rnd, i, d):
    typ = rnd.choice(['N', 'N', 'P'])
    alpha = 'ACGTRYKMSWNBDHVXacgtuUn-.~*' if typ == 'N' else 'ACDEFGHIKLMNPQRSTVWYBZXUOJacdef-.*?'
    seq = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(0, 60)))
    args = ['-name', rnd.choice(['nn', 'My seq', 'x1', 'a_b', 'seq 2', '']), '-description', rnd.choice(['d d', 'a description of it', '', 'x']), '-type', typ, '-sequence', seq]
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {}


def gen_yank(rnd, i, d):
    kind = rnd.choice(['n', 'p'])
    items = seqs_for(rnd, 'yank', kind, rnd.randint(1, 3), (3, 80))
    text, fmt = seq_file(rnd, d, items)
    name = items[rnd.randrange(len(items))][0]
    spec = rnd.choice(['in.dat', 'in.dat', 'fasta::in.dat', 'in.dat:' + name, 'in.dat[%d:%d]' % (rnd.randint(1, 5), rnd.randint(6, 30)), 'in.dat[%d:%d:r]' % (rnd.randint(1, 5), rnd.randint(6, 30)), 'asis::ACGTACGTAC', 'in.dat:' + name.lower(), 'in.dat:{%s|%s}' % (name, name)])
    args = ['-sequence', spec] + seq_quals(rnd, fmt, p=0.15)
    files = {'in.dat': text}
    if rnd.random() < 0.3: files['out.txt'] = 'earlier line\n'
    if rnd.random() < 0.2: args.append('-newfile')
    return args + ['-outfile', 'out.txt'], files


def gen_seqcount(rnd, i, d):
    kind = rnd.choice(['n', 'p'])
    text, fmt, args = common(rnd, d, 'seqcount', kind, nseq=rnd.randint(1, 7), lens=(1, 60))
    return args + ['-outfile', 'out.txt'], {'in.dat': text}


def gen_splitter(rnd, i, d):
    kind = rnd.choice(['n', 'p'])
    text, fmt, args = common(rnd, d, 'splitter', kind, nseq=rnd.randint(1, 3), lens=(1, 120))
    if rnd.random() < 0.9: args += ['-size', str(rnd.choice([1, 2, 5, 8, 10, 15, 30, 50, 100]))]
    if rnd.random() < 0.5: args += ['-overlap', str(rnd.choice([0, 1, 2, 3, 5, 9, 20]))]
    if rnd.random() < 0.4: args.append('-addoverlap')
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_maskseq(rnd, i, d):
    kind = rnd.choice(['n', 'p'])
    text, fmt = seq_file(rnd, d, seqs_for(rnd, 'maskseq', kind, 1, lens=(1, 100)))
    args = ['-sequence', 'in.dat'] + seq_quals(rnd, fmt, p=0.1)
    if '-sbegin1' in args and int(args[args.index('-sbegin1') + 1]) < 0:
        i = args.index('-sbegin1'); del args[i:i + 2]
    sb = int(args[args.index('-sbegin1') + 1]) if '-sbegin1' in args else 1
    if rnd.random() < 0.9:
        n = rnd.choice([1, 1, 2, 3])
        lo = max(1, sb) if sb > 0 else 1
        args += ['-regions', ','.join('%d-%d' % (a, a + rnd.randint(0, 25)) for a in sorted(rnd.randint(lo + 1, lo + 60) for _ in range(n)))]
    if rnd.random() < 0.3: args.append('-tolower')
    if rnd.random() < 0.4: args += ['-maskchar', rnd.choice(['N', 'X', '#', '?', ' ', '-', 'z', '*'])]
    if rnd.random() < 0.03: args += ['-maskchar', rnd.choice(['', 'NN'])]
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_maskambignuc(rnd, i, d):
    text, fmt, args = common(rnd, d, 'maskambignuc', 'n', gaps=True)
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def gen_maskambigprot(rnd, i, d):
    text, fmt, args = common(rnd, d, 'maskambigprot', 'p', gaps=True)
    return args + out_quals(rnd) + ['-outseq', 'stdout'], {'in.dat': text}


def text_lines(rnd, cr=False):
    """lines of text with blanks, tabs, runs of spaces, empty lines, trailing white space, form feeds and (cr) carriage returns"""
    words = ['alpha', 'beta', 'x', 'gene 1', 'a.b', '', ' ', '\t', '  ', 'ACGT', '12', '#c', '\x0b']
    lines = []
    for _ in range(rnd.randint(0, 8)):
        l = ''.join(rnd.choice(words) + rnd.choice(['', ' ', '  ', '\t', ' \t ', '   ']) for _ in range(rnd.randint(0, 6)))
        if cr and rnd.random() < 0.4: l = l.replace(' ', '\r', 1) if ' ' in l else l + '\r'
        lines.append(l)
    end = rnd.choice(['\n', '\n', '\n', '', '\r\n', '\n\n'])
    sep = rnd.choice(['\n', '\n', '\r\n', '\r'])
    return sep.join(lines) + end


def gen_trimspace(rnd, i, d):
    return ['-infile', 'in.txt', '-outfile', rnd.choice(['out.txt', 'stdout'])] if rnd.random() < 0.9 else ['-infile', 'in.txt'], {'in.txt': text_lines(rnd)}


def gen_nospace(rnd, i, d):
    args = ['-infile', 'in.txt']
    if rnd.random() < 0.8: args += ['-menu', rnd.choice(['all', 'end', 'excess'])]
    if rnd.random() < 0.04: args += ['-menu', rnd.choice(['x', 'e', 'a'])]
    return args + ['-outfile', rnd.choice(['out.txt', 'stdout'])], {'in.txt': text_lines(rnd)}


def gen_notab(rnd, i, d):
    return ['-infile', 'in.txt', '-outfile', rnd.choice(['out.txt', 'stdout'])], {'in.txt': text_lines(rnd)}


def gen_noreturn(rnd, i, d):
    args = ['-infile', 'in.txt']
    if rnd.random() < 0.8: args += ['-system', rnd.choice(['unix', 'pc', 'mac'])]
    return args + ['-outfile', rnd.choice(['out.txt', 'stdout'])], {'in.txt': text_lines(rnd, True)}


GENS = {'revseq': gen_revseq, 'transeq': gen_transeq, 'extractseq': gen_extractseq, 'cutseq': gen_cutseq, 'degapseq': gen_degapseq, 'trimseq': gen_trimseq,
        'skipseq': gen_skipseq, 'nthseq': gen_nthseq, 'notseq': gen_notseq, 'pasteseq': gen_pasteseq, 'trimest': gen_trimest, 'newseq': gen_newseq, 'yank': gen_yank,
        'seqcount': gen_seqcount, 'trimspace': gen_trimspace, 'nospace': gen_nospace, 'notab': gen_notab, 'noreturn': gen_noreturn, 'splitter': gen_splitter, 'maskseq': gen_maskseq, 'maskambignuc': gen_maskambignuc, 'maskambigprot': gen_maskambigprot}

if __name__ == '__main__':
    prog = sys.argv.pop(1)
    L.main(prog, GENS[prog], 400)
