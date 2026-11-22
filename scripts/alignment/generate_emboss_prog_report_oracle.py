#!/usr/bin/env python3
"""B09 WP62 (F281) oracle of the EMBOSS 6.6.0 programs that read sets of sequences or write reports about sequences: listor nthseqset (and the programs added to GENS below).
    <program> -auto ARGS ...      every row: arguments, input files, standard output, error stream, exit status, written files (see emboss_prog_oracle_lib.py)
Usage: python3 scripts/alignment/generate_emboss_prog_report_oracle.py PROGRAM OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402
import generate_seqret_oracle as S  # noqa: E402
import generate_emboss_prog_seq_oracle as G  # noqa: E402

SET_FORMATS = ['fasta'] * 4 + ['clustal', 'msf', 'clustal', 'msf', 'selex', 'stockholm', 'phylip', 'nexus', 'mega', 'selex', 'stockholm', 'gcg', 'embl', 'genbank', 'pir', 'ncbi', 'raw']


def related_seqs(rnd, kind, base=None):
    """two lists of (name, desc, residues) with common members (equal ignoring the case) and equal lengths"""
    base = base or G.seqs_for(rnd, 'x', kind, rnd.randint(2, 6), (1, 40))
    other = []
    for name, desc, s in base:
        r = rnd.random()
        if r < 0.35: other.append((name + 'b', desc, s.lower() if rnd.random() < 0.5 else s))
        elif r < 0.55: other.append((name + 'c', desc, ''.join(rnd.choice('ACGT' if kind == 'n' else 'ACDEFG') for _ in range(len(s)))))
    other += G.seqs_for(rnd, 'x', kind, rnd.randint(0, 3), (1, 40))
    if not other: other = G.seqs_for(rnd, 'x', kind, 1, (1, 40))
    rnd.shuffle(other)
    return base, other


def gen_listor(rnd, i, d):
    kind = rnd.choice(['n', 'n', 'p'])
    a, b = related_seqs(rnd, kind)
    ta, fa = G.seq_file(rnd, d, a, rnd.choice(['fasta'] * 4 + ['embl', 'genbank', 'pir', 'gcg', 'raw']))
    tb, fb = G.seq_file(rnd, d, b, rnd.choice(['fasta'] * 4 + ['embl', 'genbank', 'pir', 'gcg', 'raw']))
    args = ['-firstsequences', 'a.dat', '-secondsequences', 'b.dat']
    if rnd.random() < 0.85: args += ['-operator', rnd.choice(['OR', 'AND', 'XOR', 'NOT', 'or', 'and', 'xor', 'not', 'A', 'X', 'N', 'O', 'AN', 'XO', 'No', 'merger', 'only those in both sets', 'XOR - only those not in both sets'])]
    if rnd.random() < 0.04: args += ['-operator', rnd.choice(['foo', '', 'NOR', 'OR,AND', '*'])]
    args += G.seq_quals(rnd, fa, '1', 0.2) + G.seq_quals(rnd, fb, '2', 0.2)
    if rnd.random() < 0.9: args += ['-outfile', rnd.choice(['out.list', 'stdout'])]   # (no -outfile: the default file name BASE.program of the name of the first sequence)
    return args, {'a.dat': ta, 'b.dat': tb}


def multi_set_text(rnd, d, kind, fmt, nsets):
    """a file of several datasets: the same writer called once per set, joined"""
    parts = []
    for _ in range(nsets):
        items = G.seqs_for(rnd, 'x', kind, rnd.randint(2, 4), (4, 30))
        n = max(len(s) for _, _, s in items)
        items = [(nme, desc, (s + 'ACGT' * 10)[:n] if rnd.random() < 0.5 else s) for nme, desc, s in items]
        parts.append(L.convert(L.fasta(rnd, items), fmt, d))
    return parts


def gen_nthseqset(rnd, i, d):
    kind = rnd.choice(['n', 'n', 'p'])
    fmt = rnd.choice(SET_FORMATS)
    nsets = rnd.choice([1, 2, 2, 3]) if fmt in ('clustal', 'msf', 'selex', 'stockholm') else rnd.choice([1, 1, 1, 2, 3]) if fmt == 'phylip' else 1       # the sets of clustal / msf / selex / stockholm are one set for the reader (no Multidone): the oracle is run on several datasets in one file
    if nsets > 1:
        text = ''.join(multi_set_text(rnd, d, kind, fmt, nsets))
    else:
        text, fmt = G.seq_file(rnd, d, G.seqs_for(rnd, 'x', kind, rnd.randint(1, 5), (1, 80)), fmt)
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.12)
    if rnd.random() < 0.85: args += ['-number', str(rnd.choice([1, 1, 2, 3, 4]))]
    elif rnd.random() < 0.5: args += ['-number', '0']
    return args + G.out_quals(rnd, 0.4) + ['-outseq', 'stdout'], {'in.dat': text}, ({'vg': True} if fmt in ('nexus', 'selex', 'stockholm', 'mega', 'msf', 'clustal') else {})


def wrap(s, n=60):
    return [s[i:i + n] for i in range(0, len(s), n)] or ['']


def entry(rnd, kind, fmt, name, acc, desc, org, seq):
    """a database entry in genbank / embl / swiss format with accession, version, GI, organism"""
    if fmt == 'genbank':
        out = 'LOCUS       %-10s %9d %s    %s  linear   UNC 01-JAN-2000\n' % (name, len(seq), 'bp' if kind == 'n' else 'aa', 'DNA' if kind == 'n' else '')
        out += 'DEFINITION  %s.\nACCESSION   %s\nVERSION     %s.%d  GI:%d\n' % (desc, acc, acc, rnd.randint(1, 4), rnd.randint(1000, 99999))
        out += 'SOURCE      %s\n  ORGANISM  %s\n            Eukaryota; Metazoa.\n' % (org, org)
        out += 'ORIGIN\n' + ''.join('%9d %s\n' % (i * 60 + 1, ' '.join(seq[i * 60 + j:i * 60 + j + 10].lower() for j in range(0, 60, 10) if seq[i * 60 + j:i * 60 + j + 10])) for i in range((len(seq) + 59) // 60)) + '//\n'
        return out
    if fmt == 'embl':
        out = 'ID   %s; SV %d; linear; genomic DNA; STD; HUM; %d BP.\nXX\nAC   %s;\nXX\nDE   %s\nXX\nOS   %s\nXX\n' % (name, rnd.randint(1, 3), len(seq), acc, desc, org)
        out += '     %s %9d\n' % (' '.join(seq[j:j + 10].lower() for j in range(0, min(len(seq), 60), 10)), len(seq)) + '//\n'
        return out
    out = 'ID   %s_HUMAN    Reviewed;    %d AA.\nAC   %s;\nDE   RecName: Full=%s;\nOS   %s.\nSQ   SEQUENCE   %d AA;  1000 MW;  ABCD CRC64;\n' % (name, len(seq), acc, desc, org, len(seq))
    out += '     ' + ' '.join(seq[j:j + 10] for j in range(0, min(len(seq), 60), 10)) + '\n//\n'
    return out


def gen_infoseq(rnd, i, d):
    kind = rnd.choice(['n', 'n', 'p'])
    fmt = rnd.choice(['fasta', 'fasta', 'genbank', 'embl', 'swiss', 'ncbi', 'gcg', 'pir', 'raw', 'fasta'])
    if fmt == 'swiss' and kind == 'n': fmt = 'embl'
    items = G.seqs_for(rnd, 'x', kind, rnd.randint(1, 4), (1, 90))
    if fmt in ('genbank', 'embl', 'swiss'):
        text = ''.join(entry(rnd, kind, fmt, nme.upper()[:10], 'X%05d' % rnd.randint(1, 99999), desc or 'some description', rnd.choice(['Homo sapiens', 'Mus musculus', 'Escherichia coli', '']), s.upper()) for nme, desc, s in items)
    elif fmt == 'ncbi':
        text = ''.join('>%s|%s|%s %s\n%s\n' % (rnd.choice(['gb', 'sp', 'emb', 'gi|%d|gb' % rnd.randint(1, 9999), 'ref', 'pdb', 'gnl|db', 'lcl', 'pir']), 'AB%06d' % rnd.randint(1, 999999), nme, desc, s) for nme, desc, s in items)
        fmt = 'fasta'
    else:
        text, fmt = G.seq_file(rnd, d, items, fmt)
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '1', 0.12)
    for o in ['-nohead', '-nocolumns', '-html', '-only', '-gi', '-seqversion']:
        if rnd.random() < 0.15: args.append(o)
    if '-only' in args:
        for o in ['-usa', '-database', '-name', '-accession', '-type', '-length', '-pgc', '-organism', '-description']:
            if rnd.random() < 0.3: args.append(o)
    else:
        for o in ['-nousa', '-nodatabase', '-noname', '-noaccession', '-notype', '-nolength', '-nopgc', '-noorganism', '-nodescription']:
            if rnd.random() < 0.08: args.append(o)
    if rnd.random() < 0.3: args += ['-delimiter', rnd.choice(['|', ',', ' ', ';;', '\\t', ''])]
    args += ['-outfile', rnd.choice(['out.txt', 'stdout'])]
    return args, {'in.dat': text}


def gen_wordcount(rnd, i, d):
    kind = rnd.choice(['n', 'n', 'p'])
    lens = rnd.choice([(1, 40), (5, 200), (50, 600)])
    items = G.seqs_for(rnd, 'x', kind, rnd.choice([1, 1, 2, 3, 5]), lens)
    if rnd.random() < 0.3:   # a short first sequence makes a small table that the later sequences make grow (the rehash of ajTableResizeHashsize)
        items = [(items[0][0], items[0][1], items[0][2][:8] or 'A')] + [(nme + 'L', desc, (s * 20)[:rnd.randint(300, 1200)]) for nme, desc, s in (items[1:] or items)]
    if rnd.random() < 0.4:   # low complexity: many ties and repeated words
        alpha = 'ACGT'[:rnd.randint(1, 4)] if kind == 'n' else 'ACDEF'[:rnd.randint(1, 5)]
        items = [(nme, desc, ''.join(rnd.choice(alpha) for _ in range(rnd.randint(1, 300)))) for nme, desc, s in items]
    text, fmt = G.seq_file(rnd, d, items)
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '1', 0.1)
    if rnd.random() < 0.8: args += ['-wordsize', str(rnd.choice([1, 2, 2, 3, 3, 4, 5, 6, 8, 12, 0, -1]))]
    if rnd.random() < 0.4: args += ['-mincount', str(rnd.choice([1, 2, 3, 5, 0]))]
    if rnd.random() < 0.9: args += ['-outfile', rnd.choice(['out.txt', 'stdout'])]   # (no -outfile: the default file name BASE.program of the name of the first sequence)
    return args, {'in.dat': text}


def dr_entry(rnd, kind, fmt, name, desc, seq):
    dbs = ['EMBL', 'PDB', 'GO', 'InterPro', 'Pfam', 'PROSITE', 'MIM', 'HGNC', 'RefSeq', 'UniGene', 'embl', 'GO']
    lines = ''
    for _ in range(rnd.randint(0, 8)):
        db = rnd.choice(dbs)
        lines += 'DR   %s; %s%s.\n' % (db, rnd.choice(['P%05d' % rnd.randint(1, 99999), 'X%05d' % rnd.randint(1, 99), '1ABC', 'GO:%07d' % rnd.randint(1, 9999)]), rnd.choice(['', '; -', '; X-ray; 1.9 A; A/B=1-100', '; secondary']))
        if rnd.random() < 0.2: lines += lines.split('\n')[-2] + '\n'
    if fmt == 'embl':
        return 'ID   %s; SV 1; linear; genomic DNA; STD; HUM; %d BP.\nXX\nAC   X%05d;\nXX\nDE   %s\nXX\n%sXX\n     %s %9d\n//\n' % (name, len(seq), rnd.randint(1, 99999), desc, lines, ' '.join(seq[j:j + 10].lower() for j in range(0, min(len(seq), 60), 10)), len(seq))
    return 'ID   %s_HUMAN    Reviewed;    %d AA.\nAC   P%05d;\nDE   RecName: Full=%s;\n%sSQ   SEQUENCE   %d AA;  1000 MW;  ABCD CRC64;\n     %s\n//\n' % (name, len(seq), rnd.randint(1, 99999), desc, lines, len(seq), ' '.join(seq[j:j + 10] for j in range(0, min(len(seq), 60), 10)))


def gen_seqxref(rnd, i, d):
    kind = rnd.choice(['n', 'p', 'p'])
    fmt = rnd.choice(['embl', 'swiss', 'swiss', 'fasta', 'genbank'])
    if fmt == 'embl' and kind == 'p': kind = 'n'
    if fmt == 'swiss' and kind == 'n': kind = 'p'
    items = G.seqs_for(rnd, 'x', kind, rnd.randint(1, 4), (1, 90))
    if fmt in ('embl', 'swiss'):
        text = ''.join(dr_entry(rnd, kind, fmt, nme.upper()[:10], desc or 'some description', s.upper()) for nme, desc, s in items)
    else:
        text, fmt = G.seq_file(rnd, d, items, fmt)
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '1', 0.1)
    if rnd.random() < 0.9: args += ['-outfile', rnd.choice(['out.txt', 'stdout'])]   # (no -outfile: the default file name BASE.program of the name of the first sequence)
    return args, {'in.dat': text}


def gen_showorf(rnd, i, d):
    lens = rnd.choice([(1, 12), (1, 60), (30, 300)])
    items = G.seqs_for(rnd, 'x', 'n', 1, lens)
    text, fmt = G.seq_file(rnd, d, items)
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.15)
    if rnd.random() < 0.6: args += ['-frames', rnd.choice(['1,2,3,4,5,6', '1', '2,3', '4,5,6', '0', '1,4', '6,5,4', 'F1', 'R2,F3', '1,2,3', '*', '7', '1,1', ',1', '0,1', '1,0', 'F1,R1', 'None'])]
    if rnd.random() < 0.4: args += ['-table', rnd.choice(['0', '1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23', 'Standard', '7', 'Bacterial', '1,2'])]
    if rnd.random() < 0.5: args += ['-width', str(rnd.choice([10, 20, 30, 50, 60, 100, 5, 0, -3, 1000]))]
    for o in ['-noruler', '-noplabel', '-nonlabel']:
        if rnd.random() < 0.2: args.append(o)
    if rnd.random() < 0.9: args += ['-outfile', rnd.choice(['out.txt', 'stdout'])]   # (no -outfile: the default file name BASE.program of the name of the first sequence)
    return args, {'in.dat': text}, {'vg': 'underflow'}


def gen_prettyseq(rnd, i, d):
    lens = rnd.choice([(1, 12), (1, 60), (30, 300)])
    items = G.seqs_for(rnd, 'x', 'n', 1, lens)
    text, fmt = G.seq_file(rnd, d, items)
    n = len(items[0][2])
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.15)
    r = rnd.random()
    if r < 0.6:
        k = rnd.choice([1, 1, 1, 2, 3])
        parts = []
        for _ in range(k):
            a = rnd.randint(1, max(1, n)); b = rnd.randint(a, max(a, n))
            parts.append('%d-%d' % (a, b))
        args += ['-range', rnd.choice([',', ' ', ';']).join(parts)]
    elif r < 0.65:
        args += ['-range', rnd.choice(['0-3', '5-3', '1-1000', '3', 'a-b', '1,2,3', '1-3,2-9', '2-4,1-6', '1-2,1-2,1-3'])]
    if rnd.random() < 0.4: args += ['-table', rnd.choice(['0', '1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23', '7'])]
    if rnd.random() < 0.5: args += ['-width', str(rnd.choice([10, 20, 30, 50, 60, 100, 5, 0, -3, 1000]))]
    for o in ['-noruler', '-noplabel', '-nonlabel']:
        if rnd.random() < 0.2: args.append(o)
    if rnd.random() < 0.9: args += ['-outfile', rnd.choice(['out.txt', 'stdout'])]   # (no -outfile: the default file name BASE.program of the name of the first sequence)
    return args, {'in.dat': text}, {'vg': 'underflow'}


GENS = {'prettyseq': gen_prettyseq, 'showorf': gen_showorf, 'seqxref': gen_seqxref, 'listor': gen_listor, 'nthseqset': gen_nthseqset, 'infoseq': gen_infoseq, 'wordcount': gen_wordcount}

if __name__ == '__main__':
    prog = sys.argv.pop(1)
    L.main(prog, GENS[prog], 300)
