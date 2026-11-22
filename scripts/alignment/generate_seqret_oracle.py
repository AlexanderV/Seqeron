#!/usr/bin/env python3
"""A5-8c3g (B09 F248 ...) oracle: the sequence READERS of EMBOSS 6.6.0 (ajseqread.c, `seqinFormatDef`) through
    seqret -auto -feature -sequence in.dat [-sformat1 F] -outseq o.fa -osformat2 fasta -offormat2 gff3 -ofname2 o.gff
(the FASTA file of every sequence the reader returned and the GFF3 file of the feature tables of the sequences that have one; with a faked run date).
The input of a format is made in three ways: (a) a random set of sequences (names, accessions, descriptions, nucleotide / protein) written in the format by `seqret -osformat F` of the same installation
(a real writer output: header lines, line widths, numbering), (b) an entry with a feature table (EMBL / GenBank entries of generate_diffseq_features_oracle.py, protein GenPept entries) converted to the
format with `seqret -feature -osformat F` (the features of swiss, pir, gff2, gff3, genpept, refseq, refseqp, codata, embl, genbank inputs), (c) hand-made texts of the formats seqret cannot write (stockholm, pdb, ace,
experiment, biomart, dbid ...). A fraction of the inputs is damaged (a line deleted / duplicated / cut / replaced, CRLF line ends, blank lines) so that the failure and recovery paths of the readers are compared too.
Each input is read with `-sformat1 F`, with the `F::in.dat` form and (for the formats the program tries when none is named) with no format at all.
JSONL.GZ rows: id, fmt, how (the way the input was made), text, args (the arguments after the options of seqret), rc, fa, gff, err, hang (the program did not stop within the limit).
Usage: python3 scripts/alignment/generate_seqret_oracle.py OUT.jsonl.gz PER_FORMAT [seed] [format ...]      (needs EMBOSS 6.6.0 and gcc)
"""
import gzip, json, os, random, re, resource, shutil, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_diffseq_features_oracle as D2  # noqa: E402
import generate_featread_oracle as FR  # noqa: E402

NUC = 'ACGT'
PROT = 'ACDEFGHIKLMNPQRSTVWY'
WORDS = ['alpha', 'beta', 'kinase', 'putative', 'hypothetical protein', 'ribosomal', 'subunit', 'transport', 'regulator', 'of', 'the', 'binding', 'domain', 'zinc', 'finger', 'Homo sapiens', 'precursor']
NAMES = ['SEQA', 'AB123456', 'seq1', 'chr7', 'X56789', 'My_seq', 'sp_P12345', 'ctg.1', 'NP_000001', 'ABC_HUMAN', 'gene-9', 'a', 'S1234567890123456789012345']

# the formats and the kind of residues their writer and reader accept: n = nucleotide only, p = protein only, b = both
FORMATS = {
    'fasta': 'b', 'ncbi': 'b', 'pearson': 'b', 'gifasta': 'b', 'dbid': 'b', 'gcg': 'b', 'gcg8': 'b', 'embl': 'n', 'em': 'n', 'swiss': 'p', 'sw': 'p', 'swissprot': 'p', 'uniprot': 'p',
    'nbrf': 'b', 'pir': 'b', 'pdb': 'p', 'pdbseq': 'p', 'pdbnuc': 'n', 'pdbnucseq': 'n', 'fastq': 'n', 'fastq-sanger': 'n', 'fastq-illumina': 'n', 'fastq-solexa': 'n', 'sam': 'n',
    'genbank': 'n', 'gb': 'n', 'ddbj': 'n', 'refseq': 'n', 'refseqp': 'p', 'genpept': 'p', 'codata': 'b', 'strider': 'n', 'clustal': 'b', 'aln': 'b', 'phylip': 'b', 'phylipnon': 'b',
    'ace': 'n', 'consed': 'n', 'acedb': 'b', 'msf': 'b', 'hennig86': 'b', 'jackknifer': 'b', 'nexus': 'b', 'paup': 'b', 'treecon': 'b', 'mega': 'b', 'igstrict': 'b', 'iguspto': 'b', 'ig': 'b',
    'staden': 'b', 'textonly': 'b', 'plain': 'b', 'asis': 'b', 'text': 'b', 'raw': 'b', 'gff2': 'b', 'gff3': 'b', 'gff': 'b', 'stockholm': 'b', 'pfam': 'b', 'selex': 'b', 'fitch': 'b', 'biomart': 'b',
    'mase': 'b', 'experiment': 'n', 'gde': 'b', 'unknown': 'b',
}
# the writer name that makes the format (None: hand-made)
WRITER = {f: f for f in FORMATS}
WRITER.update({'gcg8': 'gcg', 'em': 'embl', 'sw': 'swiss', 'swissprot': 'swiss', 'uniprot': 'swiss', 'pir': 'nbrf', 'gb': 'genbank', 'ddbj': 'genbank', 'aln': 'clustal', 'consed': None, 'ace': None, 'paup': 'nexus',
               'asis': 'text', 'plain': 'text', 'textonly': 'text', 'unknown': 'text', 'gff': 'gff3', 'pfam': None, 'stockholm': None, 'pdb': None, 'pdbseq': None, 'pdbnuc': None, 'pdbnucseq': None, 'experiment': None,
               'biomart': None, 'dbid': 'fasta', 'pearson': 'fasta', 'gifasta': 'fasta', 'ncbi': 'ncbi', 'fastq-sanger': 'fastq-sanger', 'refseq': 'genbank', 'refseqp': 'genbank', 'genpept': 'genbank',
               'gcg': 'gcg', 'ig': 'ig', 'igstrict': 'ig', 'iguspto': 'ig', 'raw': 'raw'})
# the formats whose reference runs are checked with valgrind (memerr: the run reads or writes memory it does not own)
VG_FORMATS = {'pdb', 'pdbseq', 'pdbnuc', 'pdbnucseq', 'experiment', 'phylip', 'phylipnon', 'nexus', 'paup', 'mega', 'selex', 'stockholm', 'pfam', 'jackknifer', 'hennig86'}
# formats the program tries when none is named (the Try flag of seqinFormatDef)
TRY = {'gcg', 'embl', 'swiss', 'nbrf', 'pdb', 'fasta', 'fastq', 'sam', 'genbank', 'codata', 'strider', 'clustal', 'phylip', 'ace', 'consed', 'acedb', 'msf', 'hennig86', 'jackknifer', 'nexus', 'treecon', 'mega',
       'igstrict', 'iguspto', 'gff2', 'gff3', 'stockholm', 'pfam', 'fitch', 'experiment', 'gde', 'raw'}


def rand_seq(rnd, kind, n=None):
    if kind == 'n' or (kind == 'b' and rnd.random() < 0.6):
        alpha, extra = NUC, 'RYNKMSW'
    else:
        alpha, extra = PROT, 'BZX'
    n = n or rnd.choice([1, 5, 12, 33, 50, 59, 60, 61, 80, 100, 130, 200, 301])
    s = ''.join(rnd.choice(alpha) for _ in range(n))
    if rnd.random() < 0.15:
        s = ''.join(rnd.choice(extra) if rnd.random() < 0.1 else c for c in s)
    return s


def rand_name(rnd):
    nm = rnd.choice(NAMES)
    if rnd.random() < 0.3:
        nm += ''.join(rnd.choice('abcdefghijklmnopqrstuvwxyz0123456789_.-') for _ in range(rnd.choice([1, 4, 9, 17])))
    return nm


def rand_desc(rnd):
    if rnd.random() < 0.15:
        return ''
    return ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(1, 9)))


def rand_fasta(rnd, kind, nseq=None, gi=0.15):
    out = ''
    nuc = kind == 'n' or (kind == 'b' and rnd.random() < 0.6)
    for i in range(nseq or rnd.choice([1, 1, 1, 2, 3, 4])):
        s = rand_seq(rnd, 'n' if nuc else 'p')
        name = rand_name(rnd) + ('' if i == 0 else str(i))
        acc = rnd.choice(['', '', 'AB123456', 'P12345', 'NM_000001', 'U12345.2', 'X56789'])
        desc = rand_desc(rnd)
        head = '>' + name + (' ' + acc if acc else '') + (' ' + desc if desc else '')
        if rnd.random() < gi:
            head = '>gi|%d|gb|%s|%s' % (rnd.randint(1, 999999), rnd.choice(['AB123456.1', 'U12345.2']), desc)
        w = rnd.choice([60, 60, 60, 70, 80, 50, 1000])
        out += head + '\n' + ''.join(s[j:j + w] + '\n' for j in range(0, len(s), w))
    return out


def seqret(d, args, timeout=5, env=None):
    hang = False
    def limit():
        resource.setrlimit(resource.RLIMIT_FSIZE, (4 << 20, 4 << 20))
    proc = subprocess.Popen(['seqret'] + args, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=limit, env=env)
    try:
        _, errb = proc.communicate(timeout=timeout)
        rc, err = proc.returncode, errb.decode('latin-1')
        if rc == -25:
            hang = True
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.communicate()
        rc, err, hang = -9, '', True
    return rc, err, hang


def damaged(rnd, text):
    lines = text.split('\n')
    r = rnd.random()
    if r < 0.2 and len(lines) > 2:
        del lines[rnd.randrange(len(lines))]
    elif r < 0.35:
        i = rnd.randrange(len(lines)); lines.insert(i, lines[i])
    elif r < 0.5:
        lines.insert(rnd.randrange(len(lines) + 1), rnd.choice(['', 'random text', '#comment', '//', '>x y', 'ID   X; 1 BP.', ';', '\t\t', 'LOCUS       X 1 bp']))
    elif r < 0.65 and lines:
        i = rnd.randrange(len(lines)); lines[i] = lines[i][:rnd.randint(0, len(lines[i]))]
    elif r < 0.75:
        lines = [l + '\r' for l in lines]
    elif r < 0.85 and lines:
        i = rnd.randrange(len(lines)); l = lines[i]
        if l: j = rnd.randrange(len(l)); lines[i] = l[:j] + rnd.choice('0123456789 *-.?#xX') + l[j + 1:]
    elif r < 0.92:
        lines = lines[:rnd.randint(1, len(lines))]
    else:
        lines = [l for l in lines if l.strip()] if rnd.random() < 0.5 else [''] + lines + ['']
    return '\n'.join(lines)


# ---------------------------------------------------------------------------------------------------------------- hand-made texts


def hand_made(rnd, fmt):
    s = rand_seq(rnd, FORMATS[fmt] if FORMATS[fmt] != 'b' else rnd.choice('np'))
    name = rand_name(rnd)
    desc = rand_desc(rnd)
    if fmt in ('stockholm', 'pfam'):
        rows = [(name + str(i), s if rnd.random() < 0.8 else s[:-1] + '-') for i in range(rnd.choice([1, 2, 3]))]
        o = '# STOCKHOLM 1.0\n' + ('#=GF ID   %s\n#=GF DE   %s\n' % (name, desc) if rnd.random() < 0.6 else '')
        for nm, r in rows:
            o += '%-20s %s\n' % (nm, r)
        if rnd.random() < 0.4: o += '#=GC SS_cons %s\n' % ('.' * len(s))
        return o + '//\n'
    if fmt in ('pdb', 'pdbseq', 'pdbnuc', 'pdbnucseq'):
        nuc = fmt.startswith('pdbnuc')
        aa3 = {'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'}
        s = ''.join(rnd.choice(NUC if nuc else PROT) for _ in range(rnd.randint(3, 40)))
        o = 'HEADER    %-40s%-9s   %-4s\nTITLE     %s\n' % (rnd.choice(['HYDROLASE', 'TRANSFERASE/DNA', 'DNA', 'PROTEIN']), '01-JAN-00', rnd.choice(['1ABC', '2XYZ', '9ZZ9', '1A2B']), desc.upper())
        for ch in rnd.choice(['A', 'AB', 'ABC', 'B']):
            if not fmt.endswith('seq') and rnd.random() < 0.15: o += 'MODEL     %4d\n' % rnd.randint(1, 3)
            if fmt.endswith('seq'):
                toks = [('D' + c) if nuc else aa3[c] for c in s]
                o += ''.join('SEQRES %3d %s %4d  %s\n' % (k // 13 + 1, ch, len(s), ' '.join(toks[k:k + 13])) for k in range(0, len(toks), 13))
            else:
                for i, c in enumerate(s):
                    res = (' D' + c) if nuc else aa3[c]
                    o += 'ATOM  %5d  %-3s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n' % (i * 4 + 1, 'CA' if not nuc else "C1'", res, ch, i + 1, rnd.random() * 9, rnd.random() * 9, rnd.random() * 9)
                if rnd.random() < 0.9: o += 'TER   %5d      %3s %s%4d\n' % (len(s) * 4 + 1, res, ch, len(s))
        return o + 'END\n'
    if fmt in ('ace', 'consed'):
        o = 'AS 1 %d\n\nCO %s %d 1 0 U\n' % (rnd.randint(1, 3), name, len(s))
        o += ''.join(s[j:j + 50] + '\n' for j in range(0, len(s), 50)) + '\nBQ\n' + ' '.join('30' for _ in s) + '\n\nAF r1 U 1\n\nRD r1 %d 0 0\n%s\n\nQA 1 %d 1 %d\n' % (len(s), s, len(s), len(s))
        return o
    if fmt == 'experiment':
        o = ''
        for e in range(rnd.choice([1, 1, 2, 3])):
            sq = rand_seq(rnd, 'n', rnd.choice([1, 5, 12, 40, 70]))
            o += 'ID   %s\n' % (name + str(e) if e else name)
            for _ in range(rnd.randint(0, 6)):
                k = rnd.choice(['EN', 'EX', 'AC', 'SV', 'DE', 'AV', 'KW', 'OS', 'OC', 'TN', 'FT', 'CC', 'XX', 'LN'])
                if k == 'AV':
                    vals = ' '.join(rnd.choice(['30', '5', '-1', '12,30', 'x', '7,y', '99']) for _ in range(rnd.choice([len(sq), len(sq), len(sq) - 1 if len(sq) > 1 else 1, len(sq) + 2])))
                    o += 'AV   %s\n' % vals
                elif k == 'AC':
                    o += 'AC   %s;\n' % rnd.choice(['AB123456', 'X56789; Y12345', 'P12345'])
                elif k == 'SV':
                    o += 'SV   %s.%d\n' % (rnd.choice(['AB123456', 'X1']), rnd.randint(1, 3))
                elif k in ('KW', 'OS', 'OC'):
                    o += '%s   %s\n' % (k, rnd.choice(['Homo sapiens (human).', 'Eukaryota; Metazoa; Chordata.', 'kinase; transport.']))
                elif k == 'FT':
                    o += 'FT   CDS             1..%d\n' % len(sq)
                else:
                    o += '%s   %s\n' % (k, rand_desc(rnd) or 'x')
            o += rnd.choice(['SQ\n', 'SQ   %d BP\n' % len(sq)]) + ''.join('     ' + sq[j:j + 60] + '\n' for j in range(0, len(sq), 60)) + '//\n'
        return o
    if fmt == 'biomart':
        return '%s\t%s\n' % (name, s) if rnd.random() < 0.5 else '>%s\n%s\n' % (name, s)
    return '>%s\n%s\n' % (name, s)


# ---------------------------------------------------------------------------------------------------------------- feature-bearing entries

PROT_FEATS = ['Region', 'Site', 'Bond', 'Domain', 'Chain', 'Signal', 'Transit']


def genpept_entry(rnd, kind):
    """a GenPept / RefSeq-protein entry with a feature table (kind 'genpept' or 'refseqp')"""
    n = rnd.randint(30, 200)
    s = rand_seq(rnd, 'p', n)
    name = rnd.choice(['NP_000001', 'AAA12345', 'XP_00123', 'P12345'])
    o = 'LOCUS       %-16s %d aa            linear   PRI 01-JAN-2000\nDEFINITION  %s.\nACCESSION   %s\nVERSION     %s.1  GI:%d\n' % (name, n, rand_desc(rnd) or 'x', name, name, rnd.randint(1, 99999))
    if rnd.random() < 0.7: o += 'DBSOURCE    REFSEQ: accession NM_000001.1\n'
    o += 'KEYWORDS    RefSeq.\nSOURCE      Homo sapiens (human)\n  ORGANISM  Homo sapiens\n            Eukaryota; Metazoa.\n'
    if rnd.random() < 0.5: o += 'REFERENCE   1  (residues 1 to %d)\n  AUTHORS   Smith,J.\n  TITLE     A title\n  JOURNAL   Nature 1 (1), 1-2 (2000)\n' % n
    o += 'FEATURES             Location/Qualifiers\n     source          1..%d\n                     /organism="Homo sapiens"\n' % n
    o += ''.join(l + '\n' for l in FR.ft_lines(rnd, n, 'genbank', FR.TYPES['protein' if kind == 'genpept' else 'refseqp'], protein=True))
    o += 'ORIGIN\n' + ''.join('%9d %s\n' % (i + 1, ' '.join(s[i + j:i + j + 10].lower() for j in range(0, 60, 10) if s[i + j:i + j + 10])) for i in range(0, n, 60))
    return o + '//\n'


def feature_source(rnd, fmt):
    """(text, format of the text, kind of residues) of an entry with features that is converted to `fmt` by seqret"""
    k = FORMATS[fmt]
    if k == 'p' or fmt in ('nbrf', 'pir') or (k == 'b' and rnd.random() < 0.35):
        return genpept_entry(rnd, rnd.choice(['genpept', 'refseqp'])), rnd.choice(['genpept', 'refseqp']), 'p'
    n = rnd.randint(40, 330)
    seq = ''.join(rnd.choice(NUC) for _ in range(n))
    ef = rnd.choice(['embl', 'genbank'])
    D2.SPECIAL = rnd.random() < 0.3
    return D2.entry(rnd, ef, rnd.choice(['SEQA', 'AB123456', 'gb_a', 'X56789']), rand_desc(rnd) or 'test', seq, D2.feature_lines(rnd, ef, n)), ef, 'n'


def make_input(rnd, d, fmt):
    """(text, how): the input of the format; the writers of seqret 6.6 leave some texts unreadable (jackknifer: the title quote is not closed) so most of those are repaired by hand"""
    text, how = make_input0(rnd, d, fmt)
    if fmt == 'jackknifer' and rnd.random() < 0.85 and text.startswith("'"):
        first, _, rest = text.partition('\n')
        text, how = first.rstrip() + "'\n" + rest, how + '+closed'
    if fmt in ('nexus', 'paup') and rnd.random() < 0.5:
        text, how = nexus_rich(rnd, text), 'nexus-rich'
    return text, how


def nexus_rich(rnd, text):
    """a hand-made NEXUS file: taxa block with labels, comments, quoted labels, sequential or interleaved matrix, other blocks, gap / datatype / missing / format variants"""
    nuc = rnd.random() < 0.6
    alpha = NUC if nuc else PROT
    ntax = rnd.choice([1, 2, 3, 4])
    nchar = rnd.choice([5, 12, 33, 60, 80])
    names = [rand_name(rnd).replace(' ', '_') for _ in range(ntax)]
    if rnd.random() < 0.15 and ntax > 1:
        names[1] = names[0]
    seqs = [''.join(rnd.choice(alpha + '-') for _ in range(nchar)) for _ in range(ntax)]
    o = '#NEXUS\n'
    if rnd.random() < 0.5: o += '[comment line]\n'
    taxa = rnd.random() < 0.6
    lab = lambda n: ("'%s'" % n) if rnd.random() < 0.2 else n
    if taxa:
        o += 'begin taxa;\n  dimensions ntax=%d;\n  taxlabels %s;\nend;\n' % (ntax, ' '.join(lab(n) for n in names))
    blk = rnd.choice(['data', 'characters', 'DATA', 'Characters'])
    o += 'begin %s;\n' % blk
    if taxa and blk.lower() == 'data' or rnd.random() < 0.8:
        o += '  dimensions %sntax=%d nchar=%d;\n' % ('newtaxa ' if not taxa and rnd.random() < 0.3 else '', ntax, nchar)
    else:
        o += '  dimensions nchar=%d;\n' % nchar
    fm = 'datatype=%s' % rnd.choice(['dna', 'DNA', 'protein', 'rna', 'nucleotide', 'standard', 'continuous', 'mixed'] if rnd.random() < 0.3 else ['dna' if nuc else 'protein'])
    if rnd.random() < 0.7: fm += ' missing=%s' % rnd.choice('?Nx')
    gap = rnd.choice(['-', '.', '-', 'x'])
    if rnd.random() < 0.8: fm += ' gap=%s' % gap
    inter = rnd.random() < 0.5
    if inter: fm += ' interleave'
    if rnd.random() < 0.3: fm += rnd.choice([' nolabels', ' labels', ' respectcase', " symbols=\"ACGT\"", ' transpose'])
    o += '  format %s;\n' % fm
    if rnd.random() < 0.2: o += '  charlabels %s;\n' % ' '.join('c%d' % i for i in range(3))
    if rnd.random() < 0.2: o += '  [a comment\n   over two lines]\n'
    o += '  matrix\n'
    w = rnd.choice([10, 20, 50, nchar])
    if inter:
        for j in range(0, nchar, w):
            for n, sq in zip(names, seqs):
                o += '    %s %s\n' % (lab(n) if rnd.random() < 0.8 else n, sq[j:j + w])
            if rnd.random() < 0.5: o += '\n'
    else:
        for n, sq in zip(names, seqs):
            o += '    %s\n' % lab(n)
            for j in range(0, nchar, w): o += '    %s\n' % sq[j:j + w]
    o += ';\nend;\n'
    if rnd.random() < 0.4: o += 'begin sets;\n  charset a = 1-3;\n  taxset t = 1 2;\nend;\n'
    if rnd.random() < 0.3: o += 'begin assumptions;\n  options deftype=unord polytcount=minsteps;\nend;\n'
    if rnd.random() < 0.3: o += 'begin trees;\n  tree t1 = (a,b);\nend;\n'
    if rnd.random() < 0.2: o += 'begin notes;\n  text taxon=1 text="x";\nend;\n'
    if rnd.random() < 0.15: o += 'begin unknownblock;\n  foo bar;\nend;\n'
    return o


def make_input0(rnd, d, fmt):
    """(text, how)"""
    writer = WRITER[fmt]
    r = rnd.random()
    if writer is None:
        return hand_made(rnd, fmt), 'hand'
    feature_formats = ('embl', 'em', 'genbank', 'gb', 'ddbj', 'refseq', 'refseqp', 'genpept', 'swiss', 'sw', 'swissprot', 'uniprot', 'nbrf', 'pir', 'gff2', 'gff3', 'gff', 'codata', 'gcg', 'gcg8')
    if fmt in feature_formats and r < 0.55:
        src, srcfmt, kind = feature_source(rnd, fmt)
        if fmt in ('genpept', 'refseqp', 'refseq', 'genbank', 'gb', 'ddbj') and kind == ('p' if fmt in ('genpept', 'refseqp') else 'n') and srcfmt in ('genpept', 'refseqp', 'embl', 'genbank'):
            if rnd.random() < 0.5: return src, 'entry'
        open(os.path.join(d, 'src.dat'), 'w', newline='\n', encoding='latin-1').write(src)
        for nme in ('in.dat', 'x.gff'):
            if os.path.exists(os.path.join(d, nme)): os.remove(os.path.join(d, nme))
        rc, err, hang = seqret(d, ['-auto', '-feature', '-sequence', 'src.dat', '-sformat1', srcfmt, '-osformat2', writer, '-outseq', 'in.dat', '-offormat2', 'gff3', '-ofname2', 'x.gff'])
        q = os.path.join(d, 'in.dat')
        if rc == 0 and not hang and os.path.exists(q) and os.path.getsize(q) > 0 and 'Error' not in err:
            return open(q, encoding='latin-1', newline='').read(), 'features'
    kind = FORMATS[fmt]
    for _ in range(5):
        fa = rand_fasta(rnd, kind, gi=0.9 if fmt == 'gifasta' else 0.15)
        open(os.path.join(d, 'src.fa'), 'w', newline='\n', encoding='latin-1').write(fa)
        if os.path.exists(os.path.join(d, 'in.dat')): os.remove(os.path.join(d, 'in.dat'))
        if writer == 'fasta' and fmt in ('fasta', 'pearson', 'gifasta', 'dbid'):
            return fa, 'direct'
        rc, err, hang = seqret(d, ['-auto', '-sequence', 'src.fa', '-osformat2', writer, '-outseq', 'in.dat', '-ossingle2' if False else '-auto'])
        q = os.path.join(d, 'in.dat')
        if rc == 0 and not hang and os.path.exists(q) and os.path.getsize(q) > 0:
            return open(q, encoding='latin-1', newline='').read(), 'writer'
    return hand_made(rnd, fmt), 'hand'


def main(out, per, seed, formats):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        gen = os.path.join(d, 'g'); os.mkdir(gen)
        for fmt in formats:
            made = 0
            while made < per:
                text, how = make_input(rnd, gen, fmt)
                if rnd.random() < 0.3: text, how = damaged(rnd, text), how + '+damaged'
                if rnd.random() < 0.03: text = rnd.choice(['', '\n\n', 'just some text\nmore text\n', '#only a comment\n', '//\n'])
                if rnd.random() < 0.05: text = text.rstrip('\n')
                r = rnd.random()
                spec = 'in.dat'
                args = ['-sformat1', fmt]
                if r < 0.25: spec = '%s::in.dat' % fmt; args = []
                elif r < 0.35 and fmt in TRY: args = []
                elif r < 0.40: args = ['-sformat1', fmt.upper()]
                if fmt == 'asis' and rnd.random() < 0.4:
                    spec = 'asis::' + ''.join(rnd.choice('ACGTacgtNRYMK*-.~?#+ 0123456789xX') for _ in range(rnd.randint(1, 80)))
                    args = []
                if rnd.random() < 0.08: args += ['-firstonly']
                if rnd.random() < 0.06: args += [rnd.choice(['-slower1', '-supper1', '-sprotein1', '-snucleotide1'])]
                feature = rnd.random() < 0.8
                for nme in os.listdir(work): os.remove(os.path.join(work, nme))
                open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
                epoch = rnd.randint(978307200, 2051222400)
                env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
                cmd = ['-auto'] + (['-feature'] if feature else []) + ['-sequence', spec] + args + ['-outseq', 'o.fa', '-osformat2', 'fasta'] + (['-offormat2', 'gff3', '-ofname2', 'o.gff'] if feature else [])
                rc, err, hang = seqret(work, cmd, env=env)
                def rd(nme):
                    q = os.path.join(work, nme)
                    return open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
                if err.startswith('Read and write (return) sequences\n'): err = err[len('Read and write (return) sequences\n'):]
                memerr = False
                if fmt in VG_FORMATS and not hang and shutil.which('valgrind'):
                    # a reference run that reads or writes freed / uninitialised memory has no defined output: valgrind -q with an exit code (the runs of the formats whose readers do so are refused by the port)
                    vgargs = [a.replace('o.fa', 'ovg.fa').replace('o.gff', 'ovg.gff') for a in cmd]
                    try:
                        vg = subprocess.run(['valgrind', '-q', '--error-exitcode=99', 'seqret'] + vgargs, cwd=work, capture_output=True, timeout=120, env=dict(os.environ))
                        # a program that dies by itself keeps its own exit code, so the error reports of valgrind are looked for in its output
                        memerr = vg.returncode == 99 or re.search(rb'^==\d+== (Invalid|Conditional|Use of|Syscall|Mismatched|Source and destination|Process terminating)', vg.stderr, re.M) is not None
                    except subprocess.TimeoutExpired:
                        memerr = True
                rows.append({'id': len(rows), 'fmt': fmt, 'how': how, 'memerr': memerr, 'text': text, 'spec': spec, 'args': args, 'feature': feature, 'epoch': epoch, 'rc': rc, 'fa': rd('o.fa'), 'gff': rd('o.gff'), 'err': err, 'hang': hang})
                made += 1
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    c = Counter((r['fmt'], 'hang' if r['hang'] else 'fa' if r['fa'] else 'none') for r in rows)
    print(len(rows), 'rows;', sorted(c.items()))


if __name__ == '__main__':
    fm = sys.argv[4:] or list(FORMATS)
    main(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 20261622, fm)
