#!/usr/bin/env python3
"""A5-8c3e (B09 F241) oracle: the feature READERS of EMBOSS 6.6.0 (ajfeatread.c) through `featcopy -features FILE -outfeat out -offormat debug`, which reads every feature table of a file in any
input format and prints the full internal state (ajFeattablePrint: type, source, score, strand, positions, flags, group, exon, tags, sub-features).
Formats: embl, genbank (entries with up to 3 tables), refseq, genpept, refseqp (feature lines of five blanks), gff2 (versions 1 / 2 / 3, ##Type, ##sequence-region, groups, featflags), gff3 (ID / Parent /
Name ..., percent-encoded values, several lines of one ID, several tables), gff3 written by EMBOSS before 6.4 (gff3emboss), swiss, pir, dasgff (XML), and damaged / mixed files.
Each file is read with no format (the program tries every format), with the `format::file` prefix, and with -fbegin1 / -fend1.
JSONL.GZ rows: id, kind, text, ufo, extra (the other arguments), rc, debug (the file the program wrote), err (stderr), hang (the program did not stop within the time limit: GFF2 loops forever on a group
field its tag pattern cannot parse)
Usage: python3 scripts/alignment/generate_featread_oracle.py OUT.jsonl.gz [n=1000] [seed]      (needs EMBOSS 6.6.0)
"""
import gzip, json, os, random, re, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_diffseq_features_oracle as D2  # noqa: E402

DATA = '/usr/share/EMBOSS/data/'
WORDS = D2.WORDS


def vocab_types(name):
    out = []
    for raw in open(DATA + 'Efeatures.' + name, encoding='latin-1'):
        line = ' '.join(raw.split())
        if not line or line.startswith('#') or line.startswith('/') or line[0] in 'mM-+' and len(line) > 1 and line[1] == '/':
            continue
        w = line.split(' ')[0]
        if w and not w.startswith('#'):
            out.append(w)
    return out


TYPES = {n: vocab_types(n) for n in ('embl', 'gff3', 'gff2', 'swiss', 'pir', 'refseqp', 'gff3protein')}
TYPES['protein'] = TYPES['gff3protein'] + TYPES['swiss']


def seqname(rnd):
    return rnd.choice(['SEQA', 'AB123456', 'chr7', 'ctg1', 'X56789', 'seq1', 'gi|123|gb|AB1.1|'.replace('|', '_')])


# ---------------------------------------------------------------------------------------------------------------- embl / genbank / refseq / genpept / refseqp


def ft_lines(rnd, n, fmt, types, protein=False):
    """feature table lines (embl: FT; genbank / refseq / genpept / refseqp: five blanks)"""
    lines = []
    for _ in range(rnd.choice([1, 1, 2, 3, 4, 6])):
        ty = rnd.choice(types)
        if rnd.random() < 0.08:
            ty = rnd.choice([ty.upper(), ty.lower(), 'weird_key', 'Source'])
        if protein:
            a = rnd.randint(1, n); b = rnd.randint(a, n)
            loc = rnd.choice(['%d..%d' % (a, b), str(a), '<%d..%d' % (a, b), '%d..>%d' % (a, b), 'join(%d..%d,%d..%d)' % (a, a + 1, min(n, b), n), '%d^%d' % (a, a + 1), '(%d.%d)..%d' % (a, a + 1, b), 'bad', 'complement(%d..%d)' % (a, b)])
            quals = []
            for _ in range(rnd.choice([0, 1, 2, 3])):
                quals.append(rnd.choice([('note', '"%s"' % ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(1, 8)))), ('region_name', '"%s"' % rnd.choice(WORDS)), ('site_type', rnd.choice(['active', 'binding', 'other'])),
                                         ('product', '"%s"' % rnd.choice(WORDS)), ('db_xref', '"%s"' % rnd.choice(['taxon:9606', 'GI:1', 'bad'])), ('codon_start', '2'), ('experiment', '"x"'), ('weird', '"w"'), ('pseudo', None),
                                         ('calculated_mol_wt', str(rnd.randint(100, 9999))), ('bond_type', rnd.choice(['disulfide', 'glycal', 'x']))]))
        else:
            loc = D2.locations(rnd, n)
            quals = D2.qualifiers(rnd)
        key = ty if len(ty) < 16 else ty + ' '
        if fmt == 'embl':
            lines.append('FT   %-16s%s' % (key, loc))
            for name, val in quals:
                text = '/%s' % name + ('' if val is None else '=%s' % val)
                for w in D2.wrap(text, 58): lines.append('FT' + ' ' * 19 + w)
        else:
            lines.append('     %-16s%s' % (key, loc))
            for name, val in quals:
                text = '/%s' % name + ('' if val is None else '=%s' % val)
                for w in D2.wrap(text, 58): lines.append(' ' * 21 + w)
    return lines


def embl_file(rnd, fmt):
    D2.SPECIAL = rnd.random() < 0.35
    out = ''
    for k in range(rnd.choice([1, 1, 1, 2, 3])):
        n = rnd.randint(60, 400)
        name = seqname(rnd) + ('' if k == 0 else str(k))
        feats = ft_lines(rnd, n, 'embl' if fmt == 'embl' else 'genbank', TYPES['embl'])
        if fmt == 'embl':
            out += 'ID   %s; SV 1; linear; genomic DNA; STD; UNC; %d BP.\nXX\nFH   Key             Location/Qualifiers\nFH\n' % (name, n)
            if rnd.random() < 0.8: out += 'FT   source          1..%d\n' % n
            out += ''.join(l + '\n' for l in feats) + 'XX\n//\n'
        else:
            out += 'LOCUS       %s %d bp    DNA     linear   UNK 01-JAN-2000\n' % (name, n)
            if rnd.random() < 0.5: out += 'DEFINITION  test.\n'
            out += 'FEATURES             Location/Qualifiers\n'
            if rnd.random() < 0.8: out += '     source          1..%d\n' % n
            out += ''.join(l + '\n' for l in feats) + 'ORIGIN\n        1 acgt\n//\n'
    return out


def five_blank_file(rnd, kind):
    D2.SPECIAL = kind == 'refseq' and rnd.random() < 0.3
    n = rnd.randint(40, 300)
    types = {'refseq': TYPES['embl'], 'genpept': TYPES['protein'], 'refseqp': TYPES['refseqp']}[kind]
    lines = ft_lines(rnd, n, 'genbank', types, protein=kind != 'refseq')
    head = ''
    if rnd.random() < 0.5: head = 'FEATURES             Location/Qualifiers\n'
    return head + ''.join(l + '\n' for l in lines)


# ---------------------------------------------------------------------------------------------------------------- gff


def gff_group(rnd, ver):
    parts = []
    r = rnd.random()
    if ver == 1 or r < 0.1:
        return rnd.choice(['', 'Sequence "%s"' % rnd.choice(['A.1', 'B', 'x.12']), '"grp.%d"' % rnd.randint(1, 5), 'grp 3', 'A.4'])
    if r < 0.5: parts.append('Sequence "%s.%d"' % (rnd.choice(['A', 'B', 'seq']), rnd.randint(1, 4)))
    if rnd.random() < 0.5: parts.append('Note "%s"' % ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(1, 5))))
    if rnd.random() < 0.3: parts.append('gene "%s"' % rnd.choice(['abc', 'recA']))
    if rnd.random() < 0.2: parts.append('db_xref "%s"' % rnd.choice(['taxon:9606', 'bad']))
    if rnd.random() < 0.25: parts.append('FeatFlags "%s"' % rnd.choice(['0x1', '0x2 start2:5', '0x10 end2:9', '0x1000 remoteid:AB1.1', 'label:foo', '0x20', 'start2:3 end2:8', '0x200', '0x800']))
    if rnd.random() < 0.15: parts.append('ID %s' % rnd.choice(['grp.1', 'x.5']))
    if rnd.random() < 0.1: parts.append('weirdtag %s' % rnd.choice(['a', '"b c"', '1']))
    if rnd.random() < 0.05: parts.append(rnd.choice(['=', 'bad', ' ;', 'x=']))        # a group field the tag pattern cannot parse: the program loops forever
    return '; '.join(parts)


def gff2_file(rnd):
    ver = rnd.choice([1, 2, 2, 2, 3])
    lines = []
    if rnd.random() < 0.7: lines.append('##gff-version %d' % ver)
    if rnd.random() < 0.4: lines.append('##Type %s %s' % (rnd.choice(['DNA', 'Protein', 'dna', 'protein']), seqname(rnd)) if rnd.random() < 0.8 else '##Type DNA')
    if rnd.random() < 0.4: lines.append('##sequence-region %s %d %d' % (seqname(rnd).replace('_', 'x'), rnd.randint(1, 5), rnd.randint(50, 400)))
    if rnd.random() < 0.3: lines.append('# a comment line')
    types = TYPES['gff2'] + TYPES['embl'][:20] + ['weird', 'SO:0000704']
    for _ in range(rnd.choice([1, 2, 3, 5, 8])):
        s = rnd.randint(1, 300); e = rnd.randint(s, 340)
        if rnd.random() < 0.1: s, e = e, s
        fields = [rnd.choice(['SEQ', 'chr1', 'x']), rnd.choice(['src', 'EMBL', 'prog', '.']), rnd.choice(types), str(s), str(e), rnd.choice(['.', '5.5', '0', '-3.25', '1e2', 'x', '12']),
                  rnd.choice(['+', '-', '.', '+', '-']), rnd.choice(['.', '0', '1', '2', '.'])]
        g = gff_group(rnd, ver)
        if g or rnd.random() < 0.5: fields.append(g)
        line = '\t'.join(fields)
        if rnd.random() < 0.05: line = '\t'.join(fields[:rnd.randint(1, 7)])
        lines.append(line)
        if rnd.random() < 0.1: lines.append('')
    return '\n'.join(lines) + '\n'


def gff3_file(rnd, old=False):
    tables = rnd.choice([1, 1, 1, 2, 3])
    out = ''
    for t in range(tables):
        lines = []
        if old: lines += ['##gff-version 3', '#!Date 2005-01-01', '#!Type DNA', '#!Source-version EMBOSS %s' % rnd.choice(['6.0.0.0', '6.2.0.0', '6.3.1.0', '6.4.0.0'])]
        else:
            lines.append('##gff-version 3')
            if rnd.random() < 0.7: lines.append('##sequence-region %s %d %d' % ('seq%d' % t, 1, rnd.randint(100, 400)))
            if rnd.random() < 0.3: lines.append('#!Type %s' % rnd.choice(['DNA', 'Protein']))
            if rnd.random() < 0.2: lines.append('#!Source-version EMBOSS %s' % rnd.choice(['6.6.0.0', '6.1.0', '6.3.0']))
        ids = []
        for i in range(rnd.choice([1, 2, 4, 6, 9])):
            s = rnd.randint(1, 300); e = rnd.randint(s, 340)
            attrs = []
            r = rnd.random()
            if r < 0.5:
                myid = rnd.choice(['g%d' % i, 'g1', 'shared', 'm%d' % i]) if not old else 'x.%d' % i
                attrs.append('ID=%s' % myid); ids.append(myid)
            if ids and rnd.random() < 0.35 and not old: attrs.append('Parent=%s' % rnd.choice(ids + ['missing']))
            if rnd.random() < 0.4: attrs.append('Name=%s' % rnd.choice(['abc', 'recA', 'a%20b', 'x%3Bq']))
            if rnd.random() < 0.2: attrs.append('Alias=%s' % rnd.choice(['a1,a2', 'zz']))
            if rnd.random() < 0.3: attrs.append('Note=%s' % rnd.choice(['hello', 'hello%20world', 'a%2Cb', 'x,y,z']))
            if rnd.random() < 0.2: attrs.append('Dbxref=%s' % rnd.choice(['taxon:9606', 'GI:1,GI:2']))
            if rnd.random() < 0.2: attrs.append('product=%s' % rnd.choice(['foo', 'bar baz']))
            if rnd.random() < 0.15: attrs.append('featflags=%s' % rnd.choice(['0x1', 'start_before:true', 'end_after:true,between:true', 'start2:4,end2:9', 'type:gene', 'label:lab', 'start_unsure:true', '0x200', 'bogus:1']))
            if rnd.random() < 0.15: attrs.append('Is_circular=true')
            if rnd.random() < 0.1: attrs.append('Target=x %d %d +' % (s, e))
            if rnd.random() < 0.1: attrs.append(rnd.choice(['Custom=val', 'lower=v', 'weird_tag=%41b', 'Gap=M3 I1']))
            types = TYPES['gff3'] + ['gene', 'mRNA', 'exon', 'CDS', 'region', 'weird_type', 'SO:0000704']
            fields = [rnd.choice(['seq%d' % t, 'seq%d' % t, 'other']), rnd.choice(['src', 'EMBL', '.']), rnd.choice(types), str(s), str(e), rnd.choice(['.', '5.5', '0', '12']), rnd.choice(['+', '-', '.']), rnd.choice(['.', '0', '1', '2'])]
            if attrs or rnd.random() < 0.5: fields.append(';'.join(attrs) if rnd.random() < 0.9 else '; '.join(attrs))
            lines.append('\t'.join(fields))
        if rnd.random() < 0.2: lines.append('###')
        if rnd.random() < 0.1: lines += ['##FASTA', '>seq%d' % t, 'ACGT']
        out += '\n'.join(lines) + '\n'
    return out


# ---------------------------------------------------------------------------------------------------------------- swiss / pir / dasgff


def swiss_file(rnd):
    return ''.join(swiss_entry(rnd) for _ in range(rnd.choice([1, 1, 1, 2, 3])))


def swiss_entry(rnd):
    out = 'ID   %s_HUMAN   Reviewed;   %d AA.\nAC   P%05d;\n' % (rnd.choice(['ABC', 'XYZ', 'KIN1']), rnd.randint(50, 500), rnd.randint(1, 99999))
    for _ in range(rnd.choice([1, 2, 4, 7])):
        key = rnd.choice(TYPES['swiss'] + ['WEIRD', 'DOMAIN', 'CHAIN'])
        a = rnd.choice(['%d' % rnd.randint(1, 200), '<%d' % rnd.randint(1, 50), '?', '?%d' % rnd.randint(1, 9), '>%d' % rnd.randint(1, 50)])
        b = rnd.choice(['%d' % rnd.randint(1, 300), '>%d' % rnd.randint(1, 300), '?', '<%d' % rnd.randint(1, 9), '?%d' % rnd.randint(1, 9)])
        desc = rnd.choice(['', 'Kinase domain.', 'By similarity.', 'Zinc finger (By similarity).', 'N-acetylalanine.', 'Phosphoserine (by PKA) /FTId=PRO_0000012345.', 'x /FTId=VAR_000001.', 'a b c d e f g h'])
        out += 'FT   %-8s %6s %6s       %s\n' % (key, a, b, desc)
        if rnd.random() < 0.3: out += 'FT                                %s\n' % rnd.choice(['continued here.', '(Ref.2)', 'second /FTId=X_1.'])
    return out + ('//\n' if rnd.random() < 0.8 else '') + (rnd.choice(['', 'ID   OTHER_HUMAN   Reviewed;   9 AA.\nFT   DOMAIN 1 9 x.\n//\n']))


def pir_file(rnd):
    out = ''
    for _ in range(rnd.choice([1, 2, 4, 6])):
        ty = rnd.choice(TYPES['pir'] + ['weird', 'Domain', 'Disulfide bonds', 'Modified site'])
        loc = rnd.choice(['%d-%d' % (rnd.randint(1, 50), rnd.randint(51, 200)), '%d' % rnd.randint(1, 100), '%d,%d' % (rnd.randint(1, 20), rnd.randint(21, 40)), '%d-%d,%d-%d' % (1, 5, 8, 12), 'x-y', '1-', '-5', ',3'])
        note = rnd.choice(['', ' a note', ' some  extra   spaces', ' note #status predicted', ' n #status experimental #label L1'])
        out += 'F;%s/%s:%s\n' % (loc, ty.replace('_', ' '), note)
        if rnd.random() < 0.15: out += 'C;Comment line\n'
    return out


def das_file(rnd):
    def feat(i):
        o = '<FEATURE id="f%d" label="lab%d">\n' % (i, i)
        if rnd.random() < 0.9: o += ' <TYPE id="%s" category="%s" cvId="SO:0000%03d">%s</TYPE>\n' % (rnd.choice(['gene', 'exon', 'SO:0000704', 'CDS', 'weird']), rnd.choice(['transcription', 'other']), rnd.randint(1, 999), 'txt')
        if rnd.random() < 0.9: o += ' <METHOD id="%s"%s>m</METHOD>\n' % (rnd.choice(['prog', 'EMBL']), rnd.choice(['', ' cvId="ECO:0000001"']))
        if rnd.random() < 0.9: o += ' <START>%d</START>\n <END>%d</END>\n' % (rnd.randint(1, 100), rnd.randint(101, 300))
        if rnd.random() < 0.5: o += ' <SCORE>%s</SCORE>\n' % rnd.choice(['5.5', '-', '12', '0'])
        if rnd.random() < 0.5: o += ' <ORIENTATION>%s</ORIENTATION>\n' % rnd.choice(['+', '-', '0'])
        if rnd.random() < 0.4: o += ' <PHASE>%s</PHASE>\n' % rnd.choice(['0', '1', '2', '-'])
        for _ in range(rnd.choice([0, 0, 1, 2])): o += ' <NOTE>%s</NOTE>\n' % rnd.choice(['a note', 'text, with comma', ''])
        if rnd.random() < 0.2: o += ' <LINK href="http://x.org/%d">link</LINK>\n' % i
        if rnd.random() < 0.2: o += ' <TARGET id="t%d" start="%d" stop="%d">t</TARGET>\n' % (i, rnd.randint(1, 9), rnd.randint(10, 20))
        if rnd.random() < 0.15: o += ' <GROUP id="g%d"%s%s/>\n' % (i, rnd.choice(['', ' type="gene"']), rnd.choice(['', ' label="L"']))
        if rnd.random() < 0.15: o += ' <PARENT id="p%d"/>\n' % i
        if rnd.random() < 0.15: o += ' <PART id="q%d"/>\n' % i
        return o + '</FEATURE>\n'
    seg = '<SEGMENT id="%s" start="%d" stop="%d" version="1.0" label="x">\n' % (rnd.choice(['seq1', 'chr2', '']), rnd.randint(1, 5), rnd.randint(100, 500))
    body = ''.join(feat(i) for i in range(rnd.choice([1, 2, 3])))
    return '<?xml version="1.0" standalone="no"?>\n<DASGFF>\n<GFF version="1.0" href="http://x">\n' + seg + body + '</SEGMENT>\n</GFF>\n</DASGFF>\n'


def damaged(rnd, text):
    """a damaged or mixed version of a file"""
    lines = text.split('\n')
    r = rnd.random()
    if r < 0.25 and len(lines) > 2: del lines[rnd.randrange(len(lines))]
    elif r < 0.5: lines.insert(rnd.randrange(len(lines) + 1), rnd.choice(['', 'random text', '#comment', '//', 'ID   X; 1 BP.', '     key 1..2', '##gff-version 3', 'F;1-2/Domain:', '\t\t']))
    elif r < 0.7 and lines: i = rnd.randrange(len(lines)); lines[i] = lines[i][:rnd.randint(0, len(lines[i]))]
    elif r < 0.85: lines = [l.replace('\t', '  ') for l in lines]
    else: lines = [l + '\r' for l in lines]
    return '\n'.join(lines)


# ---------------------------------------------------------------------------------------------------------------- driver


def run(d, text, ufo, extra, timeout=4):
    for nme in ('out.dbg',):
        if os.path.exists(os.path.join(d, nme)): os.remove(os.path.join(d, nme))
    open(os.path.join(d, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
    hang = False
    try:
        p = subprocess.run(['featcopy', '-features', ufo] + extra + ['-outfeat', 'out.dbg', '-offormat', 'debug'], cwd=d, capture_output=True, timeout=timeout)
        rc, err = p.returncode, p.stderr.decode('latin-1')
    except subprocess.TimeoutExpired:
        rc, err, hang = -9, '', True
        subprocess.run(['pkill', '-9', 'featcopy'])
    q = os.path.join(d, 'out.dbg')
    dbg = open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
    if err.startswith('Read and write a feature table\n'): err = err[len('Read and write a feature table\n'):]
    return rc, dbg, err, hang


KINDS = ['embl', 'genbank', 'refseq', 'genpept', 'refseqp', 'gff2', 'gff3', 'gff3emboss', 'swiss', 'pir', 'dasgff']
FMTNAME = {'embl': 'embl', 'genbank': 'genbank', 'refseq': 'refseq', 'genpept': 'genpept', 'refseqp': 'refseqp', 'gff2': 'gff2', 'gff3': 'gff3', 'gff3emboss': 'gff3emboss', 'swiss': 'swiss', 'pir': 'pir', 'dasgff': 'dasgff'}


def make(rnd, kind):
    if kind in ('embl', 'genbank'): return embl_file(rnd, kind)
    if kind in ('refseq', 'genpept', 'refseqp'): return five_blank_file(rnd, kind)
    if kind == 'gff2': return gff2_file(rnd)
    if kind == 'gff3': return gff3_file(rnd)
    if kind == 'gff3emboss': return gff3_file(rnd, old=True)
    if kind == 'swiss': return swiss_file(rnd)
    if kind == 'pir': return pir_file(rnd)
    return das_file(rnd)


def main(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        i = 0
        while len(rows) < n:
            kind = KINDS[i % len(KINDS)]; i += 1
            text = make(rnd, kind)
            if rnd.random() < 0.2: text = damaged(rnd, text)
            if rnd.random() < 0.04: text = rnd.choice(['', '\n\n', 'just some text\nmore text\n', '#only a comment\n', '//\n'])
            r = rnd.random()
            if r < 0.4: ufo = 'in.dat'; extra = []
            elif r < 0.5: ufo = rnd.choice(['in.dat[1:50]', 'in.dat[5:20:r]', 'in.dat[:]', '  in.dat', 'in.dat ']); extra = rnd.choice([[], ['-fformat1', 'embl'], ['-fformat1', 'bogus']])
            elif r < 0.8: ufo = '%s::in.dat' % rnd.choice([FMTNAME[kind], FMTNAME[kind], rnd.choice(['gff3', 'embl', 'swiss', 'pir', 'gff2', 'genpept', 'refseq', 'refseqp', 'dasgff', 'GFF3', 'EMBL', 'bogus']), 'unknown']); extra = []
            else:
                ufo = 'in.dat'; extra = ['-fbegin1', str(rnd.randint(1, 120)), '-fend1', str(rnd.randint(60, 360))]
                if rnd.random() < 0.3: extra = ['-fbegin1', str(rnd.randint(-50, 80))] + (['-fend1', str(rnd.randint(-50, 300))] if rnd.random() < 0.5 else [])
            rc, dbg, err, hang = run(d, text, ufo, extra)
            rows.append({'id': len(rows), 'kind': kind, 'text': text, 'ufo': ufo, 'extra': extra, 'rc': rc, 'debug': dbg, 'err': err, 'hang': hang})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', Counter((r['kind'], r['rc'] == 0 and bool(r['debug'])) for r in rows))




# ---------------------------------------------------------------------------------------------------------------- the features of a sequence read from a separate file (-ufo, -fformat, -fopenfile)


def run_seqret(d, files, args, timeout=6):
    for nme in os.listdir(d):
        os.remove(os.path.join(d, nme))
    for nme, text in files.items():
        open(os.path.join(d, nme), 'w', newline='\n', encoding='latin-1').write(text)
    hang = False
    try:
        p = subprocess.run(['seqret'] + args + ['-feature', '-outseq', 'o.fa', '-offormat2', 'debug', '-auto'], cwd=d, capture_output=True, timeout=timeout)
        rc, err = p.returncode, p.stderr.decode('latin-1')
    except subprocess.TimeoutExpired:
        rc, err, hang = -9, '', True
        subprocess.run(['pkill', '-9', 'seqret'])
    q = os.path.join(d, 'o.debug')
    dbg = open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
    if err.startswith('Read and write (return) sequences\n'): err = err[len('Read and write (return) sequences\n'):]
    return rc, dbg, err, hang


def main_ufo(out, n, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        i = 0
        while len(rows) < n:
            kind = KINDS[i % len(KINDS)]; i += 1
            name = rnd.choice(['SEQA', 'seq1', 'AB123456'])
            seq = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(60, 300)))
            text = make(rnd, kind)
            if rnd.random() < 0.15: text = damaged(rnd, text)
            files = {'a.fa': '>%s a test sequence\n%s\n' % (name, seq), 'f.dat': text}
            r = rnd.random()
            args = ['a.fa']
            mode = 'ufo'
            if r < 0.3: args += ['-ufo1', 'f.dat']
            elif r < 0.55: args += ['-ufo1', '%s::f.dat' % rnd.choice([FMTNAME[kind], FMTNAME[kind], 'gff3', 'embl', 'swiss', 'pir', 'GFF3', 'bogus', 'gff', 'em', 'sw', 'nbrf', 'unknown'])]
            elif r < 0.7: args += ['-ufo1', 'f.dat', '-fformat1', rnd.choice([FMTNAME[kind], 'gff3', 'bogus', 'embl'])]
            elif r < 0.8: args += ['-fformat1', FMTNAME[kind]]; files['%s.%s' % (name, FMTNAME[kind])] = text; mode = 'name'
            elif r < 0.9: args += ['-fformat1', rnd.choice([FMTNAME[kind], 'pir']), '-fopenfile1', 'f.dat']
            else: args += ['-ufo1', rnd.choice(['missing.dat', 'x::f.dat', ':f.dat', 'f.dat:', 'a:b:c', '']) ]
            rc, dbg, err, hang = run_seqret(d, files, args)
            rows.append({'id': len(rows), 'kind': kind, 'name': name, 'seq': seq, 'files': files, 'args': args[1:], 'rc': rc, 'debug': dbg, 'err': err, 'hang': hang})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', Counter((r['kind'], bool(r['debug'])) for r in rows))


def main_write(out, n, seed):
    """featcopy -features FILE -outfeat out.X -offormat F: the feature tables written in an output format (dates faked with scripts/alignment/faketime.c)"""
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        work = os.path.join(d, 'w'); os.mkdir(work)
        formats = ['gff3', 'gff3', 'gff', 'unknown', 'gff2', 'gff2', 'embl', 'embl', 'em', 'genbank', 'genbank', 'gb', 'ddbj', 'GFF3', 'EMBL']
        i = 0
        while len(rows) < n:
            kind = KINDS[i % len(KINDS)]; i += 1
            text = make(rnd, kind)
            if rnd.random() < 0.15: text = damaged(rnd, text)
            fmt = rnd.choice(formats)
            ufo = 'in.dat'
            extra = []
            if rnd.random() < 0.25: extra = ['-fbegin1', str(rnd.randint(1, 100)), '-fend1', str(rnd.randint(60, 300))]
            for nme in os.listdir(work): os.remove(os.path.join(work, nme))
            open(os.path.join(work, 'in.dat'), 'w', newline='\n', encoding='latin-1').write(text)
            epoch = rnd.randint(978307200, 2051222400)
            env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
            hang = False
            try:
                p = subprocess.run(['featcopy', '-features', ufo] + extra + ['-outfeat', 'out.dat', '-offormat', fmt], cwd=work, capture_output=True, timeout=4, env=env)
                rc, err = p.returncode, p.stderr.decode('latin-1')
            except subprocess.TimeoutExpired:
                rc, err, hang = -9, '', True
                subprocess.run(['pkill', '-9', 'featcopy'])
            q = os.path.join(work, 'out.dat')
            outtext = open(q, encoding='latin-1', newline='').read() if os.path.exists(q) and not hang else ''
            if err.startswith('Read and write a feature table\n'): err = err[len('Read and write a feature table\n'):]
            rows.append({'id': len(rows), 'kind': kind, 'text': text, 'ufo': ufo, 'extra': extra, 'format': fmt, 'epoch': epoch, 'rc': rc, 'out': outtext, 'err': err, 'hang': hang})
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', Counter((r['format'].lower(), bool(r['out'])) for r in rows))


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == '--write':
    main_write(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    sys.exit(0)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == '--ufo':
    main_ufo(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    sys.exit(0)


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1000, int(sys.argv[3]) if len(sys.argv) > 3 else 20261401)
