#!/usr/bin/env python3
"""B09 WP63 (F282) oracle of the EMBOSS 6.6.0 programs of work package 63 (report layer over the stream runner: tcode, ...):
    <program> -auto ARGS ...      every row: arguments, input files, standard output, error stream, exit status, written files (see emboss_prog_oracle_lib.py)
Usage: python3 scripts/alignment/generate_emboss_prog_wp63_oracle.py PROGRAM OUT.jsonl.gz [n=300] [seed]      (needs EMBOSS 6.6.0)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402
import generate_seqret_oracle as S  # noqa: E402
import generate_emboss_prog_seq_oracle as G  # noqa: E402

REPORT_FORMATS = ['simple', 'embl', 'genbank', 'gff', 'gff3', 'pir', 'swiss', 'dasgff', 'debug', 'listfile', 'dbmotif', 'diffseq', 'draw', 'restrict', 'excel', 'feattable', 'motif', 'nametable', 'regions', 'seqtable', 'srs', 'summary', 'table', 'tagseq']


def report_quals(rnd, p=0.5, master=None, formats=None):
    """the report qualifiers (-rformat -raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq -rname -rextension -rdirectory)"""
    a = []
    if rnd.random() < p:
        a += ['-rformat', rnd.choice(formats or REPORT_FORMATS)]
    for o in ['-raccshow', '-rdesshow', '-rusashow']:
        if rnd.random() < 0.12: a.append(o)
    for o in ['-norscoreshow', '-norstrandshow']:
        if rnd.random() < 0.1: a.append(o)
    if rnd.random() < 0.1: a += ['-rmaxall', str(rnd.choice([0, 1, 2, 5, 10, 100, -1]))]
    if rnd.random() < 0.1: a += ['-rmaxseq', str(rnd.choice([0, 1, 2, 3, 50, -2]))]
    if rnd.random() < 0.03: a += [rnd.choice(['-rmaxall', '-rmaxseq']), rnd.choice(['x', '1.5', ''])]
    if rnd.random() < 0.03: a += ['-rformat', rnd.choice(['foo', '', 'TABLE', 'Excel', 'tab', 'gff2'])]
    return a


def report_file(rnd, ext, p_default=0.12):
    """-outfile and the name qualifiers of the default file"""
    r = rnd.random()
    if r < 0.55: return ['-outfile', 'stdout']
    if r < 0.8: return ['-outfile', 'out.txt']
    if r < 0.8 + p_default / 2: return []
    a = []
    if rnd.random() < 0.5: a += ['-rname', rnd.choice(['rep', 'my report', 'a:b', 'x.y', ''])]
    if rnd.random() < 0.5: a += ['-rextension', rnd.choice(['txt', 'rep', '', '.x'])]
    if rnd.random() < 0.3: a += ['-rdirectory', rnd.choice(['nodir', 'sub', '', 'sub/'])]
    if rnd.random() < 0.3: a += ['-outfile', rnd.choice(['sub/o.rep', 'nodir/o.rep', 'o.rep'])]
    return a


def etcode_variant(rnd):
    """Etcode.dat with changed numbers, extra comment and blank lines"""
    lines = open('/usr/share/EMBOSS/data/Etcode.dat').read().split('\n')
    out = []
    for ln in lines:
        if ln and not ln.startswith('#') and rnd.random() < 0.3:
            parts = ln.split()
            if len(parts) > 1:
                ln = rnd.choice([' ', '  ', '\t']).join('%.2f' % max(0.0, float(x) + rnd.uniform(-0.1, 0.1)) for x in parts)
        if rnd.random() < 0.05: out.append('# extra comment')
        out.append(ln)
    return '\n'.join(out)


def graph_quals(rnd, p_dev=1.0, devices=('data', 'data', 'data', 'none', 'bedgraph', 'null', 'text', 'DATA', 'da', 'bed', 'nu', 'te', 'n', 'bog', 'x', 'p')):
    """-graph and the associated qualifiers (the devices that write text: data, bedgraph, none / null / text; and the names the library refuses)"""
    a = []
    if rnd.random() < p_dev: a += ['-graph', rnd.choice(devices)]
    if rnd.random() < 0.25: a += ['-gtitle', rnd.choice(['my title', 't', 'A longer title', ' x ', 'tt'])]
    if rnd.random() < 0.12: a += ['-gsubtitle', rnd.choice(['sub', 's'])]
    if rnd.random() < 0.2: a += ['-gxtitle', rnd.choice(['xx', 'X axis', 'x'])]
    if rnd.random() < 0.2: a += ['-gytitle', rnd.choice(['yy', 'Y axis', 'y'])]
    if rnd.random() < 0.25: a += ['-goutfile', rnd.choice(['gg', 'out', 'a.b', 'sub/g'])]
    if rnd.random() < 0.1: a += ['-gdirectory', rnd.choice(['sub', 'nodir', 'sub/'])]
    if rnd.random() < 0.15: a += ['-gdesc', rnd.choice(['A description', 'd', 'plot'])]
    if rnd.random() < 0.05: a += [rnd.choice(['-gtitle', '-gytitle', '-gxtitle', '-goutfile']), '']
    return a


def gen_tcode(rnd, i, d):
    n = rnd.choice([1, 1, 1, 2, 3, 5])
    first = rnd.choice([250, 300, 400, 600, 900, 1500, 230])
    items = []
    for k in range(n):
        ln = first if k == 0 else (first + rnd.randint(0, 500) if rnd.random() < 0.88 else rnd.randint(1, first))     # a later sequence shorter than the window kills the reference (an array index below 0)
        items += G.seqs_for(rnd, 'x', 'n', 1, (ln, ln))
    items = [(nme + ('' if k == 0 else str(k)), de, sq) for k, (nme, de, sq) in enumerate(items)]
    if rnd.random() < 0.05:
        items = G.seqs_for(rnd, 'x', 'n', rnd.choice([1, 2]), (1, 150))
    text, fmt = G.seq_file(rnd, d, items)
    files = {'in.dat': text}
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.06)
    if rnd.random() < 0.5:
        args += ['-window', str(rnd.choice([200, 200, 210, 250, 300, 400, 500, 100, 0, -5, 1000, 'x']))]
    if rnd.random() < 0.4:
        args += ['-step', str(rnd.choice([1, 2, 3, 3, 5, 10, 30, 100, 0, -2, 1000]))]
    if rnd.random() < 0.08:
        files['dat.txt'] = etcode_variant(rnd)
        args += ['-datafile', 'dat.txt']
    elif rnd.random() < 0.02:
        args += ['-datafile', rnd.choice(['Etcode.dat', 'nonexist.dat', ''])]
    plot = rnd.random() < 0.3
    if plot:
        args.append('-plot')
        args += graph_quals(rnd)
    elif rnd.random() < 0.08:
        args += graph_quals(rnd)
    if plot and rnd.random() < 0.3:
        pass                                                  # report to the default file or none (-outfile not given: no report)
    else:
        args += report_quals(rnd)
        args += report_file(rnd, 'tcode')
    if 'sub' in ' '.join(args):
        files['sub/keep'] = 'x'
    return args, files, ({'vg': True} if plot else {})


def stop_protein(rnd, n):
    """a protein with stops (*) at random places: consecutive stops, a stop at the end, lower case"""
    res = []
    while len(res) < n:
        if rnd.random() < 0.12: res.extend('*' * rnd.choice([1, 1, 1, 2, 3]))
        else: res.extend(rnd.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(rnd.randint(1, 25)))
    t = ''.join(res)[:n]
    if rnd.random() < 0.25: t = t.rstrip('*')
    elif rnd.random() < 0.3: t += '*'
    if rnd.random() < 0.15: t = t.lower()
    return t


FEAT_FORMATS = ['gff', 'gff2', 'gff3', 'bed', 'pir', 'swiss', 'embl', 'genbank', 'dasgff', 'draw', 'debug', 'refseqp', 'unknown', 'foo', 'GFF', 'nbrf', 'sw', 'em']
SEQ_OUT_FORMATS = ['fasta', 'fasta', 'embl', 'genbank', 'gcg', 'pir', 'swiss', 'raw', 'msf', 'clustal', 'phylip', 'nexus', 'ig', 'staden', 'strider', 'fastq', 'listfile', 'acedb', 'debug', 'foo', 'sam', 'das']


def gen_checktrans(rnd, i, d):
    n = rnd.choice([1, 1, 2, 2, 3])
    items = []
    for k in range(n):
        nme = S.rand_name(rnd) + ('' if k == 0 else str(k))
        items.append((nme, S.rand_desc(rnd), stop_protein(rnd, rnd.choice([rnd.randint(1, 20), rnd.randint(10, 120), rnd.randint(50, 300)]))))
    text, fmt = G.seq_file(rnd, d, items, rnd.choice(['fasta'] * 6 + ['gcg', 'embl', 'swiss', 'pir', 'raw', 'genbank', 'ncbi']))
    files = {'in.dat': text}
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.12)
    if rnd.random() < 0.7: args += ['-orfml', str(rnd.choice([1, 2, 3, 5, 5, 8, 10, 20, 50, 100, 0, -1, 'x']))]
    if rnd.random() < 0.2: args.append('-noaddlast')
    r = rnd.random()
    args += ['-outfile', rnd.choice(['stdout', 'stdout', 'out.txt'])] if r < 0.85 else ([] if r < 0.95 else ['-outfile', 'sub/o.txt'])
    r = rnd.random()
    if r < 0.75: args += ['-outseq', rnd.choice(['stdout', 'stdout', 's.out', 'embl::s.out', 'fasta::stdout', 'msf::stdout', 'raw::q.seq'])]
    elif r < 0.9: pass
    else: args += ['-outseq', rnd.choice(['sub/s.out', 'foo::s.out'])]
    if rnd.random() < 0.3: args += ['-osformat3', rnd.choice(SEQ_OUT_FORMATS)]
    elif rnd.random() < 0.05: args += ['-osformat', rnd.choice(SEQ_OUT_FORMATS)]
    r = rnd.random()
    if r < 0.7: args += ['-outfeat', rnd.choice(['stdout', 'stdout', 'f.gff', 'f.out', 'pir::f.out', 'foo::f.out', 'gff2::stdout', 'embl::f.out'])]
    elif r < 0.9: pass
    else: args += ['-outfeat', rnd.choice(['sub/f.gff', 'nodir/f.gff', ''])]
    if rnd.random() < 0.35: args += ['-offormat4', rnd.choice(FEAT_FORMATS)]
    if rnd.random() < 0.08: args += ['-ofname4', rnd.choice(['nm', 'a b', 'x.y'])]
    if rnd.random() < 0.06: args += ['-ofdirectory4', rnd.choice(['sub', 'nodir', 'sub/'])]
    if rnd.random() < 0.04: args += ['-ofopenfile4', 'op.gff']
    if rnd.random() < 0.04: args += ['-ofextension4', 'ext']
    if rnd.random() < 0.05: args += ['-odirectory2', rnd.choice(['sub', 'nodir', 'sub/'])]
    if rnd.random() < 0.06: args += ['-osdirectory3', rnd.choice(['sub', 'nodir'])]
    if rnd.random() < 0.05: args += ['-ossingle3']
    if rnd.random() < 0.04: args += [rnd.choice(['-osname3', '-osextension3']), rnd.choice(['xx', 'y'])]
    if rnd.random() < 0.04: args += ['-outfile', 'nodir/o.txt']
    if 'sub' in ' '.join(args): files['sub/keep'] = 'x'
    return args, files


def gen_plotorf(rnd, i, d):
    lens = rnd.choice([(1, 40), (30, 300), (100, 900), (200, 2000)])
    items = G.seqs_for(rnd, 'x', 'n', rnd.choice([1, 1, 1, 2]), lens)
    text, fmt = G.seq_file(rnd, d, items)
    files = {'in.dat': text}
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.15)
    if rnd.random() < 0.4: args += ['-start', rnd.choice(['ATG', 'atg', 'ATG,GTG', 'GTG, TTG', 'ATG,CTG,TTG,GTG', 'AT', 'ATGC', '', ',', 'ATG,,GTG', ',ATG', 'A,T,G,C', 'NNN', 'ATG\nGTG'])]
    if rnd.random() < 0.4: args += ['-stop', rnd.choice(['TAA,TAG,TGA', 'taa', 'TAA', 'TAG, TGA', 'TA', '', 'TAAG,TGA', 'TAA,,TAG', 'TGA,', 'ATG', 'XXX'])]
    args += ['-graph', rnd.choice(['data', 'data', 'data', 'data', 'none', 'null', 'bedgraph', 'DATA', 'da', 'bed', 'te', 'bog', 'text'])]
    args += graph_quals(rnd, p_dev=0.0)
    if 'sub' in ' '.join(args): files['sub/keep'] = 'x'
    return args, files, {'vg': True}


AFORMATS = ['simple', 'simple', 'simple', 'pair', 'fasta', 'a2m', 'msf', 'clustal', 'mega', 'meganon', 'nexus', 'nexusnon', 'phylip', 'phylipnon', 'selex', 'treecon', 'markx0', 'markx1', 'markx2', 'markx3', 'markx10', 'match', 'multiple', 'score', 'srs', 'tcoffee', 'unknown', 'sam', 'foo', 'bam']
MATRICES_N = ['EDNAFULL', 'EDNAMAT', 'EDNASIMPLE', 'ENUC.4.4']
MATRICES_P = ['EBLOSUM62', 'EBLOSUM30', 'EBLOSUM90', 'EPAM250', 'EPAM10', 'EBLOSUM62-12', 'EBLOSUMN']


def mutate(rnd, s, alphabet, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate: out.append(rnd.choice(alphabet))
        elif r < rate * 1.3: continue
        elif r < rate * 1.6: out.append(c + rnd.choice(alphabet))
        else: out.append(c)
    return ''.join(out)


def gen_wordfinder(rnd, i, d):
    kind = rnd.choice(['n', 'n', 'n', 'p'])
    alphabet = 'ACGT' if kind == 'n' else 'ACDEFGHIKLMNPQRSTVWY'
    ntrg = rnd.choice([1, 1, 2, 3])
    nq = rnd.choice([1, 1, 2, 3])
    tl = rnd.choice([(40, 200), (100, 500), (30, 80)])
    trgs = G.seqs_for(rnd, 'x', kind, ntrg, tl)
    queries = []
    for k in range(nq):
        src = rnd.choice(trgs)[2].upper() if rnd.random() < 0.85 else ''.join(rnd.choice(alphabet) for _ in range(rnd.randint(30, 150)))
        a = rnd.randint(0, max(0, len(src) - 20)); b = rnd.randint(a + 15, len(src)) if len(src) > a + 15 else len(src)
        frag = mutate(rnd, src[a:b], alphabet, rnd.choice([0, 0.02, 0.08, 0.2]))
        flank = lambda: ''.join(rnd.choice(alphabet) for _ in range(rnd.randint(0, 30)))
        queries.append((S.rand_name(rnd) + 'q' + str(k), S.rand_desc(rnd), flank() + frag + flank()))
    tt, tfmt = G.seq_file(rnd, d, trgs)
    qt, qfmt = G.seq_file(rnd, d, queries, rnd.choice(['fasta'] * 5 + ['gcg', 'embl', 'pir', 'raw', 'genbank', 'swiss'] if kind == 'p' else ['fasta'] * 6 + ['gcg', 'embl', 'pir', 'raw', 'genbank']))
    files = {'a.dat': qt, 'b.dat': tt}
    args = ['-asequence', 'a.dat', '-bsequence', 'b.dat'] + G.seq_quals(rnd, qfmt, '1', 0.1) + G.seq_quals(rnd, tfmt, '2', 0.1)
    if rnd.random() < 0.5: args += ['-wordlen', str(rnd.choice([3, 4, 5, 6, 8, 10, 12, 2, 0, 20]))]
    if rnd.random() < 0.4: args += ['-width', str(rnd.choice([1, 2, 5, 8, 16, 16, 24, 40, 0]))]
    if rnd.random() < 0.4: args += ['-gapopen', str(rnd.choice([0, 1, 5, 10, 30, 50, -1, 2000, 7.5]))]
    if rnd.random() < 0.4: args += ['-gapextend', str(rnd.choice([0, 0.5, 1, 1.5, 3, 11, -2]))]
    if rnd.random() < 0.3: args += ['-datafile', rnd.choice(MATRICES_N if kind == 'n' else MATRICES_P + ['EDNAFULL'])]
    elif rnd.random() < 0.03: args += ['-datafile', rnd.choice(['nonexist', 'EBLOSUM62', 'BLOSUM62'])]
    for o in ['-limitmatch', '-lowmatch']:
        if rnd.random() < 0.2: args += [o, str(rnd.choice([0, 5, 10, 20, 40, 80, -3]))]
    for o in ['-limitalign', '-lowalign']:
        if rnd.random() < 0.2: args += [o, str(rnd.choice([0, 10, 20, 40, 80, 150, -3]))]
    if rnd.random() < 0.5: args += ['-aformat', rnd.choice(AFORMATS)]
    if rnd.random() < 0.15: args += ['-awidth', str(rnd.choice([10, 20, 30, 50, 100, 0, -4, 1000]))]
    for o in ['-aaccshow', '-adesshow', '-ausashow', '-aglobal']:
        if rnd.random() < 0.1: args.append(o)
    args += ['-outfile', rnd.choice(['stdout', 'stdout', 'out.txt'])]
    r = rnd.random()
    if r < 0.45: args += ['-errfile', rnd.choice(['stdout', 'err.txt'])]
    return args, files


def plant(rnd, s, site, k=1):
    """a site (IUPAC letters resolved at random; k substitutions) put at a random place of s"""
    iupac = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'W': 'AT', 'R': 'AG', 'Y': 'CT', 'N': 'ACGT'}
    t = ''.join(rnd.choice(iupac[c.upper()]) for c in site)
    t = list(t)
    for _ in range(rnd.choice([0, 0, 1, k])):
        t[rnd.randrange(len(t))] = rnd.choice('ACGT')
    t = ''.join(t)
    p = rnd.randint(0, max(0, len(s) - len(t)))
    return s[:p] + t + s[p + len(t):]


def gen_marscan(rnd, i, d):
    n = rnd.choice([1, 1, 2, 3])
    items = []
    for k in range(n):
        ln = rnd.choice([(60, 300), (200, 900), (400, 1500)])
        sq = G.seqs_for(rnd, 'x', 'n', 1, ln)[0]
        s = sq[2].upper() if rnd.random() < 0.8 else sq[2]
        for _ in range(rnd.choice([0, 1, 2, 3, 4, 6])):
            s = plant(rnd, s, rnd.choice(['AWWRTAANNWWGNNNC', 'GNNNCWWNNTTAYWWT', 'AATAAYAA', 'TTRTTATT', 'AATAAYAA', 'TTRTTATT']))
        items.append((sq[0] + ('' if k == 0 else str(k)), sq[1], s))
    text, fmt = G.seq_file(rnd, d, items)
    args = ['-sequence', 'in.dat'] + G.seq_quals(rnd, fmt, '', 0.12)
    args += report_quals(rnd, 0.5, formats=['gff', 'gff3', 'table', 'simple', 'excel', 'embl', 'genbank', 'seqtable', 'feattable', 'debug', 'tagseq', 'dbmotif', 'regions', 'motif', 'srs', 'listfile', 'draw', 'nametable', 'restrict', 'summary', 'dasgff', 'diffseq', 'pir'])
    args += report_file(rnd, 'marscan')
    files = {'in.dat': text}
    if 'sub' in ' '.join(args): files['sub/keep'] = 'x'
    return args, files


GENS = {'tcode': gen_tcode, 'checktrans': gen_checktrans, 'plotorf': gen_plotorf, 'wordfinder': gen_wordfinder, 'marscan': gen_marscan}

if __name__ == '__main__':
    prog = sys.argv.pop(1)
    L.main(prog, GENS[prog], 300)
