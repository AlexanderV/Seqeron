#!/usr/bin/env python3
"""A5-8c (B09 F192) oracle: the report FILES of EMBOSS 6.6.0 fuzznuc, fuzzpro and fuzztran for (a) FASTA headers with NCBI-style identifiers (ajSeqParseNcbi: gb|AB123456.1|locus desc, sp|P12345|NAME_HUMAN, gi|123|gb|..|,
gnl|db|id, gnl|BL_ORD_ID|n, lcl|, bbs|, pir||, prf||, pdb|, pat|, unknown prefixes with 2 - 5 barred tokens; ids with ':' ',' '/' (seqSetName)) and (b) pattern FILES (`-pattern @file`: one pattern per line or FASTA style with
`<mismatch=n>`, `-pname`, `-pformat simple|fasta`, `-pmismatch`) in every report format. The run date is faked (scripts/alignment/faketime.c).
JSONL.GZ rows as generate_emboss_report_formats_oracle.py plus: kind ('ncbi' | 'patfile' | 'slow'), patfile (the pattern file text or null), pargs (the pattern qualifiers as typed, e.g. '-pname x -pformat simple').
Usage: python3 scripts/alignment/generate_emboss_ncbi_pattern_oracle.py OUTFILE.jsonl.gz [per_program_format=8] [seed]      (needs the EMBOSS 6.6.0 programs and gcc; the 'slow' rows take ~30 s each: the program loops 2^32 times for a header
with exactly two barred tokens)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_emboss_report_formats_oracle as G  # noqa: E402
from generate_emboss_gap_oracle import NUC, AA  # noqa: E402

FORMATS = [f for f in G.FORMATS if f != 'diffseq']
DESC = ['', '', 'Homo sapiens', 'putative kinase', 'OS=Mus musculus GN=abc', 'unknown protein', 'x', '12 ab', 'long description for a test sequence']
ACCS = ['AB123456.1', 'AB123456', 'P12345', 'Q9Y6K9', 'NM_000001.3', 'X56789.2', 'NP_001234.1', 'BK000001', 'A0A024R161', '1ABC']


def ncbi_id(rnd, i):
    acc = rnd.choice(ACCS); loc = rnd.choice(['locus%d' % i, 'NAME_HUMAN', 'a_b', 'abc:def', 'x,y', 'p/q', 'LOC%d' % (i * 7), 'EF-Tu'])
    form = rnd.choice(['gb', 'emb', 'dbj', 'sp', 'ref', 'tpg', 'tpd', 'gi', 'gi1', 'gnl', 'blord', 'lcl', 'bbs', 'pir', 'prf', 'pdb', 'pat', 'unk3', 'unk4', 'unk5', 'gb2', 'trail', 'gi_gb', 'gi_sp'])
    n = rnd.randint(1, 99999)
    return {'gb': 'gb|%s|%s' % (acc, loc), 'emb': 'emb|%s|%s' % (acc, loc), 'dbj': 'dbj|%s|%s' % (acc, loc), 'sp': 'sp|%s|%s' % (acc, loc), 'ref': 'ref|%s|%s' % (acc, loc),
            'tpg': 'tpg|%s|%s' % (acc, loc), 'tpd': 'tpd|%s|' % acc, 'gi': 'gi|%d' % n, 'gi1': 'gi|%d|gb|%s|%s' % (n, acc, loc), 'gnl': 'gnl|%s|%s' % (rnd.choice(['db', 'taxon', 'LOCAL']), loc),
            'blord': 'gnl|BL_ORD_ID|%d' % i, 'lcl': 'lcl|%s' % loc, 'bbs': 'bbs|%d' % n, 'pir': 'pir||%s' % loc, 'prf': 'prf||%s' % loc, 'pdb': 'pdb|1ABC|%s' % rnd.choice('ABH'),
            'pat': 'pat|US|%d' % n, 'unk3': 'xx|%s|%s' % (acc, loc), 'unk4': 'q|xx|%s|%s' % (acc, loc), 'unk5': 'a|b|c|%s|%s' % (acc, loc), 'gb2': 'gb|%s' % acc, 'trail': 'gb|%s|' % loc,
            'gi_gb': 'gi|%d|gb|%s|' % (n, acc), 'gi_sp': 'gi|%d|sp|%s|%s' % (n, acc, loc)}[form], form


def ncbi_fasta(rnd, make_seq, n_seq):
    out = []
    for i in range(n_seq):
        hid, form = ncbi_id(rnd, i + 1)
        if form == 'blord':
            rest = ncbi_id(rnd, i + 1)[0]
            desc = rest + ((' ' + rnd.choice(DESC[2:])) if rnd.random() < 0.6 else '')
        else:
            desc = rnd.choice(DESC)
        if rnd.random() < 0.15: hid = rnd.choice(['plain%d' % i, 'a:b%d' % i, 'x,y%d' % i, 'p/q%d' % i, 'one\\two%d' % i])
        out.append('>' + hid + (' ' + desc if desc else '') + '\n' + make_seq() + '\n')
    return ''.join(out)


def patterns_nuc(rnd, base):
    out = []
    for _ in range(rnd.randint(2, 4)):
        st = rnd.randint(0, len(base) - 8)
        out.append(rnd.choice([base[st:st + rnd.randint(3, 6)], base[st:st + 4].replace(base[st], 'N', 1), base[st:st + 3] + '[AC]' + base[st + 4:st + 6], 'RYR', 'AC(1,2)N', base[st:st + 2] + 'N(1,2)' + base[st + 3], base[st:st + 3]]))
    return out


def patterns_prot(rnd, base):
    out = []
    for _ in range(rnd.randint(2, 4)):
        st = rnd.randint(0, len(base) - 6)
        out.append(rnd.choice(['-'.join(base[st:st + 2]), '-'.join(base[st:st + 3]), base[st] + '-x-' + base[st + 2], '[' + base[st] + base[st + 1] + ']-' + base[st + 2], 'x(2)-' + base[st + 3], base[st] + '-x(1,2)-' + base[st + 3]]))
    return out


def patfile(rnd, pats, simple_ok=True, mm_ok=True):
    """returns (text, pformat or None)"""
    style = rnd.choice(['simple', 'fasta', 'fasta', 'fasta_mm'] if mm_ok else ['simple', 'fasta', 'fasta'])
    if style == 'simple' and simple_ok:
        return ''.join(p + '\n' for p in pats), rnd.choice([None, 'simple'])
    names = [rnd.choice(['p%d' % i, 'pat%d' % i, 'motif %d' % i, 'a_longer_pattern_name%d' % i, '', 'x desc']) for i in range(len(pats))]
    out = []
    for nm, p in zip(names, pats):
        mm = rnd.choice([0, 1, 2]) if style == 'fasta_mm' else None
        h = '>' + nm + ('<mismatch=%d>' % mm if mm is not None and rnd.random() < 0.8 else '')
        if rnd.random() < 0.3 and len(p) > 4:
            k = rnd.randint(2, len(p) - 2); out.append(h + '\n' + p[:k] + '\n' + p[k:] + '\n')
        else:
            out.append(h + '\n' + p + '\n')
    return ''.join(out), rnd.choice([None, 'fasta'])


def run_long(prog, args, fasta, epoch, shim, d, timeout):
    for n in ('o.rep', 'o.orig'):
        p = os.path.join(d, n)
        if os.path.exists(p): os.remove(p)
    with open(os.path.join(d, 's.fa'), 'w', encoding='latin-1', newline='\n') as f: f.write(fasta)
    env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
    try:
        p = subprocess.run([prog] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, '', None
    rep = open(os.path.join(d, 'o.rep'), encoding='latin-1', newline='').read() if os.path.exists(os.path.join(d, 'o.rep')) else ''
    return p, rep, None


def make_row(rnd, prog, fmt, kind, shim, d, rows, fixed=None):
    seqs = rnd.choice([1, 2, 3])
    patfile_text = None; pargs = []
    if prog == 'fuzznuc':
        mk = lambda: ''.join(rnd.choice(NUC) for _ in range(rnd.randint(40, 100)))
    elif prog == 'fuzzpro':
        mk = lambda: ''.join(rnd.choice(AA) for _ in range(rnd.randint(40, 100)))
    else:
        mk = lambda: ''.join(rnd.choice(NUC) for _ in range(rnd.randint(60, 120)))
    fa = fixed if fixed is not None else (ncbi_fasta(rnd, mk, seqs) if kind == 'ncbi' else G.rnd_fasta(rnd, mk, seqs))
    base = fa.split('\n')[1].upper()
    extra = []
    if kind == 'patfile':
        pats = patterns_prot(rnd, base) if prog == 'fuzzpro' else (patterns_nuc(rnd, base) if prog == 'fuzznuc' else rnd.sample(['L-x-S', 'S-S', 'G-R', 'x(2)-K', '[LV]-T', 'M', 'A-A', 'T-x(1,2)-G'], rnd.randint(2, 3)))
        patfile_text, pf = patfile(rnd, pats, mm_ok=prog == 'fuzznuc')
        if pf: pargs += ['-pformat', pf]
        if rnd.random() < 0.4: pargs += ['-pname', rnd.choice(['nm', 'motif', 'zz1'])]
        if prog == 'fuzznuc' and rnd.random() < 0.4: pargs += ['-pmismatch', str(rnd.choice([0, 1, 2]))]
        with open(os.path.join(d, 'p.txt'), 'w', newline='\n') as f: f.write(patfile_text)
        patarg = ['-pattern', '@p.txt']
    else:
        if prog == 'fuzznuc':
            st = rnd.randint(0, len(base) - 6); pat = rnd.choice([base[st:st + 4], base[st:st + 3], 'RYR'])
        elif prog == 'fuzzpro':
            st = rnd.randint(0, len(base) - 5); pat = '-'.join(base[st:st + 2])
        else:
            pat = rnd.choice(['L-x-S', 'S-S', 'G-R', 'x(2)-K'])
        patarg = ['-pattern', pat]
    if prog == 'fuzznuc' and rnd.random() < 0.4: extra += ['-complement', 'Y']
    if prog == 'fuzztran': extra += ['-frame', rnd.choice(['F', '6', '1', '2', '-1', 'R']), '-table', '0']
    rq = G.report_quals(rnd); rtokens = [t for q in rq for t in q.split(' ')]
    sq = []
    args = ['s.fa'] + patarg + pargs + extra + ['-rformat', fmt] + rtokens + ['-auto', '-outfile', 'o.rep']
    epoch = rnd.randint(978307200, 2051222400)
    p, rep, orig = run_long(prog, args, fa, epoch, shim, d, 200 if kind == 'slow' else 8)
    if p is None: return False
    err = p.stderr.strip()
    if p.returncode != 0 and fmt not in G.PROTEIN_ONLY and fmt not in G.NUCLEOTIDE_ONLY:
        return False

    def echo(name, value):
        if value is None: return name
        if name in ('-complement',) and value == 'Y': return name
        if value == '' or any(c in value for c in ' \t*?[]{}|!&^'): return name + ' "' + value + '"'
        return name + ' ' + value
    entries = ['[-sequence] s.fa']; i = 1
    while i < len(args):
        a = args[i]
        if a.startswith('-') and i + 1 < len(args) and (not args[i + 1].startswith('-') or args[i + 1].lstrip('-').isdigit()): entries.append(echo(a, args[i + 1])); i += 2
        else: entries.append(echo(a, None)); i += 1
    rows.append({'id': len(rows), 'kind': kind, 'program': prog, 'format': fmt, 'params': ' '.join(patarg[1:] + extra), 'patfile': patfile_text, 'pargs': ' '.join(pargs), 'sq': '', 'rq': ' '.join(rq), 'cmd': entries,
                 'epoch': epoch, 'fasta': fa, 'report': rep, 'orig': None, 'error': err or None, 'rc': p.returncode})
    return True


def main(outfile, per, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        for kind in ('ncbi', 'patfile'):
            for prog in ('fuzznuc', 'fuzzpro', 'fuzztran'):
                for fmt in FORMATS:
                    made = 0; tries = 0
                    while made < per and tries < per * 6:
                        tries += 1
                        if make_row(rnd, prog, fmt, kind, shim, d, rows): made += 1
        # headers with exactly two barred tokens: the program steps 2^32 times (~28 s) and leaves the sequence unnamed
        for prog in ('fuzzpro', 'fuzznuc'):
            for k in range(3):
                fa = '>p1 d\n%s\n>%s|%s some desc\n%s\n' % ('MKLAAGLVHHHHKLMKLAAGLV' if prog == 'fuzzpro' else 'ACGTACGTTTGACGTAAACCCGGGTTTACGTA',
                                                           rnd.choice(['xx', 'aa', 'ab']), rnd.choice(['yy', 'AB123456', 'q']), 'MKLAAGLVHHHHKLMKLAAGLV' if prog == 'fuzzpro' else 'ACGTACGTTTGACGTAAACCCGGGTTTACGTA')
                make_row(rnd, prog, rnd.choice(['seqtable', 'table', 'gff', 'simple']), 'slow', shim, d, rows, fixed=fa)
    with gzip.open(outfile, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 8, int(sys.argv[3]) if len(sys.argv) > 3 else 20261221)
