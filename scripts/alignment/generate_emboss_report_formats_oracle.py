#!/usr/bin/env python3
"""A5-8b (B09 F184) oracle: the report FILES of the EMBOSS 6.6.0 programs fuzznuc, fuzzpro, fuzztran, etandem and equicktandem in EVERY `-rformat` the programs accept
(simple embl genbank gff gff3 pir swiss dasgff debug listfile dbmotif draw restrict excel feattable motif nametable regions seqtable srs summary table tagseq; diffseq needs 7 report tags and
is an error for these programs), with the report options (-raccshow -rdesshow -rscoreshow -rstrandshow -rusashow -rmaxall -rmaxseq), several sequences per file (names / accessions /
descriptions of FASTA headers: sp|P12345|NAME desc, gb|AB1234.1|, X56789.2, plain), the sequence qualifiers -sbegin / -send / -sreverse, and the extra output file `-origfile` of etandem /
equicktandem. The run date is faked (LD_PRELOAD time() shim, scripts/alignment/faketime.c) so that the `# Rundate` / `#!Date` lines are compared byte for byte.

JSONL.GZ (tests/.../TestData/Alignment/emboss_report_formats_oracle.jsonl.gz), one object per line:
  id, program, format, params (as generate_emboss_report_oracle.py: fuzznuc 'pattern|mismatch|complement', fuzzpro 'pattern', fuzztran 'pattern|frame', etandem 'min,max,threshold,mismatch,uniform',
  equicktandem 'max,threshold'), sq (sequence qualifiers text, may be ''), rq (report qualifiers text WITHOUT -rformat, e.g. '-raccshow -rmaxall 3'), cmd (the command-line entries of the header),
  epoch (the faked time), fasta (the input file), report (the -outfile text, bytes decoded latin-1), orig (the -origfile text of etandem / equicktandem, else null), error (the stderr when the
  program refused the format, else null: then `report` is empty).
Usage: python3 scripts/alignment/generate_emboss_report_formats_oracle.py OUTFILE.jsonl.gz [per_program_format=100] [seed]      (needs the EMBOSS 6.6.0 programs and gcc)
"""
import gzip, json, os, random, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_emboss_gap_oracle import mutate, NUC, AA, tandem_seq  # noqa: E402

FORMATS = ['simple', 'embl', 'genbank', 'gff', 'gff3', 'pir', 'swiss', 'dasgff', 'debug', 'listfile', 'dbmotif', 'diffseq', 'draw', 'restrict', 'excel', 'feattable', 'motif', 'nametable', 'regions',
           'seqtable', 'srs', 'summary', 'table', 'tagseq']
NAMES = ['seq%d', 'X%05d.2', 'a_b-c%d', 'Q%05d', 'chr%d', 'AB%06d', 'NM_%06d.3']
WORDS = ['AB123456.1', 'P12345', 'alpha', 'beta', 'gene', 'protein', 'kinase', 'putative', 'unknown', 'Homo sapiens', 'x', '12', 'test-case', 'OS=Mus musculus', 'long description for a test sequence']


def rnd_fasta(rnd, make_seq, n_seq):
    out = []
    for i in range(n_seq):
        pat = rnd.choice(NAMES)
        name = pat % tuple(rnd.randint(1000, 99999) for _ in range(pat.count('%')))
        if rnd.random() < 0.15: name = 'seq%d' % (i + 1)
        desc = ' '.join(rnd.choice(WORDS) for _ in range(rnd.randint(0, 4))) if rnd.random() < 0.6 else ''
        out.append('>' + name + (' ' + desc if desc else '') + '\n' + make_seq() + '\n')
    return ''.join(out)


def fuzznuc_case(rnd, seqs):
    def mk():
        s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(25, 90)))
        return s.lower() if rnd.random() < 0.15 else s
    fa = rnd_fasta(rnd, mk, seqs)
    base = fa.split('\n')[1].upper()
    st = rnd.randint(0, len(base) - 6)
    pat = rnd.choice([base[st:st + rnd.randint(3, 6)], base[st:st + 4].replace(base[st], 'N', 1), base[st:st + 3] + '[AC]' + base[st + 4:st + 6], 'AC' + 'N' * rnd.randint(1, 2) + 'GT', 'RYR',
                      base[st:st + 3], base[st] + '(2,3)', 'AC(1,2)N', base[st:st + 2] + 'N(1,2)' + base[st + 3], 'N(18,24)' if rnd.random() < 0.3 else base[st:st + 4]])
    mm = rnd.choice([0, 0, 1, 1, 2]); comp = rnd.choice([0, 1])
    return '%s|%d|%d' % (pat, mm, comp), fa, ['-pattern', pat, '-pmismatch', str(mm), '-complement', 'Y' if comp else 'N']


def fuzzpro_case(rnd, seqs):
    fa = rnd_fasta(rnd, lambda: ''.join(rnd.choice(AA) for _ in range(rnd.randint(25, 90))), seqs)
    s = fa.split('\n')[1]
    st = rnd.randint(0, len(s) - 5)
    pat = rnd.choice(['-'.join(s[st:st + 2]), '-'.join(s[st:st + 3]), s[st] + '-x-' + s[st + 2], '[' + s[st] + s[st + 1] + ']-' + s[st + 2], 'x(2)-' + s[st + 3], s[st] + '-x(1,2)-' + s[st + 3], 'x(1,3)-' + s[st + 3],
                      'x(20,26)' if rnd.random() < 0.3 else '-'.join(s[st:st + 2])])
    return pat, fa, ['-pattern', pat]


def fuzztran_case(rnd, seqs):
    fa = rnd_fasta(rnd, lambda: ''.join(rnd.choice(NUC) for _ in range(rnd.randint(30, 100))), seqs)
    pat = rnd.choice(['L-x-S', 'S-S', 'G-R', 'x(2)-K', '[LV]-T', 'M', 'A-A', 'T-x(1,2)-G', 'x(8,12)-L'])
    fr = rnd.choice(['F', 'F', '6', 'R', '1', '2', '3', '-1', '-2', '-3'])
    return pat + '|' + fr, fa, ['-pattern', pat, '-frame', fr, '-table', '0']


def etandem_case(rnd, seqs):
    fa = '>' + rnd.choice(NAMES).replace('%d', '1').replace('%05d', '12345').replace('%04d', '1234').replace('%06d', '123456') + (' ' + rnd.choice(WORDS) if rnd.random() < 0.5 else '') + '\n' + tandem_seq(rnd) + '\n'
    mn = rnd.choice([2, 3]); mx = mn + rnd.randint(1, 5); th = rnd.choice([8, 12, 15, 20]); mm = rnd.choice([0, 1]); un = rnd.choice([0, 1])
    return '%d,%d,%d,%d,%d' % (mn, mx, th, mm, un), fa, ['-minrepeat', str(mn), '-maxrepeat', str(mx), '-threshold', str(th)] + (['-mismatch'] if mm else []) + (['-uniform'] if un else [])


def equick_case(rnd, seqs):
    fa = '>' + rnd.choice(NAMES).replace('%d', '1').replace('%05d', '12345').replace('%04d', '1234').replace('%06d', '123456') + (' ' + rnd.choice(WORDS) if rnd.random() < 0.5 else '') + '\n' + tandem_seq(rnd) + '\n'
    mx = rnd.choice([4, 6, 8]); th = rnd.choice([5, 8, 12, 20])
    return '%d,%d' % (mx, th), fa, ['-maxrepeat', str(mx), '-threshold', str(th)]


CASES = {'fuzznuc': fuzznuc_case, 'fuzzpro': fuzzpro_case, 'fuzztran': fuzztran_case, 'etandem': etandem_case, 'equicktandem': equick_case}
PROTEIN_ONLY = {'pir', 'swiss'}
NUCLEOTIDE_ONLY = {'embl', 'genbank'}


def report_quals(rnd):
    q = []
    if rnd.random() < 0.3: q.append('-raccshow')
    if rnd.random() < 0.3: q.append('-rdesshow')
    if rnd.random() < 0.2: q.append(rnd.choice(['-norscoreshow', '-rscoreshow']))
    if rnd.random() < 0.2: q.append(rnd.choice(['-norstrandshow', '-rstrandshow']))
    if rnd.random() < 0.2: q.append('-rusashow')
    if rnd.random() < 0.2: q.append('-rmaxall %d' % rnd.randint(1, 6))
    if rnd.random() < 0.2: q.append('-rmaxseq %d' % rnd.randint(1, 3))
    return q


def seq_quals(rnd, prog):
    if rnd.random() > 0.22: return []
    q = []
    if rnd.random() < 0.7: q.append('-sbegin %d' % rnd.randint(1, 12))
    if rnd.random() < 0.7: q.append('-send %d' % rnd.randint(25, 60))
    if prog != 'fuzzpro' and rnd.random() < 0.4: q.append('-sreverse')
    return q


def run(prog, args, fasta, epoch, shim, d):
    for n in ('o.rep', 'o.orig'):
        p = os.path.join(d, n)
        if os.path.exists(p): os.remove(p)
    with open(os.path.join(d, 's.fa'), 'w', encoding='latin-1', newline='\n') as f: f.write(fasta)
    env = dict(os.environ, LD_PRELOAD=shim, FAKE_EPOCH=str(epoch))
    try:
        p = subprocess.run([prog] + args, cwd=d, capture_output=True, errors='replace', env=env, timeout=8)
    except subprocess.TimeoutExpired:
        return None, '', None           # the program does not terminate in practice (loops over ~2^32 positions of a negative fuzztran coordinate): no oracle row
    rep = open(os.path.join(d, 'o.rep'), encoding='latin-1', newline='').read() if os.path.exists(os.path.join(d, 'o.rep')) else ''
    orig = open(os.path.join(d, 'o.orig'), encoding='latin-1', newline='').read() if os.path.exists(os.path.join(d, 'o.orig')) else None
    return p, rep, orig


def main(outfile, per, seed):
    rnd = random.Random(seed)
    rows = []
    with tempfile.TemporaryDirectory() as d:
        shim = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', shim, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        for prog, fn in CASES.items():
            for fmt in FORMATS:
                made = 0; tries = 0; timeouts = 0
                while made < per and tries < per * 6:
                    tries += 1
                    seqs = 1 if prog in ('etandem', 'equicktandem') else rnd.choice([1, 1, 2, 3])
                    params, fasta, pargs = fn(rnd, seqs)
                    sq = seq_quals(rnd, prog); rq = report_quals(rnd)
                    order = [('seq', ['-sequence' if False else 's.fa'])]
                    tokens = [t for q in sq for t in q.split(' ')]
                    rtokens = [t for q in rq for t in q.split(' ')]
                    args = ['s.fa'] + pargs + ['-rformat', fmt] + rtokens + tokens + ['-auto', '-outfile', 'o.rep'] + (['-origfile', 'o.orig'] if prog in ('etandem', 'equicktandem') else [])
                    epoch = rnd.randint(978307200, 2051222400)
                    p, rep, orig = run(prog, args, fasta, epoch, shim, d)
                    if p is None:
                        timeouts += 1
                        continue
                    err = p.stderr.strip()
                    # the entries of the header (ajacd.c acdArgsParse): the positional file as `[-sequence] name`, then every qualifier as typed; a boolean given Y (or its `no` form given N) prints without value,
                    # a value with one of  * ? [ ] { } | ! & ^  (or a blank) is quoted
                    def echo(name, value):
                        if value is None: return name
                        if name in ('-complement', '-mismatch', '-uniform') and value == 'Y': return name
                        if value == '' or any(c in value for c in ' \t*?[]{}|!&^'): return name + ' "' + value + '"'
                        return name + ' ' + value
                    entries = ['[-sequence] s.fa']; i = 1
                    while i < len(args):
                        a = args[i]
                        if a.startswith('-') and i + 1 < len(args) and (not args[i + 1].startswith('-') or args[i + 1].lstrip('-').isdigit()): entries.append(echo(a, args[i + 1])); i += 2
                        else: entries.append(echo(a, None)); i += 1
                    if p.returncode != 0 and not (fmt in PROTEIN_ONLY or fmt in NUCLEOTIDE_ONLY or fmt == 'diffseq'):
                        continue
                    if rep and fmt != 'diffseq' and p.returncode == 0:
                        for e in entries:
                            if ('#    ' + e + '\n') not in rep and not fmt in ('excel', 'gff', 'gff3', 'embl', 'genbank', 'feattable', 'pir', 'swiss', 'dasgff'):
                                print('HEADER MISMATCH', prog, fmt, e, file=sys.stderr)
                    rows.append({'id': len(rows), 'program': prog, 'format': fmt, 'params': params, 'sq': ' '.join(sq), 'rq': ' '.join(rq), 'cmd': entries, 'epoch': epoch, 'fasta': fasta,
                                 'report': rep, 'orig': orig, 'error': err or None, 'rc': p.returncode})
                    made += 1
    with gzip.open(outfile, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r, ensure_ascii=True) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 100, int(sys.argv[3]) if len(sys.argv) > 3 else 20261220)
