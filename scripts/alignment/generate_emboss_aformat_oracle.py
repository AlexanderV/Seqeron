#!/usr/bin/env python3
"""ALIGN-GLOBAL/LOCAL/MULTI audit A2-27 (B09 F120) oracle: every EMBOSS 6.6.0 `-aformat` alignment output format, byte for byte.

Seqeron API: EmbossAlignmentWriter.Write(IReadOnlyList<EmbossPairAlignment>, EmbossAlignFormat, EmbossAlignOptions).

Reference (run here): the installed EMBOSS 6.6.0 `needle` / `water` / `stretcher` (integer scores) with `-auto -stdout -aformat3 <format>` (+ -awidth3, -aglobal3, -aaccshow3, -adesshow3, -ausashow3, -datafile, -gapopen, -gapextend)
on random related DNA / IUPAC-DNA / protein pairs (and one sequence against several: several alignments in one file). The alignment (gapped rows, offsets) is read back from the tool's own
`fasta` and `markx10` output, so the oracle tests the *writer*, not the aligner (needle/water ties are tested by the aligner oracles).
Output TSV (tests/.../TestData/Alignment/emboss_aformat_oracle.tsv):
  id  program  format  args  fastaA  fastaB  matrix  alignments  date  file
  args = the command line printed in the header (space separated, file names a.fa / b.fa); alignments = row1:row2:off1:offend1:off2:offend2:score joined by '|'; date = ISO UTC of the run; file = output (\\n \\t \\\\ escaped).
Usage: python3 scripts/alignment/generate_emboss_aformat_oracle.py OUTFILE [casesPerFormat] [seed]      (needs EMBOSS needle/water and gcc on PATH)
"""
import datetime, os, random, re, subprocess, sys, tempfile

FORMATS = ['unknown', 'fasta', 'a2m', 'msf', 'clustal', 'mega', 'meganon', 'nexus', 'paup', 'nexusnon', 'paupnon', 'phylip', 'phylipnon', 'selex', 'treecon',
           'markx0', 'markx1', 'markx2', 'markx3', 'markx10', 'match', 'multiple', 'pair', 'simple', 'sam', 'score', 'srs', 'srspair', 'tcoffee']
PROT_MATRICES = ['EBLOSUM62', 'EBLOSUM62', 'EBLOSUM45', 'EBLOSUM80', 'EPAM250', 'EBLOSUM30']
DNA_MATRICES = ['EDNAFULL', 'EDNAFULL', 'EDNAMAT', 'EDNASIMPLE']
GAPS = [(10, 0.5), (10, 0.5), (8, 1), (12.5, 0.25), (15, 3), (6, 0.5), (10, 1.5)]


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')


def rand_name(rnd):
    n = rnd.choice(['s', 'seq', 'My_seq', 'sp', 'x.y', 'NC_', 'A'])
    n += ''.join(rnd.choice('abcdefghijklmnopqrstuvwxyz0123456789_.-') for _ in range(rnd.randint(0, rnd.choice([3, 8, 20, 40]))))
    return n


def rand_desc(rnd):
    if rnd.random() < 0.5:
        return ''
    return ' '.join(rnd.choice(['alpha', 'beta', 'gene', 'protein', 'hypothetical', 'kinase', 'from', 'E.', 'coli', 'x1']) for _ in range(rnd.randint(1, 5)))


def mutate(rnd, s, alphabet, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate: out.append(rnd.choice(alphabet))
        elif r < rate * 1.3: continue
        elif r < rate * 1.6: out.append(c); out.append(rnd.choice(alphabet))
        else: out.append(c)
    return ''.join(out) or s[:1]


def make_pair(rnd, typ):
    alpha = 'ACGT' if typ == 'n' else 'ACDEFGHIKLMNPQRSTVWY'
    ext = alpha + ('RYSWKMBDHVN' if typ == 'n' and rnd.random() < 0.3 else ('BZX' if typ == 'p' and rnd.random() < 0.3 else ''))
    n = rnd.choice([12, 20, 35, 60, 90, 140, 220])
    base = ''.join(rnd.choice(ext) for _ in range(n))
    seqs = []
    for _ in range(rnd.choice([1, 1, 1, 2, 3])):
        s = mutate(rnd, base, alpha, rnd.choice([0.02, 0.1, 0.25, 0.4]))
        if rnd.random() < 0.3: s = s[rnd.randint(0, 5):]
        if rnd.random() < 0.3: s = s[:len(s) - rnd.randint(0, 5)] or s
        if rnd.random() < 0.15: s = ''.join(rnd.choice([c.lower(), c]) for c in s)
        seqs.append(s)
    a = base if rnd.random() < 0.5 else mutate(rnd, base, alpha, 0.05)
    if rnd.random() < 0.15: a = a.lower()
    return a, seqs


def fasta(name, desc, seq):
    return f'>{name}' + (f' {desc}' if desc else '') + '\n' + '\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60)) + '\n'


SHIM = None
EPOCH = 0


def run(prog, args, cwd):
    env = dict(os.environ, TZ='UTC', LD_PRELOAD=SHIM, FAKE_EPOCH=str(EPOCH))
    r = subprocess.run([prog] + args, cwd=cwd, env=env, capture_output=True)
    return subprocess.CompletedProcess(r.args, r.returncode, r.stdout.decode('latin-1'), r.stderr.decode('latin-1'))


def parse_aln(d, prog, base_args):
    """gapped rows and offsets of every alignment of the run, from the tool's fasta and markx10 outputs"""
    f = run(prog, base_args + ['-aformat3', 'fasta'], d).stdout
    m = run(prog, base_args + ['-aformat3', 'markx10'], d).stdout
    p = run(prog, base_args + ['-aformat3', 'pair'], d).stdout
    recs = re.findall(r'^>(\S*)[^\n]*\n((?:[^>\n]*\n?)*)', f, re.M)
    rows = [r[1].replace('\n', '') for r in recs]
    out = []
    infos = re.findall(r'>>>[^\n]*\n(?:;[^\n]*\n)*?>>#\d+\n; sw_score: (\S+)\n(?:;[^\n]*\n)*?>[^\n]*\n; sq_len: (\d+)\n; sq_type: .\n; al_start: (\d+)\n; al_stop: (\d+)\n(?:[^>\n][^\n]*\n)*>[^\n]*\n; sq_len: (\d+)\n; sq_type: .\n; al_start: (\d+)\n; al_stop: (\d+)\n', m)
    scores = re.findall(r'^# Score: (\S+)$', p, re.M)
    if len(infos) != len(rows) // 2 or len(scores) != len(infos):
        return None
    for i, (sc, l1, s1, e1, l2, s2, e2) in enumerate(infos):
        out.append(':'.join([rows[2 * i], rows[2 * i + 1], str(int(s1) - 1), str(int(l1) - int(e1)), str(int(s2) - 1), str(int(l2) - int(e2)), scores[i]]))
    return '|'.join(out)


def main(outfile, per, seed):
    global SHIM, EPOCH
    rnd = random.Random(seed)
    rows = []
    n = 0
    with tempfile.TemporaryDirectory() as d:
        # the run date is faked (LD_PRELOAD time() shim, compiled from faketime.c) so that every day / month / hour layout is tested
        SHIM = os.path.join(d, 'faketime.so')
        subprocess.run(['gcc', '-shared', '-fPIC', '-o', SHIM, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'faketime.c')], check=True)
        for case in range(per):
            EPOCH = rnd.randint(978307200, 2051222400)
            typ = 'n' if case % 2 == 0 else 'p'
            prog = ['needle', 'needle', 'needle', 'water', 'water', 'stretcher'][case % 6]
            a, bs = make_pair(rnd, typ)
            na, da = rand_name(rnd), rand_desc(rnd)
            if prog == 'stretcher': da = ''     # stretcher does not carry the description
            fa = fasta(na, da, a)
            fb = ''.join(fasta(rand_name(rnd) + str(i), '' if prog == 'stretcher' else rand_desc(rnd), b) for i, b in enumerate(bs))
            open(os.path.join(d, 'a.fa'), 'w', encoding='latin-1').write(fa)
            open(os.path.join(d, 'b.fa'), 'w', encoding='latin-1').write(fb)
            matrix = rnd.choice(DNA_MATRICES if typ == 'n' else PROT_MATRICES)
            go, ge = rnd.choice(GAPS)
            if prog == 'stretcher':
                go, ge = int(go), int(ge)
            opts = ['-datafile', matrix, '-gapopen', str(go), '-gapextend', str(ge)] if case % 3 else []
            if matrix == 'EDNAMAT' or matrix == 'EDNASIMPLE' or prog == 'stretcher':
                opts = ['-datafile', matrix, '-gapopen', str(go), '-gapextend', str(ge)]
            fmt_opts = []
            if rnd.random() < 0.4: fmt_opts += ['-awidth3', str(rnd.choice([10, 20, 30, 40, 60, 80, 100]))]
            if rnd.random() < 0.25: fmt_opts += [rnd.choice(['-aglobal3', '-noaglobal3'])]
            if rnd.random() < 0.2: fmt_opts += ['-adesshow3']
            if rnd.random() < 0.15: fmt_opts += ['-aaccshow3']
            if rnd.random() < 0.1: fmt_opts += ['-ausashow3']
            base_args = ['-asequence', 'a.fa', '-bsequence', 'b.fa', '-auto', '-stdout'] + opts + fmt_opts
            aln = parse_aln(d, prog, base_args)
            if aln is None:
                continue
            mx = matrix if opts else ('EDNAFULL' if typ == 'n' else 'EBLOSUM62')
            for fmt in FORMATS:
                args = base_args + ['-aformat3', fmt]
                res = run(prog, args, d)
                if res.returncode != 0 or not res.stdout:
                    continue
                out = res.stdout
                if fmt == 'sam' and prog == 'water' and '-aglobal3' in args and ':' in aln and any(int(x) > 0 for x in aln.split('|')[0].split(':')[2:5:2]):
                    continue    # alignCigar reads outside the row there (undefined behaviour in the tool)
                if re.search(r'[^\x00-\x7f]', out):
                    continue    # uninitialised memory printed by the tool (not reproducible)
                date = datetime.datetime.fromtimestamp(EPOCH, datetime.timezone.utc).replace(tzinfo=None)
                rows.append('\t'.join([str(n), prog, fmt, ' '.join(args), esc(fa), esc(fb), mx, aln, date.strftime('%Y-%m-%dT%H:%M:%S'), esc(out)]))
                n += 1
    with open(outfile, 'w', encoding='utf-8') as f:
        f.write('id\tprogram\tformat\targs\tfastaA\tfastaB\tmatrix\talignments\tdate\tfile\n')
        f.write('\n'.join(rows) + '\n')
    print(len(rows), 'rows')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 130, int(sys.argv[3]) if len(sys.argv) > 3 else 20261119)
