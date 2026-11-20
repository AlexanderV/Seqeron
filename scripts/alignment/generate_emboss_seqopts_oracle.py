#!/usr/bin/env python3
"""A5-7 (B09 F172) oracle: the EMBOSS 6.6.0 sequence-input qualifiers `-sbegin -send -sreverse -slower -supper -scircular -snucleotide -sprotein` (ACD type sequence / seqall / seqset) on the 16 ported programs.

Every program of the repository's EMBOSS ports reads its input sequences through the ACD sequence types, which accept the qualifiers (ajacd.c acdSetSeq / acdSetSeqall / acdSetSeqset, ajseq.c ajSeqSetRange /
ajSeqReverseDo / ajSeqTrim, ajseqread.c seqRead: -slower / -supper are applied at read time, -sbegin / -send are 1-based positions with negative values counted from the end, -sreverse reverse-complements a
NUCLEOTIDE sequence and mirrors the range, -snucleotide / -sprotein force the sequence type). This script runs the installed programs with random qualifier sets (begin / end inside, at and beyond the ends,
negative, zero, reversed ranges; reverse; case options; circular; type overrides, also those the type check refuses) and records what the program REPORTS, so that the port is compared with the binary on the
coordinates of reverse-complemented and ranged input, not with our reading of the C code.

TSV (tests/.../TestData/Alignment/emboss_seqopts_oracle.tsv):  id  program  kind(n|p)  params  seq1  seq2  seq3  quals1  quals2  expected
  quals = begin|end|reverse|lower|upper|circular|nucleotide|protein (0/1 flags); programs with one -sequence use quals1 only.
  needle / water / stretcher / matcher / supermatcher / needleall   expected = the alignment file (-aformat3 pair) after the program header (\\n, \\t, \\\\ escaped), or ERROR
  merger / megamerger          expected = merged sequence (outseq) '@@' report body, or ERROR
  seqmatchall                  expected = sorted 'seq1:seq2:pos1:pos2:len:strand1:strand2' of the -aformat match report or ERROR
  fuzznuc                      expected = sorted 'start:end:strand' (gff) or ERROR;  fuzzpro 'start:end';  fuzztran 'frame:start:end'
  etandem / equicktandem       expected = '|'-joined rows of -origfile (fields ','), '-' = none, or ERROR
  dotmatcher                   expected = 'x1:y1:x2:y2' plot segments, '-' = none, or ERROR
  distmat                      expected = upper triangle (x100, two decimals) joined by ',' or ERROR
Usage: python3 scripts/alignment/generate_emboss_seqopts_oracle.py OUTFILE [per_program=180] [seed]   (needs the EMBOSS 6.6.0 programs on PATH)
"""
import os, random, re, subprocess, sys, tempfile, collections

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from generate_emboss_gap_oracle import Runner, mutate, NUC, AA, esc, clean, gff_hits, tandem_seq, parse_orig  # noqa: E402

HASHES = '#' * 40


def body(text):
    """The alignment file after the program header block (the first two #### lines)."""
    parts = text.split(HASHES + '\n')
    return HASHES.join([''] + parts[2:]) if len(parts) > 2 else text


def rand_pos(rnd, n, allow_zero=True):
    c = rnd.random()
    if c < 0.30 and allow_zero: return 0
    if c < 0.60: return rnd.randint(1, max(1, n))
    if c < 0.78: return -rnd.randint(1, max(1, n))
    if c < 0.88: return rnd.randint(n + 1, n + 20)
    if c < 0.95: return -rnd.randint(n + 1, n + 20)
    return rnd.choice([1, -1, n, -n])


def rand_q(rnd, n, kind, plain=False):
    q = dict(begin=0, end=0, rev=False, lower=False, upper=False, circ=False, nuc=False, prot=False)
    if plain: return q
    if rnd.random() < 0.75: q['begin'] = rand_pos(rnd, n)
    if rnd.random() < 0.75: q['end'] = rand_pos(rnd, n)
    q['rev'] = rnd.random() < 0.40
    q['lower'] = rnd.random() < 0.12
    q['upper'] = rnd.random() < 0.12
    q['circ'] = rnd.random() < 0.06
    if rnd.random() < 0.10:
        q['nuc' if rnd.random() < 0.5 else 'prot'] = True
    return q


def qstr(q):
    return '|'.join(str(int(v)) for v in (q['begin'], q['end'], q['rev'], q['lower'], q['upper'], q['circ'], q['nuc'], q['prot']))


def qargs(q, sfx):
    a = []
    if q['begin']: a += ['-sbegin' + sfx, str(q['begin'])]
    if q['end']: a += ['-send' + sfx, str(q['end'])]
    for key, opt in (('rev', 'sreverse'), ('lower', 'slower'), ('upper', 'supper'), ('circ', 'scircular'), ('nuc', 'snucleotide'), ('prot', 'sprotein')):
        if q[key]: a += ['-' + opt + sfx]
    return a


def pair_inputs(rnd, kind):
    alpha = NUC if kind == 'n' else AA
    n = rnd.choice([12, 20, 30, 45, 60])
    a = ''.join(rnd.choice(alpha) for _ in range(n))
    b = mutate(rnd, a, alpha, rnd.choice([0.0, 0.05, 0.2])) if rnd.random() < 0.8 else ''.join(rnd.choice(alpha) for _ in range(rnd.choice([8, 20, 40])))
    if rnd.random() < 0.3: a = a[rnd.randint(0, 5):]
    if rnd.random() < 0.3: b = b[:max(6, len(b) - rnd.randint(0, 6))]
    return a, b


def finish(p, ok_expected):
    return 'ERROR' if p.returncode else ok_expected


def two(prog, extra_fn=None, multi=False):
    """needle / water / stretcher / matcher / supermatcher / needleall: the pair-format file."""
    def f(rn, rnd, kind):
        alpha = NUC if kind == 'n' else AA
        a, b = pair_inputs(rnd, kind)
        bs = [b]
        if multi:
            bs = [mutate(rnd, a, alpha, rnd.choice([0.02, 0.1])) for _ in range(rnd.choice([1, 2]))]
        qa, qb = rand_q(rnd, len(a), kind), rand_q(rnd, len(b), kind)
        rn.w('a.fa', '>a\n' + a + '\n')
        rn.w('b.fa', ''.join('>b%d\n%s\n' % (i + 1, t) if multi else '>b\n%s\n' % t for i, t in enumerate(bs)))
        extra, params = (extra_fn(rnd) if extra_fn else ([], ''))
        # the matrix: half the rows let the program choose it (default by the type of the first sequence), half give the matching one
        base = [prog, '-asequence', 'a.fa', '-bsequence', 'b.fa', '-auto', '-stdout', '-aformat3', 'pair'] + (['-errfile', 'err.txt'] if multi or prog == 'needleall' else []) + extra
        p = rn.run(base + qargs(qa, '1') + qargs(qb, '2'))
        exp = finish(p, esc(body(p.stdout)))
        return params, [a] + bs + [''] * (2 - len(bs)) if multi else [a, b, ''], qa, qb, exp
    return f


def gaps(rnd):
    g = rnd.choice([(10, 0.5), (10, 0.5), (12, 1), (8, 2), (15, 3)])
    return ['-gapopen', str(g[0]), '-gapextend', str(g[1])], '%s,%s' % g


def igaps(rnd):
    g = rnd.choice([(16, 4), (12, 2), (10, 1), (20, 5)])
    return ['-gapopen', str(g[0]), '-gapextend', str(g[1])], '%s,%s' % g


def matcher_extra(rnd):
    k = rnd.choice([1, 2, 3])
    return ['-alternatives', str(k)], str(k)


def super_extra(rnd):
    w = rnd.choice(['3', '4', '5'])
    return ['-wordlen', w], w


def merger(rn, rnd, kind):
    base = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(50, 100)))
    cut1 = rnd.randint(len(base) // 2, 3 * len(base) // 4)
    cut2 = rnd.randint(len(base) // 4, len(base) // 3)
    a, b = base[:cut1 + 10], base[cut2:]
    qa, qb = rand_q(rnd, len(a), 'n'), rand_q(rnd, len(b), 'n')
    rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n'); rn.rm('o.txt', 'o.fa')
    p = rn.run(['merger', '-asequence', 'a.fa', '-bsequence', 'b.fa', '-outfile', 'o.txt', '-outseq', 'o.fa', '-osformat', 'fasta', '-auto'] + qargs(qa, '1') + qargs(qb, '2'))
    if p.returncode: return '', [a, b, ''], qa, qb, 'ERROR'
    merged = ''.join(l for l in rn.r('o.fa').split('\n') if l and not l.startswith('>'))
    return '', [a, b, ''], qa, qb, esc(merged + '@@' + body(rn.r('o.txt')))


def megamerger(rn, rnd, kind):
    core = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(60, 140)))
    s1 = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(0, 40))) + core
    s2 = core + ''.join(rnd.choice(NUC) for _ in range(rnd.randint(0, 40)))
    a, b = mutate(rnd, s1, NUC, 0.01), mutate(rnd, s2, NUC, 0.01)
    w = rnd.choice([6, 8, 10, 12])
    qa, qb = rand_q(rnd, len(a), 'n'), rand_q(rnd, len(b), 'n')
    rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n'); rn.rm('o.txt', 'o.fa')
    p = rn.run(['megamerger', 'a.fa', 'b.fa', '-wordsize', str(w), '-outseq', 'o.fa', '-osformat', 'fasta', '-outfile', 'o.txt', '-auto'] + qargs(qa, '1') + qargs(qb, '2'))
    if p.returncode: return str(w), [a, b, ''], qa, qb, 'ERROR'
    merged = ''.join(l for l in rn.r('o.fa').split('\n') if l and not l.startswith('>'))
    return str(w), [a, b, ''], qa, qb, esc(merged + '@@' + body(rn.r('o.txt')))


def seqmatchall(rn, rnd, kind):
    alpha = NUC if kind == 'n' else AA
    base = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(20, 50)))
    seqs = [mutate(rnd, base, alpha, 0.1) for _ in range(rnd.choice([2, 3]))]
    k = rnd.choice([3, 4, 5]) if kind == 'n' else rnd.choice([2, 3])
    q = rand_q(rnd, len(base), kind)
    rn.w('s.fa', ''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs))); rn.rm('s.out')
    p = rn.run(['seqmatchall', '-sequence', 's.fa', '-wordsize', str(k), '-outfile', 's.out', '-aformat', 'match', '-auto'] + qargs(q, ''))
    if p.returncode or not os.path.exists(os.path.join(rn.d, 's.out')): return str(k), seqs + [''] * (3 - len(seqs)), q, None, 'ERROR'
    matches = []; cur = None; a = b = None
    for line in rn.r('s.out').split('\n'):
        if line.startswith('# 1: s'): a = int(line.split('s')[-1])
        elif line.startswith('# 2: s'): b = int(line.split('s')[-1]); cur = (a, b)
        m = re.match(r'\s*(\d+) s(\d+)\s+([+-])\s+(\d+)\.\.(\d+)\s+s(\d+)\s+([+-])\s+(\d+)\.\.(\d+)', line)
        if m and cur: matches.append((cur[0], cur[1], int(m.group(4)), int(m.group(5)), int(m.group(8)), int(m.group(9)), int(m.group(1)), m.group(3), m.group(7)))
    matches.sort()
    return str(k), seqs + [''] * (3 - len(seqs)), q, None, ' '.join('%d:%d:%d..%d:%d..%d:%d:%s%s' % t for t in matches) or '-'


def fuzznuc(rn, rnd, kind):
    s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(25, 80)))
    st = rnd.randint(0, len(s) - 6)
    pat = rnd.choice([s[st:st + rnd.randint(3, 6)], s[st:st + 4].replace(s[st], 'N', 1), s[st:st + 3] + '[AC]' + s[st + 4:st + 6], 'AC' + 'N' * rnd.randint(1, 2) + 'GT', 'RYR', s[st:st + 3]])
    q = rand_q(rnd, len(s), 'n')
    mm = rnd.choice([0, 0, 1]); comp = rnd.choice([0, 1])
    rn.w('a.fa', '>a\n' + s + '\n')
    p = rn.run(['fuzznuc', '-sequence', 'a.fa', '-pattern', pat, '-pmismatch', str(mm), '-complement', 'Y' if comp else 'N', '-rformat', 'gff', '-stdout', '-auto'] + qargs(q, ''))
    params = '%s|%d|%d' % (pat, mm, comp)
    if p.returncode: return params, [s, '', ''], q, None, 'ERROR'
    hits = sorted(gff_hits(p.stdout, lambda t: '%s:%s:%s' % (t[3], t[4], t[6])), key=lambda x: (int(x.split(':')[0]), int(x.split(':')[1]), x.split(':')[2]))
    return params, [s, '', ''], q, None, ' '.join(hits) or '-'


def fuzzpro(rn, rnd, kind):
    s = ''.join(rnd.choice(AA) for _ in range(rnd.randint(25, 80)))
    st = rnd.randint(0, len(s) - 5)
    pat = rnd.choice(['-'.join(s[st:st + 2]), '-'.join(s[st:st + 3]), s[st] + '-x-' + s[st + 2], '[' + s[st] + s[st + 1] + ']-' + s[st + 2], 'x(2)-' + s[st + 3]])
    q = rand_q(rnd, len(s), 'p')
    rn.w('a.fa', '>a\n' + s + '\n')
    p = rn.run(['fuzzpro', '-sequence', 'a.fa', '-pattern', pat, '-rformat', 'gff', '-stdout', '-auto'] + qargs(q, ''))
    if p.returncode: return pat, [s, '', ''], q, None, 'ERROR'
    hits = sorted(gff_hits(p.stdout, lambda t: '%s:%s' % (t[3], t[4])), key=lambda x: tuple(int(v) for v in x.split(':')))
    return pat, [s, '', ''], q, None, ' '.join(hits) or '-'


def fuzztran(rn, rnd, kind):
    s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(30, 90)))
    pat = rnd.choice(['L-x-S', 'S-S', 'G-R', 'x(2)-K', '[LV]-T', 'M', 'A-A', 'T-x(1,2)-G'])
    q = rand_q(rnd, len(s), 'n')
    fr = rnd.choice(['F', '6', '1', '-2'])
    rn.w('a.fa', '>a\n' + s + '\n')
    p = rn.run(['fuzztran', '-sequence', 'a.fa', '-pattern', pat, '-frame', fr, '-table', '0', '-rformat', 'gff', '-stdout', '-auto'] + qargs(q, ''))
    params = pat + '|' + fr
    if p.returncode: return params, [s, '', ''], q, None, 'ERROR'
    hits = []
    for line in p.stdout.split('\n'):
        if not line or line.startswith('#'): continue
        t = line.split('\t')
        if len(t) >= 9:
            f = re.search(r'\*frame (-?\d+)', t[8]).group(1)
            hits.append('%s:%s:%s' % (f, t[3], t[4]))
    hits.sort(key=lambda x: tuple(int(v) for v in x.split(':')))
    return params, [s, '', ''], q, None, ' '.join(hits) or '-'


def tandem(prog):
    def f(rn, rnd, kind):
        d = tandem_seq(rnd)
        q = rand_q(rnd, len(d), 'n')
        rn.w('s.fa', '>s\n' + d + '\n'); rn.rm('o.rep', 'o.orig')
        if prog == 'etandem':
            mn = rnd.choice([2, 3]); mx = mn + rnd.randint(1, 5); th = rnd.choice([8, 12, 15, 20]); mm = rnd.choice([0, 1]); un = rnd.choice([0, 1])
            args = ['-minrepeat', str(mn), '-maxrepeat', str(mx), '-threshold', str(th)] + (['-mismatch'] if mm else []) + (['-uniform'] if un else [])
            params, ncol = '%d,%d,%d,%d,%d' % (mn, mx, th, mm, un), 7
        else:
            mx = rnd.choice([4, 6, 8]); th = rnd.choice([5, 8, 12, 20])
            args = ['-maxrepeat', str(mx), '-threshold', str(th)]
            params, ncol = '%d,%d' % (mx, th), 5
        p = rn.run([prog, '-sequence', 's.fa', '-outfile', 'o.rep', '-origfile', 'o.orig', '-auto'] + args + qargs(q, ''))
        if p.returncode: return params, [d, '', ''], q, None, 'ERROR'
        return params, [d, '', ''], q, None, parse_orig(rn.r('o.orig'), ncol) if os.path.exists(os.path.join(rn.d, 'o.orig')) else '-'
    return f


def dotmatcher(rn, rnd, kind):
    a, b = pair_inputs(rnd, kind)
    la, lb = len(a), len(b)
    w = rnd.randint(3, min(8, la - 1, lb - 1))
    mat = rnd.choice(['EBLOSUM62', 'EPAM250']) if kind == 'p' else rnd.choice(['EDNAFULL', 'EDNASIMPLE'])
    t = rnd.randint(max(1, w * 5 // 3), w * 5) if kind == 'n' else rnd.randint(w * 2, w * 5)
    qa, qb = rand_q(rnd, la, kind), rand_q(rnd, lb, kind)
    rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n')
    for f in os.listdir(rn.d):
        if f.startswith('dm'): os.remove(os.path.join(rn.d, f))
    p = rn.run(['dotmatcher', 'a.fa', 'b.fa', '-windowsize', str(w), '-threshold', str(t), '-matrixfile', mat, '-graph', 'data', '-goutfile', 'dm', '-auto'] + qargs(qa, '1') + qargs(qb, '2'))
    files = [f for f in os.listdir(rn.d) if f.startswith('dm') and f.endswith('.dat')]
    params = '%s,%d,%d' % (mat, w, t)
    if p.returncode or not files: return params, [a, b, ''], qa, qb, 'ERROR'
    segs = []
    for line in rn.r(files[0]).split('\n'):
        if line.startswith('Rectangle'): break
        m = re.match(r'Line x1 (\S+) y1 (\S+) x2 (\S+) y2 (\S+) colour', line)
        if m: segs.append(':'.join(str(int(float(g))) for g in m.groups()))
    return params, [a, b, ''], qa, qb, ' '.join(segs) or '-'


def distmat(rn, rnd, kind):
    alpha = NUC if kind == 'n' else AA
    n = rnd.randint(10, 30)
    base = ''.join(rnd.choice(alpha) for _ in range(n))
    rows = []
    for _ in range(rnd.choice([2, 3])):
        r = list(mutate(rnd, base, alpha, 0.1))[:n]
        r += [rnd.choice(alpha)] * (n - len(r))
        for _ in range(rnd.randint(0, 3)): r[rnd.randrange(n)] = rnd.choice('-')
        rows.append(''.join(r))
    q = rand_q(rnd, n, kind)
    meth = rnd.choice([0, 1, 2]) if kind == 'p' else rnd.choice([0, 1, 2, 3, 4, 5])
    rn.w('a.fa', ''.join('>s%d\n%s\n' % (i + 1, r) for i, r in enumerate(rows))); rn.rm('o.txt')
    p = rn.run(['distmat', '-sequence', 'a.fa', '-outfile', 'o.txt', '-auto', '-nucmethod', str(meth), '-protmethod', str(meth if meth <= 2 else 0)] + qargs(q, ''))
    if p.returncode or not os.path.exists(os.path.join(rn.d, 'o.txt')): return str(meth), rows + [''] * (3 - len(rows)), q, None, 'ERROR'
    lines = rn.r('o.txt').split('\n')
    starts = [i for i, l in enumerate(lines) if re.match(r'^\s+1\s+2', l)]
    if not starts: return str(meth), rows + [''] * (3 - len(rows)), q, None, 'ERROR'
    nn = len(rows); vals = []
    for j in range(nn - 1):
        f = lines[starts[0] + 1 + j].split('\t')
        cols = f[1:nn + 1]
        vals += [cols[i].strip() for i in range(j + 1, nn)]
    return str(meth), rows + [''] * (3 - len(rows)), q, None, ','.join(vals)


PROGRAMS = [
    ('needle', two('needle', gaps), 'np'), ('water', two('water', gaps), 'np'), ('stretcher', two('stretcher', igaps), 'np'), ('matcher', two('matcher', matcher_extra), 'np'),
    ('supermatcher', two('supermatcher', super_extra, multi=True), 'np'), ('needleall', two('needleall', lambda r: ([], ''), multi=True), 'np'),
    ('merger', merger, 'n'), ('megamerger', megamerger, 'n'), ('seqmatchall', seqmatchall, 'np'),
    ('fuzznuc', fuzznuc, 'n'), ('fuzzpro', fuzzpro, 'p'), ('fuzztran', fuzztran, 'n'), ('etandem', tandem('etandem'), 'n'), ('equicktandem', tandem('equicktandem'), 'n'),
    ('dotmatcher', dotmatcher, 'np'), ('distmat', distmat, 'np'),
]


def main(outfile, per, seed, only=None):
    rnd = random.Random(seed)
    rows = []
    counts = collections.Counter(); errors = collections.Counter()
    with tempfile.TemporaryDirectory() as d:
        rn = Runner(d)
        for name, fn, kinds in PROGRAMS:
            if only and name not in only: continue
            made = attempts = 0
            while made < per and attempts < per * 6:
                attempts += 1
                kind = kinds[made % len(kinds)] if len(kinds) > 1 else kinds
                res = fn(rn, rnd, kind)
                if res is None: continue
                params, seqs, qa, qb, expected = res
                if any('\t' in s or '\n' in s for s in seqs): continue
                if expected == 'ERROR' and errors[name] > per // 6: continue      # keep the refused combinations a minority
                rows.append('\t'.join([str(len(rows)), name, kind, params] + seqs + [qstr(qa), qstr(qb) if qb else qstr(dict(begin=0, end=0, rev=False, lower=False, upper=False, circ=False, nuc=False, prot=False)), expected]))
                made += 1
                if expected == 'ERROR': errors[name] += 1
            counts[name] = made
    with open(outfile, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('id\tprogram\tkind\tparams\tseq1\tseq2\tseq3\tquals1\tquals2\texpected\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows', dict(counts), 'errors', dict(errors))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 180, int(sys.argv[3]) if len(sys.argv) > 3 else 20261127, sys.argv[4].split(',') if len(sys.argv) > 4 else None)
