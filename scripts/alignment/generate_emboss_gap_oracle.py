#!/usr/bin/env python3
"""A4-3 (B09 F153) oracle: what the installed EMBOSS 6.6.0 programs do with a PLAIN input sequence that contains gap characters and other non-letters.

Measured (scripts/alignment/measure_emboss_reader_chars.py): EMBOSS drops every character that is not a letter or `*` from a plain sequence before any length check, coordinate or score (the file reader
keeps letters and `* . ~ ? # + -`, then ajSeqTypeCheckIn of every ungapped ACD type calls ajStrRemoveGap, which keeps `A-Za-z*`; ajax/core/ajseqread.c seqAppend, ajax/core/ajseqtype.c
ajSeqTypeCheckIn, ajax/core/ajstr.c ajStrRemoveGap). The alignment readers of consambig / showalign (ACD type gapany) keep `. ~ -` as gaps and read `? # +` as an unknown residue.
This file records, for EVERY ported program that reads plain sequences, the program's output on sequences that DO contain `- . ~ ? # + _`, digits, blanks, quotes, brackets, `*` (protein) and
some non-ASCII letters, so that the tests lock the shared reader (Seqeron.Genomics.Infrastructure.EmbossSequenceReader) against the binary and not against our own reading of its source.

TSV (tests/.../TestData/Alignment/emboss_gap_oracle.tsv):  id  program  kind(n|p)  params  seq1  seq2  seq3  expected
  needle / water / stretcher   expected = score|row1|row2   (-aformat3 pair score, -aformat3 fasta rows; stretcher: integer-valued score text as printed)
  matcher                      params K; expected = ';'-joined score|start1|start2|end1|end2|row1|row2 (0-based starts and inclusive ends) or ABORT
  merger / megamerger          params W (megamerger word size, 0 for merger); expected = merged sequence (outseq) or '-' (megamerger: the program wrote no report)
  seqmatchall                  params W; seq1..seq3 (seq3 may be empty); expected = sorted 'i:j:pos1:pos2:len' (0-based) or '-'
  supermatcher / needleall     seq1 = query, seq2 / seq3 = targets (seq3 may be empty); expected = ';'-joined row1|row2 of the alignments in the order the program writes them, or '-'
  fuzznuc                      params pattern|mismatch|complement(0/1); expected = L=<cleaned length>;sorted 'start:end:strand' (1-based, as printed)
  fuzzpro                      params pattern; expected = L=<cleaned length>;sorted 'start:end'
  fuzztran                     params pattern|frame(F); expected = L=<cleaned length>;sorted 'frame:start:end'
  etandem / equicktandem       params as generate_emboss_tandem_oracle.py; expected = '|'-joined rows (fields ','), '-' = none
  dotmatcher                   params matrix,window,threshold; expected = 'x1:y1:x2:y2' segments of the plot, ' '-joined, '-' = none
  consambig                    seq1..seq3 = aligned rows (seq3 may be empty); expected = the consensus (outseq), '-' when the program stops
  showalign                    seq1..seq3 rows; expected = the text of -outfile with escapes (\\n, \\t, \\\\)
  distmat                      params = nucleotide / protein method (-nucmethod / -protmethod); seq1..seq3 rows; expected = the upper triangle of the printed matrix (x100, two decimals) joined by ','
Never aborts: a failing program call is reported on stderr and the row skipped (the programs' own errors are not the subject here).
Usage: python3 scripts/alignment/generate_emboss_gap_oracle.py OUTFILE [per_program=240] [seed]      (needs the EMBOSS 6.6.0 programs on PATH)
       python3 scripts/alignment/measure_emboss_reader_chars.py                  the per-character measurement behind the rule (which of the 252 characters every program drops)
"""
import os, random, re, subprocess, sys, tempfile, collections

NUC = 'ACGT'
AA = 'ACDEFGHIKLMNPQRSTVWY'
DIRTY = list("-.~?_# +0123456789;:'\"()[]{}|!@$%^&=<\\/`,") + ['\x0b', '\x0c', '\x7f', 'é', 'Ω', 'ß']
DIRTY_NOTAB = [c for c in DIRTY if c not in '\t\n\r>']


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def clean(s):
    return ''.join(c for c in s if (c.isascii() and c.isalpha()) or c == '*')


def mutate(rnd, s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate: out.append(rnd.choice(alpha))
        elif r < rate * 1.4: continue
        elif r < rate * 1.8: out.append(c); out.extend(rnd.choice(alpha) for _ in range(rnd.randint(1, 3)))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def dirt(rnd, s, star=False, lower=True):
    """Inserts 1..8 non-sequence characters (also at both ends) and sometimes lower-cases a stretch."""
    chars = list(s)
    pool = DIRTY_NOTAB
    n = rnd.choice([1, 1, 2, 3, 4, 6, 8, 12])
    for _ in range(n):
        pos = rnd.choice([0, len(chars), rnd.randint(0, len(chars))])
        chars.insert(pos, rnd.choice(pool))
    if star and rnd.random() < 0.5:
        chars.insert(rnd.randint(0, len(chars)), '*')
    out = ''.join(chars)
    if lower and rnd.random() < 0.25:
        out = out.lower() if rnd.random() < 0.5 else ''.join(c.lower() if rnd.random() < 0.5 else c for c in out)
    return out


class Runner:
    def __init__(self, d):
        self.d = d

    def w(self, name, text):
        with open(os.path.join(self.d, name), 'w', encoding='utf-8', newline='\n') as f:
            f.write(text)

    def r(self, name):
        p = os.path.join(self.d, name)
        return open(p, encoding='utf-8', errors='replace').read() if os.path.exists(p) else ''

    def rm(self, *names):
        for n in names:
            p = os.path.join(self.d, n)
            if os.path.exists(p):
                os.remove(p)

    def run(self, args):
        return subprocess.run(args, cwd=self.d, capture_output=True, text=True, errors='replace')


# ---------- the programs: each function returns (params, seqs, expected) or None -----------------------------------------------------------------------

def pair_inputs(rnd, kind, star):
    alpha = NUC if kind == 'n' else AA
    n = rnd.choice([8, 12, 20, 30, 45, 60])
    a = ''.join(rnd.choice(alpha) for _ in range(n))
    b = mutate(rnd, a, alpha, rnd.choice([0.0, 0.05, 0.2])) if rnd.random() < 0.8 else ''.join(rnd.choice(alpha) for _ in range(rnd.choice([5, 20, 40])))
    return dirt(rnd, a, star and kind == 'p'), dirt(rnd, b, star and kind == 'p')


def aligner(prog, extra=()):
    def f(rn, rnd, kind):
        a, b = pair_inputs(rnd, kind, star=True)
        if len(clean(a)) < 2 or len(clean(b)) < 2: return None
        rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n')
        base = [prog, '-asequence', 'a.fa', '-bsequence', 'b.fa', '-auto', '-stdout'] + list(extra)
        p = rn.run(base + ['-aformat3', 'pair']); q = rn.run(base + ['-aformat3', 'fasta'])
        m = re.search(r'^# Score: (\S+)$', p.stdout, re.M)
        recs = re.findall(r'>\S+[^\n]*\n([^>]*)', q.stdout)
        if p.returncode or q.returncode or not m or len(recs) != 2: return None
        return '', [a, b, ''], '%s|%s|%s' % (m.group(1), recs[0].replace('\n', ''), recs[1].replace('\n', ''))
    return f


def matcher(rn, rnd, kind):
    a, b = pair_inputs(rnd, kind, star=True)
    if len(clean(a)) < 6 or len(clean(b)) < 6: return None
    k = rnd.choice([1, 2, 3])
    rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n')
    mat, typ = (('EDNAFULL', ['-snucleotide1', '-snucleotide2']) if kind == 'n' else ('EBLOSUM62', ['-sprotein1', '-sprotein2']))
    rn.rm('o.txt')
    p = rn.run(['matcher', 'a.fa', 'b.fa', '-gapopen', '14', '-gapextend', '4', '-datafile', mat, '-alternatives', str(k), '-aformat3', 'srspair', '-outfile', 'o.txt', '-auto'] + typ)
    if p.returncode != 0:
        return str(k), [a, b, ''], 'ABORT'
    out = []
    for part in re.split(r'#=======================================\n#\n# Aligned_sequences', rn.r('o.txt'))[1:]:
        sc = int(re.search(r'# Score: (\d+)', part).group(1))
        ra = rb = ''; sa = sb = None
        for line in part.split('\n'):
            m_ = re.match(r'^a\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sa = sa or int(m_.group(1)); ra += m_.group(2); continue
            m_ = re.match(r'^b\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sb = sb or int(m_.group(1)); rb += m_.group(2)
        ea = sa - 1 + sum(1 for x in ra if x != '-'); eb = sb - 1 + sum(1 for x in rb if x != '-')
        out.append(f'{sc}|{sa - 1}|{sb - 1}|{ea - 1}|{eb - 1}|{ra}|{rb}')
    return str(k), [a, b, ''], ';'.join(out)


def overlapping(rnd, kind='n'):
    base = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(50, 110)))
    cut1 = rnd.randint(len(base) // 2, 3 * len(base) // 4)
    cut2 = rnd.randint(len(base) // 4, len(base) // 3)
    return (mutate(rnd, base[:cut1 + 10], NUC, 0.01) if rnd.random() < 0.3 else base[:cut1 + 10]), base[cut2:]


def merger(rn, rnd, kind):
    a, b = overlapping(rnd)
    a, b = dirt(rnd, a), dirt(rnd, b)
    rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n'); rn.rm('o.txt', 'o.fa')
    p = rn.run(['merger', '-asequence', 'a.fa', '-bsequence', 'b.fa', '-outfile', 'o.txt', '-outseq', 'o.fa', '-auto'])
    if p.returncode or not rn.r('o.fa'): return None
    merged = ''.join(l for l in rn.r('o.fa').split('\n') if l and not l.startswith('>'))
    return '0', [a, b, ''], merged


def megamerger(rn, rnd, kind):
    core = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(60, 160)))
    s1 = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(0, 40))) + core
    s2 = core + ''.join(rnd.choice(NUC) for _ in range(rnd.randint(0, 40)))
    if rnd.random() < 0.2: s2 = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(30, 90)))
    a, b = dirt(rnd, mutate(rnd, s1, NUC, 0.01)), dirt(rnd, mutate(rnd, s2, NUC, 0.01))
    w = rnd.choice([6, 8, 10, 12])
    rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n'); rn.rm('o.txt', 'o.fa')
    p = rn.run(['megamerger', 'a.fa', 'b.fa', '-wordsize', str(w), '-outseq', 'o.fa', '-osformat', 'fasta', '-outfile', 'o.txt', '-auto'])
    if p.returncode: return None
    if not rn.r('o.txt').strip(): return str(w), [a, b, ''], '-'
    return str(w), [a, b, ''], ''.join(l for l in rn.r('o.fa').split('\n') if l and not l.startswith('>'))


def seqmatchall(rn, rnd, kind):
    alpha = NUC if kind == 'n' else AA
    base = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(20, 50)))
    seqs = [dirt(rnd, mutate(rnd, base, alpha, 0.1), star=kind == 'p') for _ in range(rnd.choice([2, 3]))]
    k = rnd.choice([3, 4, 5, 6]) if kind == 'n' else rnd.choice([2, 3, 4])
    if any(len(clean(s)) <= k for s in seqs): return None
    rn.w('s.fa', ''.join('>s%d\n%s\n' % (i, s) for i, s in enumerate(seqs))); rn.rm('s.out')
    p = rn.run(['seqmatchall', '-sequence', 's.fa', '-wordsize', str(k), '-outfile', 's.out', '-aformat', 'match', '-auto'])
    if p.returncode or not os.path.exists(os.path.join(rn.d, 's.out')): return None
    matches = []; cur = None; a = b = None
    for line in rn.r('s.out').split('\n'):
        if line.startswith('# 1: s'): a = int(line.split('s')[-1])
        elif line.startswith('# 2: s'): b = int(line.split('s')[-1]); cur = (a, b)
        m = re.match(r'\s*(\d+) s(\d+)\s+\+\s+(\d+)\.\.(\d+)\s+s(\d+)\s+\+\s+(\d+)\.\.(\d+)', line)
        if m and cur: matches.append((cur[0], cur[1], int(m.group(3)) - 1, int(m.group(6)) - 1, int(m.group(1))))
    matches.sort()
    return str(k), seqs + [''] * (3 - len(seqs)), ' '.join('%d:%d:%d:%d:%d' % t for t in matches) or '-'


def pairs_fasta(text):
    recs = re.findall(r'>(\S+)[^\n]*\n([^>]*)', text)
    rows = [r[1].replace('\n', '') for r in recs]
    return ';'.join('%s|%s' % (rows[i], rows[i + 1]) for i in range(0, len(rows) - 1, 2)) or '-'


def multi(prog, extra):
    def f(rn, rnd, kind):
        alpha = NUC if kind == 'n' else AA
        core = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(30, 70)))
        q = dirt(rnd, mutate(rnd, core, alpha, 0.03), star=kind == 'p')
        ts = [dirt(rnd, mutate(rnd, core, alpha, rnd.choice([0.02, 0.1])), star=kind == 'p') for _ in range(rnd.choice([1, 2]))]
        if len(clean(q)) < 8 or any(len(clean(t)) < 8 for t in ts): return None
        rn.w('a.fa', '>q1\n' + q + '\n'); rn.w('b.fa', ''.join('>t%d\n%s\n' % (i + 1, t) for i, t in enumerate(ts))); rn.rm('err.txt')
        p = rn.run([prog, '-asequence', 'a.fa', '-bsequence', 'b.fa', '-auto', '-stdout', '-errfile', 'err.txt', '-aformat3', 'fasta'] + extra)
        if p.returncode: return None
        return '', [q] + ts + [''] * (2 - len(ts)), pairs_fasta(p.stdout)
    return f


def gff_hits(text, cols):
    hits = []
    for line in text.split('\n'):
        if not line or line.startswith('#'): continue
        t = line.split('\t')
        if len(t) >= 8: hits.append(cols(t))
    return hits


def fuzznuc(rn, rnd, kind):
    s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(25, 80)))
    st = rnd.randint(0, len(s) - 6)
    pat = rnd.choice([s[st:st + rnd.randint(3, 6)], s[st:st + 4].replace(s[st], 'N', 1), s[st:st + 3] + '[AC]' + s[st + 4:st + 6], 'AC' + 'N' * rnd.randint(1, 2) + 'GT', 'RYR', s[st:st + 3]])
    d = dirt(rnd, s)
    mm = rnd.choice([0, 0, 1]); comp = rnd.choice([0, 1])
    rn.w('a.fa', '>a\n' + d + '\n')
    p = rn.run(['fuzznuc', '-sequence', 'a.fa', '-pattern', pat, '-pmismatch', str(mm), '-complement', 'Y' if comp else 'N', '-rformat', 'gff', '-stdout', '-auto'])
    if p.returncode: return None
    hits = sorted(gff_hits(p.stdout, lambda t: '%s:%s:%s' % (t[3], t[4], t[6])), key=lambda x: (int(x.split(':')[0]), int(x.split(':')[1]), x.split(':')[2]))
    return '%s|%d|%d' % (pat, mm, comp), [d, '', ''], 'L=%d;%s' % (len(clean(d)), ' '.join(hits) or '-')


def fuzzpro(rn, rnd, kind):
    s = ''.join(rnd.choice(AA) for _ in range(rnd.randint(25, 80)))
    st = rnd.randint(0, len(s) - 5)
    pat = rnd.choice(['-'.join(s[st:st + 2]), '-'.join(s[st:st + 3]), s[st] + '-x-' + s[st + 2], '[' + s[st] + s[st + 1] + ']-' + s[st + 2], 'x(2)-' + s[st + 3]])
    d = dirt(rnd, s)
    rn.w('a.fa', '>a\n' + d + '\n')
    p = rn.run(['fuzzpro', '-sequence', 'a.fa', '-pattern', pat, '-rformat', 'gff', '-stdout', '-auto'])
    if p.returncode: return None
    hits = sorted(gff_hits(p.stdout, lambda t: '%s:%s' % (t[3], t[4])), key=lambda x: tuple(int(v) for v in x.split(':')))
    return pat, [d, '', ''], 'L=%d;%s' % (len(clean(d)), ' '.join(hits) or '-')


def fuzztran(rn, rnd, kind):
    s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(30, 90)))
    pat = rnd.choice(['L-x-S', 'S-S', 'G-R', 'x(2)-K', '[LV]-T', 'M', 'A-A', 'T-x(1,2)-G'])
    d = dirt(rnd, s)
    rn.w('a.fa', '>a\n' + d + '\n')
    p = rn.run(['fuzztran', '-sequence', 'a.fa', '-pattern', pat, '-frame', 'F', '-table', '0', '-rformat', 'gff', '-stdout', '-auto'])
    if p.returncode: return None
    hits = []
    for line in p.stdout.split('\n'):
        if not line or line.startswith('#'): continue
        t = line.split('\t')
        if len(t) >= 9:
            fr = re.search(r'\*frame (-?\d+)', t[8]).group(1)
            hits.append('%s:%s:%s' % (fr, t[3], t[4]))
    hits.sort(key=lambda x: tuple(int(v) for v in x.split(':')))
    return pat + '|F', [d, '', ''], 'L=%d;%s' % (len(clean(d)), ' '.join(hits) or '-')


def tandem_seq(rnd):
    unit = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(2, 6)))
    s = ''.join(rnd.choice(NUC) for _ in range(rnd.randint(3, 25))) + mutate(rnd, unit * rnd.randint(3, 9), NUC, 0.05) + ''.join(rnd.choice(NUC) for _ in range(rnd.randint(3, 25)))
    return s


def parse_orig(text, ncol):
    rows = []
    for line in text.split('\n'):
        t = line.split()
        if len(t) >= ncol: rows.append(','.join(t[:ncol]))
    return '|'.join(rows) if rows else '-'


def tandem(prog):
    def f(rn, rnd, kind):
        d = dirt(rnd, tandem_seq(rnd))
        rn.w('s.fa', '>s\n' + d + '\n'); rn.rm('o.rep', 'o.orig')
        if prog == 'etandem':
            mn = rnd.choice([2, 3]); mx = mn + rnd.randint(1, 5); th = rnd.choice([8, 12, 15, 20]); mm = rnd.choice([0, 1]); un = rnd.choice([0, 1])
            args = ['-minrepeat', str(mn), '-maxrepeat', str(mx), '-threshold', str(th)] + (['-mismatch'] if mm else []) + (['-uniform'] if un else [])
            params, ncol = '%d,%d,%d,%d,%d' % (mn, mx, th, mm, un), 7
        else:
            mx = rnd.choice([4, 6, 8]); th = rnd.choice([5, 8, 12, 20])
            args = ['-maxrepeat', str(mx), '-threshold', str(th)]
            params, ncol = '%d,%d' % (mx, th), 5
        p = rn.run([prog, '-sequence', 's.fa', '-outfile', 'o.rep', '-origfile', 'o.orig', '-auto'] + args)
        if p.returncode or not os.path.exists(os.path.join(rn.d, 'o.orig')): return None
        return params, [d, '', ''], parse_orig(rn.r('o.orig'), ncol)
    return f


def dotmatcher(rn, rnd, kind):
    a, b = pair_inputs(rnd, kind, star=True)
    la, lb = len(clean(a)), len(clean(b))
    if la < 8 or lb < 8: return None
    w = rnd.randint(3, min(8, la - 1, lb - 1))
    mat = rnd.choice(['EBLOSUM62', 'EPAM250']) if kind == 'p' else rnd.choice(['EDNAFULL', 'EDNASIMPLE'])
    t = rnd.randint(max(1, w * 5 // 3), w * 5) if kind == 'n' else rnd.randint(w * 2, w * 5)
    rn.w('a.fa', '>a\n' + a + '\n'); rn.w('b.fa', '>b\n' + b + '\n')
    for f in os.listdir(rn.d):
        if f.startswith('dm'): os.remove(os.path.join(rn.d, f))
    p = rn.run(['dotmatcher', 'a.fa', 'b.fa', '-windowsize', str(w), '-threshold', str(t), '-matrixfile', mat, '-graph', 'data', '-goutfile', 'dm', '-auto'])
    files = [f for f in os.listdir(rn.d) if f.startswith('dm') and f.endswith('.dat')]
    if p.returncode or not files: return None
    segs = []
    for line in rn.r(files[0]).split('\n'):
        if line.startswith('Rectangle'): break
        m = re.match(r'Line x1 (\S+) y1 (\S+) x2 (\S+) y2 (\S+) colour', line)
        if m: segs.append(':'.join(str(int(float(g))) for g in m.groups()))
    return '%s,%d,%d' % (mat, w, t), [a, b, ''], ' '.join(segs) or '-'


def alignment_rows(rnd, kind, dirty_chars):
    alpha = NUC if kind == 'n' else AA
    n = rnd.randint(8, 30)
    base = ''.join(rnd.choice(alpha) for _ in range(n))
    rows = []
    for _ in range(rnd.choice([2, 3])):
        r = list(mutate(rnd, base, alpha, 0.1))
        r = [c for c in r][:n]
        r += [rnd.choice(alpha)] * (n - len(r))
        for _ in range(rnd.randint(0, 4)):
            r[rnd.randrange(n)] = rnd.choice('-.~')
        rows.append(''.join(r))
    # the same number of non-sequence characters at the same places in every row keeps the cleaned rows aligned; some rows get the unknown-residue characters ? # + in place
    out = [list(r) for r in rows]
    for _ in range(rnd.choice([0, 1, 2, 3])):
        pos = rnd.randint(0, n)
        for r in out: r.insert(pos, rnd.choice(dirty_chars))
    for r in out:
        if rnd.random() < 0.3 and r: r[rnd.randrange(len(r))] = rnd.choice('?#+')
    return [''.join(r) for r in out]


ALIGN_DIRTY = [c for c in DIRTY_NOTAB if c.isascii() and c not in '-.~?#+*']      # non-ASCII bytes: the binary's isalpha() keeps some of them under some locales


def consambig(rn, rnd, kind):
    rows = alignment_rows(rnd, kind, ALIGN_DIRTY)
    rn.w('a.fa', ''.join('>r%d\n%s\n' % (i, r) for i, r in enumerate(rows))); rn.rm('o.fa')
    p = rn.run(['consambig', '-sequence', 'a.fa', '-outseq', 'o.fa', '-osformat', 'fasta', '-auto'])
    if p.returncode or not rn.r('o.fa'): return None
    cons = ''.join(l for l in rn.r('o.fa').split('\n') if l and not l.startswith('>'))
    return '', rows + [''] * (3 - len(rows)), cons or '-'


def showalign(rn, rnd, kind):
    rows = alignment_rows(rnd, kind, ALIGN_DIRTY)
    rn.w('a.fa', ''.join('>r%d\n%s\n' % (i, r) for i, r in enumerate(rows))); rn.rm('o.txt')
    p = rn.run(['showalign', '-sequence', 'a.fa', '-outfile', 'o.txt', '-width', '30', '-auto'])
    if p.returncode or not os.path.exists(os.path.join(rn.d, 'o.txt')): return None
    return '', rows + [''] * (3 - len(rows)), esc(rn.r('o.txt'))


def distmat(rn, rnd, kind):
    rows = alignment_rows(rnd, kind, ALIGN_DIRTY)
    rn.w('a.fa', ''.join('>s%d\n%s\n' % (i + 1, r) for i, r in enumerate(rows))); rn.rm('o.txt')
    meth = rnd.choice([0, 1, 2]) if kind == 'p' else rnd.choice([0, 1, 2, 3, 4, 5])
    p = rn.run(['distmat', '-sequence', 'a.fa', '-outfile', 'o.txt', '-auto', '-nucmethod', str(meth), '-protmethod', str(meth if meth <= 2 else 0)])
    if p.returncode != 0 or not os.path.exists(os.path.join(rn.d, 'o.txt')): return None
    n = len(rows)
    lines = rn.r('o.txt').split('\n')
    starts = [i for i, l in enumerate(lines) if re.match(r'^\s+1\s+2', l)]
    if not starts: return None
    start = starts[0] + 1
    vals = []
    for j in range(n - 1):
        f = lines[start + j].split('\t')
        cols = f[1:n + 1]
        vals += [cols[i].strip() for i in range(j + 1, n)]
    return str(meth), rows + [''] * (3 - len(rows)), ','.join(vals)


PROGRAMS = [
    ('needle', aligner('needle'), 'np'), ('water', aligner('water'), 'np'), ('stretcher', aligner('stretcher'), 'np'), ('matcher', matcher, 'np'),
    ('merger', merger, 'n'), ('megamerger', megamerger, 'n'), ('seqmatchall', seqmatchall, 'np'),
    ('supermatcher', multi('supermatcher', ['-wordlen', '4']), 'np'), ('needleall', multi('needleall', []), 'np'),
    ('fuzznuc', fuzznuc, 'n'), ('fuzzpro', fuzzpro, 'p'), ('fuzztran', fuzztran, 'n'), ('etandem', tandem('etandem'), 'n'), ('equicktandem', tandem('equicktandem'), 'n'),
    ('dotmatcher', dotmatcher, 'np'), ('consambig', consambig, 'np'), ('showalign', showalign, 'np'), ('distmat', distmat, 'np'),
]


def main(outfile, per, seed):
    rnd = random.Random(seed)
    rows = []
    counts = collections.Counter()
    with tempfile.TemporaryDirectory() as d:
        rn = Runner(d)
        for name, fn, kinds in PROGRAMS:
            made = attempts = 0
            while made < per and attempts < per * 8:
                attempts += 1
                kind = kinds[made % len(kinds)] if len(kinds) > 1 else kinds
                res = fn(rn, rnd, kind)
                if res is None:
                    continue
                params, seqs, expected = res
                if any('\t' in s or '\n' in s or '\r' in s for s in seqs): continue
                rows.append('\t'.join([str(len(rows)), name, kind, params] + seqs + [expected]))
                made += 1
            counts[name] = made
    with open(outfile, 'w', encoding='utf-8', newline='\n') as fh:
        fh.write('id\tprogram\tkind\tparams\tseq1\tseq2\tseq3\texpected\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'rows', dict(counts))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 240, int(sys.argv[3]) if len(sys.argv) > 3 else 20261126)
