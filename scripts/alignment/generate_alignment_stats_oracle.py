#!/usr/bin/env python3
"""ALIGN-STATS-001 oracle: EMBOSS needle/water 6.6.0 statistics + markup, cross-checked with Biopython 1.88.

For random DNA / IUPAC-DNA / protein pairs (and the published HBA_HUMAN vs HBB_HUMAN example) this runs
`needle` / `water` (apt package `emboss`), parses the aligned rows, the srspair markup line and the
header (Length / Identity / Similarity / Gaps), and records Biopython `Alignment.counts(matrix)`
(mismatches, open_gaps, identities, positives, gaps) for the same rows. Output: a TSV embedded in
tests/Seqeron/Seqeron.Genomics.Tests/TestData/Alignment/needle_water_stats.tsv

Usage: python3 scripts/alignment/generate_alignment_stats_oracle.py [n_random=300] > needle_water_stats.tsv
"""
import random, subprocess, re, sys, tempfile
import numpy as np
from Bio.Align import Alignment, substitution_matrices as sm

rnd = random.Random(20261104)
PROT = 'ACDEFGHIKLMNPQRSTVWY'
IUPAC = 'ACGTRYSWKMBDHVN'
HBA = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR'
HBB = 'MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRFFESFGDLSTPDAVMGNPKVKAHGKKVLGAFSDGLAHLDNLKGTFATLSELHCDKLHVDPENFRLLGNVLVCVLAHHFGKEFTPPVQAAYQKVVAGVANALAHKYH'
BPNAME = {'EBLOSUM62': 'BLOSUM62', 'EPAM250': 'PAM250', 'EDNAFULL': 'NUC.4.4'}


def mutate(s, alpha, rate):
    out = []
    for c in s:
        r = rnd.random()
        if r < rate * 0.5: out.append(rnd.choice(alpha))
        elif r < rate * 0.75: pass
        elif r < rate: out.append(c); out.append(rnd.choice(alpha))
        else: out.append(c)
    return ''.join(out) or rnd.choice(alpha)


def gen(kind):
    if kind == 'dna': alpha, mat = 'ACGT', 'EDNAFULL'
    elif kind == 'iupac': alpha, mat = IUPAC, 'EDNAFULL'
    else:
        alpha = PROT + ('BZX' if kind == 'protx' else '')
        mat = rnd.choice(['EBLOSUM62', 'EPAM250'])
    a = ''.join(rnd.choice(alpha) for _ in range(rnd.randint(8, 90)))
    return a, mutate(a, alpha, rnd.choice([0.1, 0.3, 0.6])), mat


def coords(r1, r2):
    steps = [(x != '-', y != '-') for x, y in zip(r1, r2)]
    c1, c2, i = [0], [0], 0
    while i < len(steps):
        j = i
        while j < len(steps) and steps[j] == steps[i]: j += 1
        c1.append(c1[-1] + ((j - i) if steps[i][0] else 0)); c2.append(c2[-1] + ((j - i) if steps[i][1] else 0))
        i = j
    return np.array([c1, c2])


def biopython_counts(r1, r2, mat):
    al = Alignment([r1.replace('-', ''), r2.replace('-', '')], coords(r1, r2))
    return al.counts(sm.load(BPNAME[mat]))


def run(prog, a, b, mat, go, ge, kind, d):
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind in ('dna', 'iupac') else ['-sprotein1', '-sprotein2']
    cmd = [prog, '-asequence', d + '/a.fa', '-bsequence', d + '/b.fa', '-gapopen', str(go), '-gapextend', str(ge),
           '-datafile', mat, '-awidth', '1000', '-outfile', d + '/o.txt', '-auto'] + typ
    if subprocess.run(cmd, capture_output=True, text=True).returncode != 0: return None
    t = open(d + '/o.txt').read()
    g = lambda p: int(re.search(p, t).group(1))
    lines = t.split('\n')
    for i, l in enumerate(lines):
        m = re.match(r'^a\s+\d+ ([A-Za-z\-]+)\s+\d+$', l)
        if m:
            r1 = m.group(1); r2 = re.match(r'^b\s+\d+ ([A-Za-z\-]+)\s+\d+$', lines[i + 2]).group(1)
            start = l.index(r1, 20)
            return r1, r2, lines[i + 1][start:start + len(r1)].ljust(len(r1)), g(r'# Length: (\d+)'), g(r'# Identity:\s+(\d+)/'), g(r'# Similarity:\s+(\d+)/'), g(r'# Gaps:\s+(\d+)/')
    return None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    d = tempfile.mkdtemp()
    todo = [('needle', HBA, HBB, 'EBLOSUM62', 10, 0.5, 'prot')]
    for _ in range(n):
        kind = rnd.choice(['dna', 'iupac', 'prot', 'protx', 'dna', 'prot'])
        a, b, mat = gen(kind)
        todo.append((rnd.choice(['needle', 'water']), a, b, mat, rnd.choice([10, 8, 5, 12, 3]), rnd.choice([0.5, 1, 2, 0.1]), kind))
    print('# id\tprogram\tmatrix\trow1\trow2\tmarkup\tlength\tidentity\tsimilarity\tgaps\tbp_mismatches\tbp_open_gaps')
    for k, (prog, a, b, mat, go, ge, kind) in enumerate(todo):
        r = run(prog, a, b, mat, go, ge, kind, d)
        if r is None: continue
        r1, r2, mk, L, I, S, G = r
        c = biopython_counts(r1.upper(), r2.upper(), mat)
        assert (c.identities, c.positives, c.gaps) == (I, S, G), (r1, r2)
        assert len(r1) == L
        print('\t'.join(map(str, [k, prog, mat, r1, r2, mk, L, I, S, G, c.mismatches, c.open_gaps])))

main()
