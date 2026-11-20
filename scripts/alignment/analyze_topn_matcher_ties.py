#!/usr/bin/env python3
"""Why do 7 of the 397 EMBOSS matcher score lists of topn_oracle.tsv differ from the Waterman-Eggert lists of SequenceAligner.LocalAlignTopN? (audit A3-14c, B09 F143)

Method (everything below is measured, nothing assumed):
  * The generator (generate_topn_oracle.py) already proved that EVERY matcher round is a legal Waterman-Eggert round: its score is the optimum of the matrix with the pairs of
    matcher's own earlier rounds forbidden. So a difference of two lists can only come from a TIE: two different alignments of the same optimal score in an earlier round.
  * This script re-runs `matcher -alternatives k` on every matcher row of the committed topn_oracle.tsv, takes the aligned pairs of each matcher round, runs the independent
    full-recompute Waterman-Eggert reference of the generator (the program that `SequenceAligner.LocalAlignTopN` equals on every row, SequenceAligner_TopN_Tests) and reports, per row,
      score_round   first round whose SCORES differ (none = lists equal)
      tie_round     first round whose ALIGNMENT (set of aligned pairs) differs although the score is equal = the first tie broken differently
    A list difference is explained iff tie_round exists and tie_round <= score_round (a different alignment of an equal score in an earlier or the same round).
Usage: python3 scripts/alignment/analyze_topn_matcher_ties.py [path to topn_oracle.tsv]       (needs Biopython and `matcher` on PATH, EMBOSS 6.6.0)
Output: one line per differing row, then the counts. Run result recorded in docs/Evidence/ALIGN-LOCAL-001-Evidence.md.
"""
import os, re, subprocess, sys, tempfile

here = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(here, 'generate_topn_oracle.py')).read()
assert src.rstrip().endswith('main()')
exec(compile(src.rstrip()[:-len('main()')], 'generate_topn_oracle_nomain', 'exec'), globals())

path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, '..', '..', 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'topn_oracle.tsv')


def matcher_alignments(a, b, k, G, E, matrix, kind, d):
    open(d + '/a.fa', 'w').write('>a\n' + a + '\n'); open(d + '/b.fa', 'w').write('>b\n' + b + '\n')
    typ = ['-snucleotide1', '-snucleotide2'] if kind == 'dna' else ['-sprotein1', '-sprotein2']
    r = subprocess.run(['matcher', d + '/a.fa', d + '/b.fa', '-gapopen', str(G), '-gapextend', str(E), '-datafile', matrix, '-alternatives', str(k), '-aformat3', 'srspair',
                        '-outfile', d + '/o.txt', '-auto'] + typ, capture_output=True, text=True)
    if r.returncode != 0:
        return None
    txt = open(d + '/o.txt').read()
    out = []
    for part in re.split(r'#=======================================\n#\n# Aligned_sequences', txt)[1:]:
        sc = int(re.search(r'# Score: (\d+)', part).group(1))
        ra = rb = ''; sa = sb = None
        for line in part.split('\n'):
            m_ = re.match(r'^a\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sa = sa or int(m_.group(1)); ra += m_.group(2); continue
            m_ = re.match(r'^b\s+(\d+) (\S+)\s+(\d+)$', line)
            if m_: sb = sb or int(m_.group(1)); rb += m_.group(2)
        i, j = sa - 1, sb - 1; pairs = set()
        for x, y in zip(ra, rb):
            if x != '-': i += 1
            if y != '-': j += 1
            if x != '-' and y != '-': pairs.add((i, j))
        out.append((sc, frozenset(pairs)))
    return out


def reference_alignments(a, b, k, s, o, e, affine):
    forb = set(); out = []
    while len(out) < k:
        if affine:
            M, X, Y = dp_affine(a, b, s, o, e, forb); H = M
        else:
            H = dp_linear(a, b, s, e, forb)
        best, bi, bj = best_cell(H, len(a), len(b))
        if best <= 0: break
        r = tb_affine(a, b, M, X, Y, bi, bj, s, o, e) if affine else tb_linear(a, b, H, bi, bj, s, e, forb)
        pairs = set(); i, j = r[3], r[4]
        for x, y in zip(r[0], r[1]):
            if x != '-': i += 1
            if y != '-': j += 1
            if x != '-' and y != '-': pairs.add((i, j))
        out.append((r[2], frozenset(pairs))); forb.update(r[7])
    return out


rows = [l.rstrip('\n').split('\t') for l in open(path)][1:]
total = equal = benign = explained = 0
with tempfile.TemporaryDirectory() as d:
    for idx, f in enumerate(rows):
        if f[0].startswith('#') or f[12] == 'NA':
            continue
        kind_s, model, matrix = f[0], f[1], f[2]
        a, b, k = f[7], f[8], len(f[12].split(','))
        affine = model == 'affine'
        if kind_s == 'scalar':
            ma, mi, go, ge = int(f[3]), int(f[4]), float(f[5]), float(f[6]); s = lambda x, y, ma=ma, mi=mi: ma if x == y else mi; scale = 1
            o_s, e_s = int(round(go)) if affine else 0, int(round(ge)); G, E = int(-(go + ge)), int(-ge)
            kind, mname = 'dna', 'EDNAFULL'
        else:
            scale = 4; go, ge = float(f[5]), float(f[6]); s = lambda x, y: BTAB[(x, y)] * 4
            o_s, e_s = int(round(go * scale)), int(round(ge * scale)); G, E = int(-(go + ge)), int(-ge); kind, mname = 'prot', 'EBLOSUM62'
        mt = matcher_alignments(a, b, k, G, E, mname, kind, d)
        ref = reference_alignments(a, b, k, s, o_s, e_s, affine)
        total += 1
        ms = [x[0] for x in mt]; rs = [x[0] // scale for x in ref]
        score_round = next((r for r in range(min(len(ms), len(rs))) if ms[r] != rs[r]), None)
        if score_round is None and len(ms) != len(rs):
            score_round = min(len(ms), len(rs))
        tie_round = next((r for r in range(min(len(mt), len(ref))) if mt[r][1] != ref[r][1] and mt[r][0] == ref[r][0] // scale), None)
        if score_round is None:
            equal += 1
            if tie_round is not None: benign += 1
            continue
        ok = tie_round is not None and tie_round <= score_round
        explained += ok
        print(f'line {idx + 2}: {kind_s}/{model} {a} {b} k={k}  matcher={ms} reference={rs}  first score difference at round {score_round + 1}, first tie broken differently at round {None if tie_round is None else tie_round + 1}  {"EXPLAINED (tie)" if ok else "NOT EXPLAINED"}')
print(f'rows with a matcher list: {total}; lists equal: {equal} (of which {benign} broke a tie differently without changing a score); lists differ: {total - equal}, explained by a tie broken differently in an earlier or the same round: {explained}')
