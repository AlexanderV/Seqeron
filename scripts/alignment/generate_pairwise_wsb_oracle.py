#!/usr/bin/env python3
"""B09 WP67 (F307, R14-10b): oracle of Bio.Align.PairwiseAligner 1.88 for the settings the library's own aligners refused: gap-score functions of (position, length) (algorithm 'Waterman-Smith-Beyer global / local alignment algorithm':
insertion_score / deletion_score / gap_score, strand '+' and '-') and gap scores above zero (open / extend / internal / left / right, algorithms Needleman-Wunsch, Smith-Waterman, Gotoh).
JSONL.GZ rows: id, a, b, mode, strand, match, mismatch, cfg (the settings: gap functions as {kind, params}, gap scores as {property: value}), algorithm, score, alignments (coordinates of every alignment, at most 400), total (the count; larger than the list when truncated), error.
Usage: python3 scripts/alignment/generate_pairwise_wsb_oracle.py OUT.jsonl.gz [n=900] [seed]      (needs Biopython 1.88)"""
import gzip, json, math, random, sys
from Bio import Align

FAMILIES = ['lin', 'log', 'sq', 'pos', 'step', 'tab', 'zero', 'bonus']


def make(kind, p):
    """the gap-score function of a family (the C# test builds the same one)"""
    if kind == 'lin': return lambda i, l: -(p[0] + p[1] * l)
    if kind == 'log': return lambda i, l: -(p[0] + p[1] * math.log(l + 1))
    if kind == 'sq': return lambda i, l: -(p[0] + p[1] * l * l)
    if kind == 'pos': return lambda i, l: -(p[0] + p[1] * i + p[2] * l)
    if kind == 'step': return lambda i, l: -p[0] if l <= int(p[2]) else -p[1]
    if kind == 'tab': return lambda i, l: -(((i * 7 + l * 3) % int(p[1])) * 0.5 + p[0])
    if kind == 'zero': return lambda i, l: 0.0
    if kind == 'bonus': return lambda i, l: p[0] * l - p[1]            # a gap of positive length scores above zero
    raise ValueError(kind)


SPECIAL = [
    {'a': 'CACTC', 'b': 'GTCTG', 'mode': 'local', 'strand': '+', 'match': 1, 'mismatch': -1, 'cfg': {'scores': {'gap_score': 0.5}}},
]


def rand_func(rnd):
    kind = rnd.choice(FAMILIES)
    p = [rnd.choice([0, 0.5, 1, 1.5, 2, 3]), rnd.choice([0.5, 1, 1.5, 2]), rnd.choice([1, 2, 3, 4])]
    if kind == 'tab': p[1] = rnd.choice([2, 3, 4, 5])
    return {'kind': kind, 'p': p}


def seq(rnd, alphabet, lo=1, hi=8):
    return ''.join(rnd.choice(alphabet) for _ in range(rnd.randint(lo, hi)))


def coords(a):
    return [[int(x) for x in row] for row in a.coordinates.tolist()]


def run(row):
    al = Align.PairwiseAligner()
    al.mode = row['mode']
    al.match_score = row['match']
    al.mismatch_score = row['mismatch']
    cfg = row['cfg']
    if 'gap' in cfg:
        f = make(cfg['gap']['kind'], cfg['gap']['p'])
        al.gap_score = f
    if 'ins' in cfg: al.insertion_score = make(cfg['ins']['kind'], cfg['ins']['p'])
    if 'del' in cfg: al.deletion_score = make(cfg['del']['kind'], cfg['del']['p'])
    for k, v in cfg.get('scores', {}).items():
        setattr(al, k, v)
    row['algorithm'] = al.algorithm
    try:
        row['score'] = al.score(row['a'], row['b'], strand=row['strand'])
        als = al.align(row['a'], row['b'], strand=row['strand'])
        row['total'] = len(als)
        row['alignments'] = []
        for k, a in enumerate(als):
            if k >= 400: break
            row['alignments'].append(coords(a))
        row['aligned_score'] = als.score
    except Exception as e:
        row['error'] = '%s: %s' % (type(e).__name__, e)
    return row


def main():
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 900
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20262401
    rnd = random.Random(seed)
    rows = []
    while len(rows) < n:
        kind = rnd.random()
        alphabet = rnd.choice(['ACGT', 'ACGT', 'AC', 'ACDE'])
        row = {'id': len(rows), 'a': seq(rnd, alphabet), 'b': seq(rnd, alphabet), 'mode': rnd.choice(['global', 'local']), 'strand': rnd.choice(['+', '+', '+', '-']) if alphabet == 'ACGT' else '+',
               'match': rnd.choice([1, 2, 3]), 'mismatch': rnd.choice([0, -1, -2]), 'cfg': {}}
        if kind < 0.6:
            r = rnd.random()
            if r < 0.4: row['cfg']['gap'] = rand_func(rnd)
            elif r < 0.7: row['cfg']['ins'] = rand_func(rnd)
            elif r < 0.9: row['cfg']['del'] = rand_func(rnd)
            else: row['cfg']['ins'], row['cfg']['del'] = rand_func(rnd), rand_func(rnd)
            if rnd.random() < 0.25:
                # a number set on the other direction after the function
                row['cfg']['scores'] = {rnd.choice(['open_gap_score', 'extend_gap_score']): -rnd.choice([1, 2])} if 'gap' not in row['cfg'] else {}
                if row['cfg']['scores'] and ('ins' in row['cfg'] and 'del' in row['cfg']): row['cfg']['scores'] = {}
        else:
            # gap scores above zero: all combinations of the twelve scores in a small set
            sc = {}
            names = rnd.choice([['open_gap_score', 'extend_gap_score'], ['open_gap_score'], ['gap_score'], ['open_internal_gap_score', 'extend_internal_gap_score', 'open_end_gap_score', 'extend_end_gap_score'], ['open_left_insertion_score', 'open_right_deletion_score', 'extend_gap_score']])
            for nm in names:
                sc[nm] = rnd.choice([0.5, 1, 2, -1, -2, 0, 1.5])
            if not any(v > 0 for v in sc.values()): sc[names[0]] = 1
            row['cfg']['scores'] = sc
            row['a'] = row['a'][:5]; row['b'] = row['b'][:5]
        rows.append(run(row))
    # rows on which the reference itself fails (its traceback meets a cell without a move: RuntimeError)
    for spec in SPECIAL:
        rows.append(run(dict(spec, id=len(rows))))
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    from collections import Counter
    print(len(rows), 'rows;', Counter(r['algorithm'] for r in rows), 'errors', sum(1 for r in rows if 'error' in r), 'truncated', sum(1 for r in rows if r.get('total', 0) > 400))


if __name__ == '__main__':
    main()
