#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/parasail_oracle.jsonl (B09 R14-11, F304): the score-only API of parasail 1.3.4 (nw, sw and the sg variants; plain, _stats, _table, _rowcol, _stats_table, _stats_rowcol)
on random DNA and protein pairs, with the score, end_query, end_ref, matches, similar, length and the tables / rows / columns.

One JSON object per line: {"id", "fn" (nw / sw / sg_*), "s1", "s2", "open", "gap", "matrix" (name: "dna:<match>:<mismatch>" or "blosum62"), "results": {variant: {...}}}.
Run: python3 scripts/alignment/generate_parasail_oracle.py
"""
import json
import os
import random

import parasail

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'parasail_oracle.jsonl')
rng = random.Random(20261406)
VARIANTS = ['', '_stats', '_table', '_rowcol', '_stats_table', '_stats_rowcol']
FUNCTIONS = ['nw', 'sw', 'sg', 'sg_qb', 'sg_qe', 'sg_db', 'sg_de', 'sg_qx', 'sg_dx', 'sg_qb_db', 'sg_qb_de', 'sg_qe_db', 'sg_qe_de']
AA = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(seq, alphabet):
    out = []
    for ch in seq:
        r = rng.random()
        if r < 0.15:
            out.append(rng.choice(alphabet))
        elif r < 0.22:
            continue
        elif r < 0.30:
            out.append(ch + rng.choice(alphabet))
        else:
            out.append(ch)
    return ''.join(out) or alphabet[0]


def grab(result, variant):
    d = {'score': int(result.score), 'end_query': int(result.end_query), 'end_ref': int(result.end_ref)}
    if 'stats' in variant:
        d.update({'matches': int(result.matches), 'similar': int(result.similar), 'length': int(result.length)})
    if 'table' in variant:
        d['score_table'] = result.score_table.tolist()
        if 'stats' in variant:
            d['matches_table'] = result.matches_table.tolist()
            d['similar_table'] = result.similar_table.tolist()
            d['length_table'] = result.length_table.tolist()
    if 'rowcol' in variant:
        d['score_row'] = result.score_row.tolist()
        d['score_col'] = result.score_col.tolist()
        if 'stats' in variant:
            for k in ('matches', 'similar', 'length'):
                d[k + '_row'] = getattr(result, k + '_row').tolist()
                d[k + '_col'] = getattr(result, k + '_col').tolist()
    return d


records = []
for i in range(520):
    fn = FUNCTIONS[i % len(FUNCTIONS)]
    protein = i % 3 == 2
    alphabet = AA if protein else 'ACGT'
    n = rng.randrange(1, 28)
    s1 = ''.join(rng.choice(alphabet) for _ in range(n))
    s2 = mutate(s1, alphabet) if rng.random() < 0.8 else ''.join(rng.choice(alphabet) for _ in range(rng.randrange(1, 30)))
    if rng.random() < 0.2:
        s2 = s2 + ''.join(rng.choice(alphabet) for _ in range(rng.randrange(0, 12)))
    if rng.random() < 0.2:
        s1, s2 = s2, s1
    if protein:
        matrix, mname = parasail.blosum62, 'blosum62'
        open_, gap = rng.choice([(11, 1), (10, 2), (8, 3), (5, 5)])
    else:
        match, mismatch = rng.choice([(2, -3), (1, -1), (5, -4), (3, -2)])
        matrix, mname = parasail.matrix_create('ACGT', match, mismatch), 'dna:%d:%d' % (match, mismatch)
        open_, gap = rng.choice([(5, 2), (3, 1), (2, 2), (6, 1), (4, 3)])
    rec = {'id': 'p%03d' % i, 'fn': fn, 's1': s1, 's2': s2, 'open': open_, 'gap': gap, 'matrix': mname, 'results': {}}
    for v in VARIANTS:
        f = getattr(parasail, fn + v, None)
        if f is None:
            continue
        rec['results'][v or 'plain'] = grab(f(s1, s2, open_, gap, matrix), v)
    records.append(rec)

with open(out_path, 'w', newline='\n') as fh:
    for r in records:
        fh.write(json.dumps(r) + '\n')
print(len(records), 'records ->', out_path, ', variants per record', sorted({len(r['results']) for r in records}))
