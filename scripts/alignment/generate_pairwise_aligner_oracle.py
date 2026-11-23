#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/pairwise_aligner_oracle.jsonl (B09 R14-10, F303): Bio.Align.PairwiseAligner of Biopython 1.88 (modes global / local / fogsaa, scalar and matrix scoring, the wildcard letter,
linear / affine / end-specific gap scores, strand + / -): its score, algorithm name, warnings and alignments.

One JSON object per line: {"id", "config": {"mode", "props": {name: value}, "matrix", "wildcard"}, "a", "b", "strand", "algorithm", "score", "alignments": [{"coordinates": [[...], [...]]}], "warnings": [...], "error"}.
Run: python3 scripts/alignment/generate_pairwise_aligner_oracle.py
"""
import json
import os
import random
import warnings

from Bio import Align
from Bio.Align import substitution_matrices

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'pairwise_aligner_oracle.jsonl')
rng = random.Random(20261405)


def dna(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def with_n(s):
    if len(s) > 4 and rng.random() < 0.5:
        i = rng.randrange(len(s))
        return s[:i] + 'N' + s[i + 1:]
    return s


def protein(n):
    return ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


def mutate(seq, alphabet='ACGT'):
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


def gap_props(kind):
    if kind == 'linear':
        return {'gap_score': -rng.choice([1, 2, 1.5])}
    if kind == 'affine':
        return {'open_gap_score': -rng.choice([2, 3, 4, 2.5]), 'extend_gap_score': -rng.choice([0.5, 1, 0.25])}
    if kind == 'free_ends':
        return {'open_gap_score': -rng.choice([2, 3]), 'extend_gap_score': -rng.choice([0.5, 1]), 'end_gap_score': 0.0}
    if kind == 'free_left':
        return {'gap_score': -2.0, 'left_gap_score': 0.0}
    if kind == 'free_insertion_ends':
        return {'gap_score': -1.0, 'end_insertion_score': 0.0}
    if kind == 'end_scores':
        return {'open_gap_score': -3.0, 'extend_gap_score': -1.0, 'left_gap_score': -1.5, 'right_open_gap_score': -0.5, 'right_extend_gap_score': -0.25}
    if kind == 'directional':
        return {'target_open_gap_score': -2.0, 'target_extend_gap_score': -0.5, 'query_open_gap_score': -4.0, 'query_extend_gap_score': -1.0}
    raise ValueError(kind)


def make(cfg):
    a = Align.PairwiseAligner()
    a.mode = cfg['mode']
    for k, v in cfg['props'].items():
        setattr(a, k, v)
    if cfg.get('matrix'):
        a.substitution_matrix = substitution_matrices.load(cfg['matrix'])
    if cfg.get('wildcard'):
        a.wildcard = cfg['wildcard']
    return a


records = []
kinds = ['linear', 'affine', 'free_ends', 'free_left', 'free_insertion_ends', 'end_scores', 'directional']
for i in range(700):
    mode = ['global', 'local', 'fogsaa'][i % 3]
    protein_case = rng.random() < 0.2
    gen = protein if protein_case else dna
    a_seq = gen(rng.randrange(1, 22))
    b_seq = mutate(a_seq, 'ACDEFGHIKLMNPQRSTVWY' if protein_case else 'ACGT')
    kind = rng.choice(kinds)
    if mode == 'fogsaa':
        kind = rng.choice(['linear', 'affine', 'free_ends', 'end_scores', 'directional', 'free_left'])
    cfg = {'mode': mode, 'props': gap_props(kind), 'matrix': None, 'wildcard': None, 'kind': kind}
    if protein_case:
        cfg['matrix'] = rng.choice(['BLOSUM62', 'BLOSUM80', 'PAM250'])
        cfg['props'] = {k: v * 4 for k, v in cfg['props'].items()}
    else:
        r = rng.random()
        if r < 0.55:
            cfg['props'].update({'match_score': rng.choice([1, 2, 5, 1.5]), 'mismatch_score': rng.choice([0, -1, -2, -0.5, -3])})
        elif r < 0.75:
            cfg['matrix'] = rng.choice(['NUC.4.4', 'BLASTN'])
        else:
            cfg['props'].update({'match_score': rng.choice([1, 2, 3]), 'mismatch_score': rng.choice([-1, -2])})
            cfg['wildcard'] = 'N'
            a_seq, b_seq = with_n(a_seq), with_n(b_seq)
    strand = '-' if (not protein_case and rng.random() < 0.3) else '+'
    rec = {'id': 'q%03d' % i, 'config': cfg, 'a': a_seq, 'b': b_seq, 'strand': strand}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            aligner = make(cfg)
            rec['algorithm'] = aligner.algorithm
            rec['score'] = aligner.score(a_seq, b_seq, strand=strand)
            alignments = aligner.align(a_seq, b_seq, strand=strand)
            n = len(alignments)
            if n > 400:
                continue
            rec['alignments'] = [{'coordinates': al.coordinates.tolist(), 'score': al.score} for al in alignments]
            rec['count'] = n
        except Exception as e:
            rec['error'] = type(e).__name__ + ': ' + str(e)
        rec['warnings'] = sorted({str(w.message) for w in caught if 'renamed' not in str(w.message) and 'deprecated' not in str(w.message)})
    records.append(rec)

with open(out_path, 'w', newline='\n') as fh:
    for r in records:
        fh.write(json.dumps(r) + '\n')
print(len(records), 'records ->', out_path, ', errors', sum('error' in r for r in records))
