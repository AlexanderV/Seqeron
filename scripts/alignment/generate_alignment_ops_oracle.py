#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/alignment_ops_oracle.jsonl (B09 R14-10, F303): the Bio.Align.Alignment operations of Biopython 1.88 on random alignments.

One JSON object per line: {"id", "kind", "spec": {"sequences": [{"id", "seq", "weight"}], "coordinates": [[...]]}, "ops": {...}} where ops holds the value of every operation
(shape, length, rows, aligned, indices, inverse_indices, frequencies, substitutions, counts, reverse_complement, sort, str, slices, add, comparison ...) or {"error": "..."} when
Biopython raises. Second alignments (add, map) are stored in "spec2"; mapall / from_alignments_with_same_reference are separate records (kind = mapall / same_reference / sorting).
Run: python3 scripts/alignment/generate_alignment_ops_oracle.py
"""
import json
import os
import random
import warnings

import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq, reverse_complement
from Bio.SeqRecord import SeqRecord

warnings.simplefilter('ignore')
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'alignment_ops_oracle.jsonl')
rng = random.Random(20261404)


def dna(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def protein(n):
    return ''.join(rng.choice('ACDEFGHIKLMNPQRSTVWY') for _ in range(n))


def mutate(seq, alphabet='ACGT'):
    out = []
    for ch in seq:
        r = rng.random()
        if r < 0.12:
            out.append(rng.choice(alphabet))
        elif r < 0.18:
            continue
        elif r < 0.24:
            out.append(ch + rng.choice(alphabet))
        else:
            out.append(ch)
    return ''.join(out) or alphabet[0]


def safe(fn):
    try:
        return fn()
    except Exception as e:
        return {'error': type(e).__name__}


def to_list(x):
    if isinstance(x, np.ndarray):
        return x.tolist()
    return x


def spec_of(alignment, records=None):
    seqs = []
    for i, s in enumerate(alignment.sequences):
        if isinstance(s, SeqRecord):
            seqs.append({'id': s.id, 'seq': str(s.seq), 'weight': s.annotations.get('weight')})
        else:
            seqs.append({'id': None, 'seq': str(s), 'weight': None})
    return {'sequences': seqs, 'coordinates': alignment.coordinates.tolist()}


def make(spec):
    seqs = []
    for s in spec['sequences']:
        if s['id'] is not None:
            r = SeqRecord(Seq(s['seq']), id=s['id'])
            if s['weight'] is not None:
                r.annotations['weight'] = s['weight']
            seqs.append(r)
        else:
            seqs.append(s['seq'])
    return Align.Alignment(seqs, np.array(spec['coordinates'], dtype=np.intp))


def random_pairwise(alphabet='ACGT', strand='+'):
    gen = dna if alphabet == 'ACGT' else protein
    a = gen(rng.randrange(5, 25))
    b = mutate(a, alphabet)
    aligner = Align.PairwiseAligner()
    aligner.match_score, aligner.mismatch_score, aligner.gap_score = 2, -1, -2
    aligner.mode = rng.choice(['global', 'local'])
    if alphabet != 'ACGT':
        strand = '+'
    alignments = aligner.align(a, b, strand=strand)
    al = alignments[0]
    return al


def random_msa(n, alphabet='ACGT', rc_rows=()):
    gen = dna if alphabet == 'ACGT' else protein
    cols = rng.randrange(6, 20)
    rows = ['' for _ in range(n)]
    for _ in range(cols):
        present = [rng.random() < 0.8 for _ in range(n)]
        if sum(present) == 0:
            present[rng.randrange(n)] = True
        for i in range(n):
            rows[i] += rng.choice(alphabet) if present[i] else '-'
    lines = [r.encode() for r in rows]
    sequences, coordinates = Align.Alignment.parse_printed_alignment(lines)
    sequences = [s.decode() for s in sequences]
    coordinates = coordinates.copy()
    for i in rc_rows:
        length = len(sequences[i])
        sequences[i] = reverse_complement(sequences[i])
        coordinates[i, :] = length - coordinates[i, :]
    return Align.Alignment(sequences, coordinates)


def with_records(al):
    """Same alignment with SeqRecord sequences (ids s0.., weights)."""
    seqs = []
    for i, s in enumerate(al.sequences):
        r = SeqRecord(Seq(str(s)), id='s%d' % i)
        if rng.random() < 0.5:
            r.annotations['weight'] = rng.choice([0.5, 1.0, 2.0])
        seqs.append(r)
    return Align.Alignment(seqs, al.coordinates.copy())


COUNT_NAMES = ['open_left_insertions', 'extend_left_insertions', 'open_left_deletions', 'extend_left_deletions',
               'open_internal_insertions', 'extend_internal_insertions', 'open_internal_deletions', 'extend_internal_deletions',
               'open_right_insertions', 'extend_right_insertions', 'open_right_deletions', 'extend_right_deletions',
               'aligned', 'identities', 'mismatches', 'positives', 'gap_score', 'substitution_score', 'score',
               'left_gaps', 'internal_gaps', 'right_gaps', 'gaps', 'left_insertions', 'left_deletions',
               'internal_insertions', 'internal_deletions', 'right_insertions', 'right_deletions', 'insertions', 'deletions',
               'open_gaps', 'extend_gaps']

AMBIGUOUS = {'wildcard': 'N'}


def counts_dict(c):
    d = {}
    for name in COUNT_NAMES:
        try:
            v = getattr(c, name)
        except AttributeError:
            continue
        d[name] = v
    d['str'] = str(c)
    return d


def aligner_configs(protein_case):
    configs = [
        {'name': 'default', 'props': {}},
        {'name': 'scores', 'props': {'match_score': 2, 'mismatch_score': -1, 'open_gap_score': -3.0, 'extend_gap_score': -0.5}},
        {'name': 'ends', 'props': {'match_score': 1, 'mismatch_score': -1, 'open_gap_score': -2.0, 'extend_gap_score': -1.0, 'left_gap_score': -0.25, 'right_extend_gap_score': -0.75, 'right_open_gap_score': -1.5, 'target_gap_score': -1.25}},
        {'name': 'function', 'props': {'match_score': 1, 'mismatch_score': -1}, 'function': True},
    ]
    if protein_case:
        configs.append({'name': 'blosum', 'props': {'open_gap_score': -10.0, 'extend_gap_score': -0.5}, 'matrix': 'BLOSUM62'})
    else:
        configs.append({'name': 'wildcard', 'props': {'match_score': 2, 'mismatch_score': -1, 'open_gap_score': -2.0, 'extend_gap_score': -1.0}, 'wildcard': 'N'})
        configs.append({'name': 'nuc44', 'props': {'open_gap_score': -3.0, 'extend_gap_score': -1.0}, 'matrix': 'NUC.4.4'})
    return configs


def make_aligner(cfg):
    a = Align.PairwiseAligner()
    for k, v in cfg['props'].items():
        setattr(a, k, v)
    if 'matrix' in cfg:
        a.substitution_matrix = substitution_matrices.load(cfg['matrix'])
    if 'wildcard' in cfg:
        a.wildcard = cfg['wildcard']
    if cfg.get('function'):
        a.insertion_score = lambda position, length: -1.0 - 0.5 * length - 0.01 * position
        a.deletion_score = lambda position, length: -2.0 - 0.25 * length
    return a


def counts_all(al, protein_case):
    result = {}
    result['plain'] = safe(lambda: counts_dict(al.counts()))
    result['wildcard'] = safe(lambda: counts_dict(al.counts('N')))
    mname = 'BLOSUM62' if protein_case else 'NUC.4.4'
    result['matrix'] = safe(lambda: counts_dict(al.counts(substitution_matrices.load(mname))))
    result['aligners'] = {}
    for cfg in aligner_configs(protein_case):
        result['aligners'][cfg['name']] = safe(lambda cfg=cfg: counts_dict(al.counts(make_aligner(cfg))))
    result['aligner_configs'] = aligner_configs(protein_case)
    return result


def slices_of(al):
    out = []
    n, m = al.shape
    for _ in range(12):
        r1 = rng.randrange(-n, n)
        c1 = rng.randrange(-m, m)
        rs = [rng.choice([None, rng.randrange(-n, n + 1)]) for _ in range(2)]
        step = rng.choice([None, None, 1, 2, -1])
        cs = [rng.choice([None, rng.randrange(-m, m + 1)]) for _ in range(2)]
        cstep = rng.choice([None, None, None, 2, -1])
        cols = [rng.randrange(-m, m) for _ in range(rng.randrange(1, 5))]
        row_slice = slice(rs[0], rs[1], step)
        col_slice = slice(cs[0], cs[1], cstep)
        entry = {'row': r1, 'col': c1, 'rowsl': [rs[0], rs[1], step], 'colsl': [cs[0], cs[1], cstep], 'collist': cols}
        entry['row_col'] = safe(lambda: al[r1, c1])
        entry['row_cols'] = safe(lambda: al[r1, col_slice])
        entry['row_collist'] = safe(lambda: al[r1, cols])
        entry['rows_col'] = safe(lambda: al[row_slice, c1])
        entry['rows'] = safe(lambda: spec_of(al[row_slice]))
        entry['rows_cols'] = safe(lambda: spec_of(al[row_slice, col_slice]))
        entry['rows_collist'] = safe(lambda: spec_of(al[row_slice, cols]))
        entry['row_int'] = safe(lambda: al[r1])
        out.append(entry)
    return out


def simple_spec(alignment_spec):
    return alignment_spec


def main():
    records = []
    cid = 0
    for i in range(260):
        kind = ['pairwise', 'pairwise_rc', 'msa', 'msa_rc', 'protein_pairwise', 'protein_msa'][i % 6]
        protein_case = kind.startswith('protein')
        if kind == 'pairwise':
            al = random_pairwise('ACGT', '+')
        elif kind == 'pairwise_rc':
            al = random_pairwise('ACGT', '-')
        elif kind == 'protein_pairwise':
            al = random_pairwise('ACDEFGHIKLMNPQRSTVWY')
        elif kind == 'protein_msa':
            al = random_msa(rng.randrange(2, 5), 'ACDEFGHIKLMNPQRSTVWY')
        elif kind == 'msa':
            al = random_msa(rng.randrange(2, 5))
        else:
            n = rng.randrange(2, 5)
            al = random_msa(n, rc_rows=[rng.randrange(n)])
        if i % 3 == 2:
            al = with_records(al)
        spec = spec_of(al)
        ops = {}
        ops['shape'] = safe(lambda: list(al.shape))
        ops['length'] = safe(lambda: al.length)
        ops['rows'] = safe(lambda: [al[k] for k in range(len(al))])
        ops['array'] = safe(lambda: [''.join(c.decode() for c in row) for row in np.array(al)])
        ops['aligned'] = safe(lambda: al.aligned.tolist())
        ops['indices'] = safe(lambda: al.indices.tolist())
        ops['inverse_indices'] = safe(lambda: [x.tolist() for x in al.inverse_indices])
        ops['frequencies'] = safe(lambda: {k: v.tolist() for k, v in al.frequencies.items()})
        ops['substitutions'] = safe(lambda: {'letters': ''.join(al.substitutions.alphabet), 'values': al.substitutions.tolist()})
        ops['reverse_complement'] = safe(lambda: spec_of(al.reverse_complement()))
        ops['str'] = safe(lambda: str(al))
        mname = 'BLOSUM62' if protein_case else 'NUC.4.4'
        ops['str_matrix'] = safe(lambda: al.format('', substitution_matrices.load(mname)) if False else al.format(scoring=substitution_matrices.load(mname)))
        ops['counts'] = counts_all(al, protein_case)

        def sorted_alignment(key_reverse):
            c = Align.Alignment(list(al.sequences), al.coordinates.copy())
            c.sort(reverse=key_reverse)
            return spec_of(c)
        ops['sort'] = safe(lambda: sorted_alignment(False))
        ops['sort_reverse'] = safe(lambda: sorted_alignment(True))
        ops['slices'] = slices_of(al)
        # add: the same sequences, a second part (same number of rows)
        if not kind.endswith('rc'):
            other = random_msa(len(al), 'ACDEFGHIKLMNPQRSTVWY' if protein_case else 'ACGT')
            ops['spec2'] = spec_of(other)
            ops['add'] = safe(lambda: spec_of(al + other))
            ops['compare'] = {'lt': safe(lambda: bool(al < other)), 'le': safe(lambda: bool(al <= other)), 'gt': safe(lambda: bool(al > other)),
                              'ge': safe(lambda: bool(al >= other)), 'eq': safe(lambda: bool(al == other)), 'eq_self': bool(al == al)}
        records.append({'id': 'a%03d' % cid, 'kind': kind, 'spec': spec, 'ops': ops})
        cid += 1

    # map: alignments of target-query and query-target', both strands
    for i in range(120):
        strand1 = rng.choice(['+', '+', '-'])
        strand2 = rng.choice(['+', '+', '-'])
        t = dna(rng.randrange(20, 50))
        q = mutate(t)
        aligner = Align.PairwiseAligner()
        aligner.match_score, aligner.mismatch_score, aligner.gap_score = 2, -1, -2
        aligner.mode = rng.choice(['global', 'local'])
        a1 = aligner.align(t, q, strand=strand1)[0]
        qq = q if strand1 == '+' else q
        r = mutate(q)
        aligner.mode = rng.choice(['global', 'local'])
        a2 = aligner.align(q, r, strand=strand2)[0]
        rec = {'id': 'm%03d' % i, 'kind': 'map', 'spec': spec_of(a1), 'spec2': spec_of(a2), 'strands': [strand1, strand2]}
        rec['map'] = safe(lambda: spec_of(a1.map(a2)))
        records.append(rec)

    # mapall (factor 1) and with protein / DNA alignments (factor 3)
    for i in range(40):
        n = rng.randrange(2, 5)
        msa = random_msa(n, 'ACGT')
        lengths = [len(s) for s in msa.sequences]
        alignments = []
        for k in range(n):
            genome = dna(lengths[k] + rng.randrange(0, 6))
            aligner = Align.PairwiseAligner()
            aligner.match_score, aligner.mismatch_score, aligner.gap_score = 2, -1, -2
            aligner.mode = 'local'
            seq = str(msa.sequences[k])
            # a transcript of the sequence inside a longer genome sequence
            genome = dna(rng.randrange(0, 5)) + seq + dna(rng.randrange(0, 5))
            al = aligner.align(genome, seq)[0]
            alignments.append(al)
        rec = {'id': 'ma%03d' % i, 'kind': 'mapall', 'msa': spec_of(msa), 'alignments': [spec_of(a) for a in alignments]}
        rec['mapall'] = safe(lambda: spec_of(msa.mapall([a for a in alignments])))
        records.append(rec)

    # from_alignments_with_same_reference
    for i in range(60):
        ref = dna(rng.randrange(5, 20))
        alignments = []
        for k in range(rng.randrange(1, 4)):
            q = mutate(ref)
            aligner = Align.PairwiseAligner()
            aligner.match_score, aligner.mismatch_score, aligner.gap_score = 2, -1, -2
            aligner.mode = 'global'
            alignments.append(aligner.align(ref, q)[0])
        rec = {'id': 'r%03d' % i, 'kind': 'same_reference', 'alignments': [spec_of(a) for a in alignments]}
        rec['result'] = safe(lambda: spec_of(Align.Alignment.from_alignments_with_same_reference(alignments)))
        records.append(rec)

    # error cases of from_alignments_with_same_reference
    a = Align.PairwiseAligner().align('ACGT', 'ACT')[0]
    b = Align.PairwiseAligner().align('ACGA', 'ACT')[0]
    records.append({'id': 'r_bad1', 'kind': 'same_reference', 'alignments': [spec_of(a), spec_of(b)], 'result': safe(lambda: spec_of(Align.Alignment.from_alignments_with_same_reference([a, b])))})
    records.append({'id': 'r_bad2', 'kind': 'same_reference', 'alignments': [], 'result': safe(lambda: spec_of(Align.Alignment.from_alignments_with_same_reference([])))})

    # parse_printed_alignment
    for i in range(60):
        n = rng.randrange(1, 5)
        cols = rng.randrange(1, 15)
        rows = [''.join(rng.choice('ACGT-') for _ in range(cols)) for _ in range(n)]
        res = safe(lambda: [[s.decode() for s in Align.Alignment.parse_printed_alignment([r.encode() for r in rows])[0]],
                            Align.Alignment.parse_printed_alignment([r.encode() for r in rows])[1].tolist()])
        records.append({'id': 'p%03d' % i, 'kind': 'parse_printed', 'rows': rows, 'result': res})

    # sorted(alignments)
    for i in range(30):
        pool = []
        base = dna(rng.randrange(6, 14))
        for k in range(rng.randrange(3, 7)):
            q = mutate(base)
            aligner = Align.PairwiseAligner()
            aligner.match_score, aligner.mismatch_score, aligner.gap_score = 1, -1, -1
            aligner.mode = rng.choice(['global', 'local'])
            als = aligner.align(base, q)
            pool.extend(list(als)[:2])
        order = sorted(range(len(pool)), key=lambda k: pool[k]) if False else [pool.index(x) for x in sorted(pool)]
        records.append({'id': 's%03d' % i, 'kind': 'sorting', 'alignments': [spec_of(a) for a in pool], 'order': order})

    with open(out_path, 'w', newline='\n') as fh:
        for r in records:
            fh.write(json.dumps(r) + '\n')
    print(len(records), 'records ->', out_path)


main()
