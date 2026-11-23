#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/alignment_result_bridge_oracle.tsv (B09 R14-9, F302): the coordinates and the Bio.Align text formats of the
alignments Biopython 1.88 itself finds (PairwiseAligner: global, fitting, local, overlap, affine), described by the ROWS of the alignment the way an
Seqeron AlignmentResult carries them, so that the C# bridge (AlignmentResultCoordinates) is compared with Biopython on the same alignment.

Columns: id, mode, tname, qname, tseq, qseq, row1, row2, start1, start2, score, coordinates (JSON), formats (JSON: name -> text or {"error": ...}), trimmed (JSON: the same
for alignment[:, first:last] = the leading / trailing gap columns dropped; null when the alignment has no aligned column).
Escaping of the text columns: backslash -> \\\\, newline -> \\n, tab -> \\t, carriage return -> \\r.
Run: python3 scripts/alignment/generate_alignment_result_bridge_oracle.py
"""
import json
import os
import random
import warnings

import numpy as np
from Bio import Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

warnings.simplefilter('ignore')
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'alignment_result_bridge_oracle.tsv')
rng = random.Random(20261402)

FORMATS = [('bed', {}), ('bed3', {}), ('bed6', {}), ('psl', {}), ('chain', {}), ('maf', {}), ('sam', {}), ('sam:md', {}), ('exonerate', {}), ('exonerate:cigar', {})]


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def fmt(alignment, name):
    try:
        if name.startswith('bed'):
            n = int(name[3:]) if len(name) > 3 else 12
            return alignment.format('bed', bedN=n)
        if name == 'sam:md':
            return alignment.format('sam', md=True)
        if name == 'exonerate:cigar':
            return alignment.format('exonerate', 'cigar')
        return alignment.format(name)
    except Exception as e:  # Biopython refuses: the bridge must refuse as well
        return {'error': type(e).__name__ + ': ' + str(e)}


def rows_of(alignment):
    """The two gapped rows of a pairwise alignment (target, query) from its coordinates."""
    t = str(alignment.target.seq if hasattr(alignment.target, 'seq') else alignment.target)
    q = str(alignment.query.seq if hasattr(alignment.query, 'seq') else alignment.query)
    c = alignment.coordinates
    r1, r2 = '', ''
    for i in range(1, c.shape[1]):
        dt, dq = int(c[0, i] - c[0, i - 1]), int(c[1, i] - c[1, i - 1])
        r1 += t[c[0, i - 1]:c[0, i]] if dt else '-' * dq
        r2 += q[c[1, i - 1]:c[1, i]] if dq else '-' * dt
    return r1, r2


def trimmed(alignment):
    c = alignment.coordinates
    n = c.shape[1] - 1
    first, last = 0, n
    while first < last and (c[0, first + 1] == c[0, first] or c[1, first + 1] == c[1, first]):
        first += 1
    while last > first and (c[0, last] == c[0, last - 1] or c[1, last] == c[1, last - 1]):
        last -= 1
    if first >= last:
        return None
    # the columns of the blocks first .. last-1
    ncols = [int(max(c[0, i + 1] - c[0, i], c[1, i + 1] - c[1, i])) for i in range(n)]
    lo, hi = sum(ncols[:first]), sum(ncols[:last])
    cut = alignment[:, lo:hi]
    cut.score = alignment.score
    return cut


def random_dna(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def mutate(seq):
    out = []
    for ch in seq:
        r = rng.random()
        if r < 0.10:
            out.append(rng.choice('ACGT'))
        elif r < 0.15:
            continue
        elif r < 0.20:
            out.append(ch + rng.choice('ACGT'))
        else:
            out.append(ch)
    return ''.join(out)


def aligner_for(mode):
    a = Align.PairwiseAligner()
    a.match_score, a.mismatch_score = 2, -1
    if mode in ('global_affine', 'fitting_affine', 'local_affine'):
        a.open_gap_score, a.extend_gap_score = -3, -1
    else:
        a.gap_score = -2
    if mode.startswith('local'):
        a.mode = 'local'
    elif mode.startswith('fitting'):
        a.end_deletion_score = 0
    elif mode.startswith('overlap'):
        a.end_deletion_score = 0
        a.end_insertion_score = 0
    return a


def main():
    rows = []
    modes = ['global', 'global_affine', 'local', 'local_affine', 'fitting', 'fitting_affine', 'overlap']
    case = 0
    for i in range(210):
        mode = modes[i % len(modes)]
        base = random_dna(rng.randrange(6, 40))
        t = random_dna(rng.randrange(0, 8)) + base + random_dna(rng.randrange(0, 8))
        q = mutate(base) if mode != 'fitting' and mode != 'fitting_affine' else mutate(base)
        if not q:
            q = 'A'
        if rng.random() < 0.15:
            t, q = q, t
        if rng.random() < 0.1:
            t = t[:3] + 'N' + t[4:] if len(t) > 5 else t
        a = aligner_for(mode)
        alignments = a.align(Seq(t), Seq(q))
        al0 = alignments[0]
        # the same alignment with records (identifiers for the SAM / MAF / PSL writers) and the score as Biopython sets it
        al = Align.Alignment([SeqRecord(Seq(t), id='ref%d' % i), SeqRecord(Seq(q), id='read%d' % i)], al0.coordinates)
        al.score = al0.score
        if al.coordinates.shape[1] < 2:
            continue
        r1, r2 = rows_of(al)
        c = al.coordinates
        if mode.startswith('local'):
            start1, start2 = int(c[0, 0]), int(c[1, 0])
        else:
            # the first index of each sequence consumed by the path: the residues of the leading columns up to the first aligned pair
            start1 = start2 = 0
            for x, y in zip(r1, r2):
                if x != '-' and y != '-':
                    break
                start1 += x != '-'
                start2 += y != '-'
            if mode.startswith('global'):
                start1 = start2 = 0
        formats = {name: fmt(al, name) for name, _ in FORMATS}
        tr = trimmed(al)
        tr_formats = None if tr is None else {name: fmt(tr, name) for name, _ in FORMATS}
        rows.append('\t'.join([
            'b%03d' % case, mode, 'ref%d' % i, 'read%d' % i, t, q, r1, r2, str(start1), str(start2), repr(float(al.score)),
            esc(json.dumps(np.asarray(c).tolist())), esc(json.dumps(formats)), esc(json.dumps(tr_formats))]))
        case += 1
    with open(out_path, 'w', newline='\n') as f:
        f.write('id\tmode\ttname\tqname\ttseq\tqseq\trow1\trow2\tstart1\tstart2\tscore\tcoordinates\tformats\ttrimmed\n' + '\n'.join(rows) + '\n')
    print(len(rows), 'cases ->', out_path)


main()
