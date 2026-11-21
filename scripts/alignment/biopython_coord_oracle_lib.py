"""Shared helpers of the generators of the Bio.Align coordinate-format oracles (B09 A2-27c): random pairwise alignments, the JSON description of an Alignment, and its reconstruction.

An alignment is described as {"sequences": [{"id", "length", "parts": {start: text} | null, "description"?, "annotations"?}], "coordinates": [[...], ...], "attrs": {...}}.
`parts` is the defined part of the sequence (one entry at 0 for a complete sequence, none = `Seq(None, length)`).
"""
import json
import random

import numpy as np
from Bio.Align import Alignment
from Bio.Seq import Seq, UndefinedSequenceError
from Bio.SeqRecord import SeqRecord

COMP = str.maketrans('ACGTacgtNn', 'TGCAtgcaNn')


def revcomp(s):
    return s.translate(COMP)[::-1]


def make_seq(part):
    if part is None:
        return None
    return part


def record_from(spec):
    parts = spec.get('parts')
    length = spec['length']
    if parts is None:
        seq = Seq(None, length=length)
    elif len(parts) == 1 and '0' in parts and len(parts['0']) == length:
        seq = Seq(parts['0'])
    else:
        seq = Seq({int(k): v for k, v in parts.items()}, length=length)
    rec = SeqRecord(seq, id=spec['id'], description=spec.get('description', ''))
    for k, v in spec.get('annotations', {}).items():
        rec.annotations[k] = v
    return rec


def alignment_from(spec):
    records = [record_from(s) for s in spec['sequences']]
    a = Alignment(records, np.array(spec['coordinates'], np.intp))
    for k, v in spec.get('attrs', {}).items():
        if k == 'annotations' and 'empty' in v:
            v = dict(v)
            v['empty'] = [(record_from(e['record']), tuple(e['segment']), e['status']) for e in v['empty']]
        setattr(a, k, v)
    return a


def record_spec(rec):
    seq = rec.seq
    n = len(seq)
    try:
        parts = {'0': str(seq)}
    except UndefinedSequenceError:
        data = getattr(seq, '_data', None)
        inner = getattr(data, '_data', None)
        if isinstance(inner, dict):
            parts = {str(k): (v.decode() if isinstance(v, bytes) else str(v)) for k, v in sorted(inner.items())}
        else:
            parts = None
    d = {'id': rec.id, 'length': n, 'parts': parts}
    if rec.description:
        d['description'] = rec.description
    if rec.annotations:
        d['annotations'] = dict(rec.annotations)
    return d


def alignment_spec(a):
    attrs = {}
    for k, v in a.__dict__.items():
        if k in ('sequences', 'coordinates'):
            continue
        if k == 'annotations' and 'empty' in v:
            v = dict(v)
            v['empty'] = [{'record': record_spec(r), 'segment': [int(x) for x in seg], 'status': st} for r, seg, st in v['empty']]
        if isinstance(v, (np.integer,)):
            v = int(v)
        if isinstance(v, (np.floating,)):
            v = float(v)
        attrs[k] = v
    return {'sequences': [record_spec(r) for r in a.sequences], 'coordinates': [[int(x) for x in row] for row in a.coordinates], 'attrs': attrs}


def random_dna(rng, n, lower=0.2):
    s = ''.join(rng.choice('ACGT') for _ in range(n))
    if rng.random() < 0.5:
        s = ''.join(c.lower() if rng.random() < lower else c for c in s)
    if rng.random() < 0.3:
        s = ''.join('N' if rng.random() < 0.05 else c for c in s)
    return s


def random_blocks(rng, max_blocks=6, max_len=25, max_gap=8):
    """A list of (kind, length): 'M' aligned block, 'T' gap in the query (target bases), 'Q' gap in the target."""
    ops = []
    for i in range(rng.randrange(1, max_blocks + 1)):
        if i:
            ops.append((rng.choice('TQ'), rng.randrange(1, max_gap + 1)))
        ops.append(('M', rng.randrange(1, max_len + 1)))
    return ops


def coordinates_from_blocks(ops, t0, q0):
    t, q = t0, q0
    ts, qs = [t], [q]
    for kind, n in ops:
        if kind == 'M':
            t += n
            q += n
        elif kind == 'T':
            t += n
        else:
            q += n
        ts.append(t)
        qs.append(q)
    return [ts, qs]


def random_pair(rng, with_contents=True, strand='+', protein=False):
    """A random nucleotide pairwise alignment (target, query) with the query on the given strand."""
    ops = random_blocks(rng)
    tspan = sum(n for k, n in ops if k in 'MT')
    qspan = sum(n for k, n in ops if k in 'MQ')
    tpre, tpost = rng.randrange(0, 15), rng.randrange(0, 15)
    qpre, qpost = rng.randrange(0, 15), rng.randrange(0, 15)
    tlen, qlen = tpre + tspan + tpost, qpre + qspan + qpost
    co = coordinates_from_blocks(ops, tpre, qpre)
    if strand == '-':
        co[1] = [qlen - x for x in co[1]]
    tname = rng.choice(['chr1', 'chrX', 'target', 'NC_000001.11', 'ctg7'])
    qname = rng.choice(['read1', 'query', 'tx_5', 'q|1'])
    if with_contents:
        tseq = random_dna(rng, tlen)
        qseq = random_dna(rng, qlen)
        tparts, qparts = {'0': tseq}, {'0': qseq}
    else:
        tparts = qparts = None
    return {'sequences': [{'id': tname, 'length': tlen, 'parts': tparts}, {'id': qname, 'length': qlen, 'parts': qparts}], 'coordinates': co, 'attrs': {}}


def random_msa(rng, n=None, with_strands=True):
    """A random multiple alignment of nucleotide rows (printed rows -> sequences + coordinates, flanks, reverse-strand rows)."""
    n = n or rng.randrange(2, 5)
    cols = rng.randrange(6, 50)
    base = ''.join(rng.choice('ACGT') for _ in range(cols))
    rows = []
    for _ in range(n):
        r = ''.join(c if rng.random() > 0.15 else '-' for c in base)
        rows.append(r)
    for c in range(cols):
        if all(r[c] == '-' for r in rows) and rng.random() < 0.9:
            rows[0] = rows[0][:c] + base[c] + rows[0][c + 1:]
    seqs, co = Alignment.parse_printed_alignment([r.encode() for r in rows])
    co = [[int(x) for x in row] for row in co]
    seqs = [s.decode() for s in seqs]
    sequences = []
    for i in range(n):
        pre = ''.join(rng.choice('ACGT') for _ in range(rng.randrange(0, 12)))
        post = ''.join(rng.choice('ACGT') for _ in range(rng.randrange(0, 12)))
        full = pre + seqs[i] + post
        co[i] = [x + len(pre) for x in co[i]]
        if rng.random() < 0.4:
            full = ''.join(c.lower() if rng.random() < 0.2 else c for c in full)
        sequences.append([rng.choice(['hg38.chr1', 'mm10.chr2', 'rn6.chrX', 'sp%d' % i, 'a|b']), full])
    # distinct names
    names = set()
    for i, sq in enumerate(sequences):
        while sq[0] in names:
            sq[0] += '_'
        names.add(sq[0])
    if with_strands:
        for i in range(n):
            if i > 0 and rng.random() < 0.35:
                full = sequences[i][1]
                sequences[i][1] = revcomp(full)
                co[i] = [len(full) - x for x in co[i]]
    return {'sequences': [{'id': a, 'length': len(b), 'parts': {'0': b}} for a, b in sequences], 'coordinates': co, 'attrs': {}}
