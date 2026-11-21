#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/biopython_coord_oracle.tsv (B09 audit A2-27c): the Bio.Align coordinate formats written and read by the installed Biopython (1.88).

Rows: id, format, kind (`write` / `read`), options (JSON: the keyword arguments of the writer), data1, data2:
  write: data1 = JSON list of alignment descriptions (see biopython_coord_oracle_lib.py), data2 = the text `Bio.Align.write` produces;
  read : data1 = the text, data2 = JSON list of the alignments `Bio.Align.parse` returns (same description; `error` = the text of the exception when Biopython refuses the text).
Escaping of the text columns: backslash -> \\\\, newline -> \\n, tab -> \\t, carriage return -> \\r.
Run: python3 scripts/alignment/generate_biopython_coord_oracle.py [format ...]
"""
import io
import json
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from Bio import Align
from biopython_coord_oracle_lib import alignment_from, alignment_spec, random_blocks, random_dna, random_pair

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'biopython_coord_oracle.tsv')
rng = random.Random(20261303)


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def dnax_pair(strand):
    ops = random_blocks(rng, 4, 12, 3)
    ts, qs = [0], [0]
    t = q = 0
    for kind, n in ops:
        if kind == 'M':
            t += 3 * n
            q += n
        elif kind == 'T':
            t += 3 * n
        else:
            q += n
        ts.append(t)
        qs.append(q)
    tpre, tpost = 3 * rng.randrange(0, 4), 3 * rng.randrange(0, 4)
    qpre, qpost = rng.randrange(0, 4), rng.randrange(0, 4)
    tlen, qlen = tpre + t + tpost, qpre + q + qpost
    ts = [x + tpre for x in ts]
    qs = [x + qpre for x in qs]
    if strand == '-':
        ts = [tlen - x for x in ts]
    return {'sequences': [{'id': 'genomic', 'length': tlen, 'parts': None}, {'id': 'prot', 'length': qlen, 'parts': None}], 'coordinates': [ts, qs], 'attrs': {}}


def pslx_pair(strand):
    p = random_pair(rng, True, strand)
    return p


def run_write(fmt, specs, **kw):
    aligns = [alignment_from(s) for s in specs]
    stream = io.StringIO()
    Align.write(aligns, stream, fmt, **kw)
    return stream.getvalue()


def run_read(fmt, text):
    try:
        aligns = list(Align.parse(io.StringIO(text), fmt))
    except Exception as e:  # noqa: BLE001
        return {'error': type(e).__name__}
    return [alignment_spec(a) for a in aligns]


rows = []
wanted = set(sys.argv[1:]) or {'bed', 'psl', 'chain'}


def add_write(fmt, specs, **kw):
    try:
        text = run_write(fmt, specs, **kw)
    except Exception as e:  # noqa: BLE001
        rows.append((fmt, 'write', json.dumps(kw), json.dumps(specs), 'ERROR:' + type(e).__name__))
        return
    rows.append((fmt, 'write', json.dumps(kw), json.dumps(specs), text))


def add_read(fmt, text):
    rows.append((fmt, 'read', '{}', text, json.dumps(run_read(fmt, text))))


def attrs(spec, **kw):
    spec['attrs'].update(kw)
    return spec


if 'bed' in wanted:
    for i in range(110):
        spec = random_pair(rng, rng.random() < 0.7, rng.choice('++-'))
        if rng.random() < 0.6:
            attrs(spec, score=rng.choice([0, 1, 42, 100.5, 7.25, 1e6, -3]))
        if rng.random() < 0.3:
            attrs(spec, thickStart=rng.randrange(0, 20), thickEnd=rng.randrange(20, 60), itemRgb=rng.choice(['0', '255,0,0', '1,2,3']))
        if i % 11 == 0:
            spec['sequences'][1]['id'] = None
        n = rng.choice([12, 12, 3, 4, 5, 6, 7, 8, 9, 10, 11])
        specs = [spec] if i % 5 else [spec, random_pair(rng, False, '-')]
        add_write('bed', specs, bedN=n)
    # texts of other writers / hand-made
    for i in range(70):
        spec = random_pair(rng, False, rng.choice('++-'))
        attrs(spec, score=rng.choice([0, 5, 11.5]), thickStart=3, thickEnd=40, itemRgb='0')
        n = rng.choice([3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
        text = run_write('bed', [spec], bedN=n)
        add_read('bed', text)
        add_read('bed', text + text)
    add_read('bed', 'chr1\t10\t20\n')
    add_read('bed', 'chr1\t10\t20\tname\t0\t-\t10\t20\t0\t2\t4,6,\t0,4,\n')
    add_read('bed', 'chr1\t10\t20\tname\t0\t+\t10\t20\t0\t2\t4,6\t0,4\n')
    add_read('bed', 'chr1\t10\t20\tname\t5.5\t+\n')
    add_read('bed', 'chr1\t10\n')
    add_read('bed', 'chr1\t10\t20\tn\t0\t+\t10\t20\t0\t2\t4\t0,4,\n')
    add_read('bed', 'chr1\t10\t20\tn\t0\t+\t10\t20\t0\t1\t4,\t0,\n')
    add_read('bed', 'chr1 10 20 n x +\n')

if 'psl' in wanted:
    for i in range(150):
        mode = i % 6
        if mode in (0, 1):
            spec = random_pair(rng, rng.random() < 0.8, '+' if mode == 0 else '-')
        elif mode == 2:
            spec = random_pair(rng, rng.random() < 0.5, rng.choice('+-'))
        elif mode == 3:
            spec = dnax_pair('+')
        elif mode == 4:
            spec = dnax_pair('-')
        else:
            spec = random_pair(rng, True, rng.choice('+-'))
        if rng.random() < 0.15:
            attrs(spec, matches=rng.randrange(0, 50), misMatches=rng.randrange(0, 5), repMatches=rng.randrange(0, 5), nCount=rng.randrange(0, 3))
        mask = rng.choice([None, None, 'lower', 'upper'])
        header = rng.random() < 0.8
        kw = {'header': header}
        if mask:
            kw['mask'] = mask
        if rng.random() < 0.15:
            kw['wildcard'] = rng.choice(['N', 'n', None, 'X'])
        specs = [spec] if i % 4 else [spec, random_pair(rng, True, '+')]
        add_write('psl', specs, **kw)
    for i in range(70):
        mode = i % 5
        spec = random_pair(rng, i % 2 == 0, '+' if mode < 2 else '-') if mode < 4 else dnax_pair(rng.choice('+-'))
        text = run_write('psl', [spec], header=(i % 3 != 0))
        add_read('psl', text)
        add_read('psl', text + text.split('\n')[-2] + '\n' if text.strip() else text)
    # pslx files: sequences in the last two columns
    for i in range(30):
        spec = random_pair(rng, True, '+' if i % 2 == 0 else '-')
        text = run_write('psl', [spec], header=True)
        lines = [l for l in text.split('\n') if l]
        head, body = lines[:5], lines[5:]
        out = head[:]
        for l in body:
            w = l.split('\t')
            qs = [int(x) for x in w[19].rstrip(',').split(',')]
            ts = [int(x) for x in w[20].rstrip(',').split(',')]
            bs = [int(x) for x in w[18].rstrip(',').split(',')]
            tseq = spec['sequences'][0]['parts']['0']
            qseq = spec['sequences'][1]['parts']['0']
            if w[8] == '-':
                qseq = qseq[::-1].translate(str.maketrans('ACGTacgtNn', 'TGCAtgcaNn'))
            q_blocks = ','.join(qseq[a:a + n] for a, n in zip(qs, bs)) + ','
            t_blocks = ','.join(tseq[a:a + n] for a, n in zip(ts, bs)) + ','
            out.append('\t'.join(w + [q_blocks, t_blocks]))
        add_read('psl', '\n'.join(out) + '\n')
    add_read('psl', 'psLayout version 3\n\nmatch\tmis- \trep. \tN\'s\tQ gap\tQ gap\tT gap\tT gap\tstrand\tQ        \tQ   \tQ    \tQ  \tT        \tT   \tT    \tT  \tblock\tblockSizes \tqStarts\t tStarts\n     \tmatch\tmatch\t   \tcount\tbases\tcount\tbases\t      \tname     \tsize\tstart\tend\tname     \tsize\tstart\tend\tcount\n---------------------------------------------------------------------------------------------------------------------------------------------------------------\n')
    add_read('psl', '1\t2\t3\n')
    add_read('psl', '')

if 'chain' in wanted:
    for i in range(130):
        spec = random_pair(rng, False, rng.choice('+-'))
        if rng.random() < 0.5:
            t = spec['sequences'][0]
            t['length'] = t['length']
        attrs(spec, score=rng.choice([0, 1, 100, 1234.5, 99999, 3.5, 1e7]))
        if rng.random() < 0.6:
            spec['attrs']['annotations'] = {'id': str(rng.randrange(1, 9999))}
        if rng.random() < 0.3:
            # target on the reverse strand: decreasing target coordinates
            tl = spec['sequences'][0]['length']
            spec['coordinates'][0] = [tl - x for x in spec['coordinates'][0]]
        specs = [spec] if i % 3 else [spec, random_pair(rng, False, '+')]
        for s in specs:
            s['attrs'].setdefault('score', 11)
        add_write('chain', specs)
    for i in range(80):
        spec = random_pair(rng, False, rng.choice('+-'))
        attrs(spec, score=rng.choice([0, 100, 56.5]))
        if i % 2:
            spec['attrs']['annotations'] = {'id': str(i)}
        text = run_write('chain', [spec])
        add_read('chain', text)
        add_read('chain', '# comment\n# another\n' + text + text.replace('chain', 'chain', 1))
        add_read('chain', text.rstrip('\n') + '\n' + text.rstrip('\n') + '\n')   # no blank line between the blocks
    add_read('chain', 'chain 10 t 100 + 0 10 q 100 + 0 10\n10\n\n')
    add_read('chain', 'chain 10 t 100 + 0 10 q 100 + 0 10 7\n10\n')
    add_read('chain', 'chain 10 t 100 x 0 10 q 100 + 0 10\n10\n')
    add_read('chain', 'chain 10 t 100 + 0 10 q 100\n10\n')
    add_read('chain', 'foo 10 t 100 + 0 10 q 100 + 0 10\n10\n')
    add_read('chain', '')

os.makedirs(os.path.dirname(out_path), exist_ok=True)
existing = []
if os.path.exists(out_path) and len(sys.argv) > 1:
    for line in open(out_path):
        if line.startswith('id\t'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[1] not in wanted:
            existing.append(f[1:])
allrows = existing + [list(map(str, r)) if False else [r[0], r[1], esc(r[2]), esc(r[3]), esc(r[4])] for r in rows]
with open(out_path, 'w', newline='\n') as fh:
    fh.write('id\tformat\tkind\toptions\tdata1\tdata2\n')
    for i, r in enumerate(allrows):
        if len(r) == 5 and not isinstance(r[1], str):
            continue
        fh.write('C%04d\t%s\n' % (i + 1, '\t'.join(r)))
print(len(allrows), 'rows;', sum(1 for r in allrows if r[1] == 'write'), 'writer rows;', sum(1 for r in allrows if r[1] == 'read'), 'reader rows;',
      sum(1 for r in allrows if r[4].startswith('ERROR') or r[4].startswith('[{\\"error') or '"error"' in r[4]), 'refusals')
