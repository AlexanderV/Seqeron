#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/alignment_result_cigar_oracle.jsonl (B09 R14-9, F302): the extended CIGAR, the MD tag and the cs string that
parasail 1.3.4 and minimap2 2.26 write for alignments, described by the rows / the CIGAR the tool wrote.

kind = "parasail": ref, query (the sequences), fn (nw / sg / sg_qx / sg_dx / sg_qb ...), row1 (the aligned reference row), row2 (the aligned query row), cigar (extended CIGAR of the library), score.
kind = "minimap2": ref, query, tstart, qstart, cg (PAF `cg` with --eqx), md, cs_short, cs_long (PAF tags of `minimap2 -c --eqx --MD --cs=short|long`).
Run: python3 scripts/alignment/generate_alignment_result_cigar_oracle.py
"""
import json
import os
import random
import subprocess
import tempfile

import parasail

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'alignment_result_cigar_oracle.jsonl')
rng = random.Random(20261403)


def dna(n):
    return ''.join(rng.choice('ACGT') for _ in range(n))


def mutate(seq, rate=0.08):
    out = []
    i = 0
    while i < len(seq):
        r = rng.random()
        if r < rate / 2:
            out.append(rng.choice('ACGT'))
        elif r < rate * 0.7:
            i += rng.randrange(1, 4)
            continue
        elif r < rate:
            out.append(seq[i] + dna(rng.randrange(1, 4)))
        else:
            out.append(seq[i])
        i += 1
    return ''.join(out)


def parasail_rows():
    rows = []
    dna_matrix = parasail.matrix_create('ACGT', 2, -3)
    fns = ['nw', 'sg', 'sg_qx', 'sg_dx', 'sg_qb', 'sg_db', 'sg_qe', 'sg_de', 'sg_qb_de', 'sg_qe_db']
    for i in range(200):
        if i % 4 == 3:
            alphabet, matrix = 'ACDEFGHIKLMNPQRSTVWY', parasail.blosum62
            ref = ''.join(rng.choice(alphabet) for _ in range(rng.randrange(8, 40)))
            query = ''.join(c for c in mutate_protein(ref, alphabet))
            if not query:
                query = 'A'
            gap_open, gap_ext = 11, 1
        else:
            ref = dna(rng.randrange(8, 50))
            query = mutate(ref) or 'A'
            matrix = dna_matrix
            gap_open, gap_ext = 5, 2
        for fn in fns[i % 3::3]:
            f = getattr(parasail, fn + '_trace')
            res = f(query, ref, gap_open, gap_ext, matrix)
            rows.append({'kind': 'parasail', 'id': 'p%03d_%s' % (i, fn), 'fn': fn, 'ref': ref, 'query': query,
                         'row1': res.traceback.ref, 'row2': res.traceback.query,
                         'cigar': res.cigar.decode.decode(), 'score': res.score})
    return rows


def mutate_protein(seq, alphabet):
    out = []
    for ch in seq:
        r = rng.random()
        if r < 0.1:
            out.append(rng.choice(alphabet))
        elif r < 0.15:
            continue
        elif r < 0.2:
            out.append(ch + rng.choice(alphabet))
        else:
            out.append(ch)
    return out


def minimap2_rows():
    rows = []
    tmp = tempfile.mkdtemp()
    for i in range(120):
        ref = dna(rng.randrange(200, 500))
        a = rng.randrange(0, 60)
        b = len(ref) - rng.randrange(0, 60)
        query = mutate(ref[a:b], rng.choice([0.02, 0.05, 0.08]))
        tf, qf = os.path.join(tmp, 't.fa'), os.path.join(tmp, 'q.fa')
        open(tf, 'w').write('>ref\n' + ref + '\n')
        open(qf, 'w').write('>read\n' + query + '\n')
        def run(*opts):
            out = subprocess.run(['minimap2', '-c', *opts, '-k11', '-w5', '--secondary=no', tf, qf], capture_output=True, text=True, check=True).stdout
            lines = [l for l in out.split('\n') if l]
            return lines[0].split('\t') if lines else None

        # minimap2 writes MD instead of cs when both are given: three runs
        eqx, short, long_ = run('--eqx', '--MD'), run('--cs=short'), run('--cs=long')
        if eqx is None or eqx[4] != '+' or short is None or long_ is None:
            continue
        f = eqx
        tags = {t[:2]: t[5:] for t in eqx[12:]}
        short_tags = {t[:2]: t[5:] for t in short[12:]}
        long_tags = {t[:2]: t[5:] for t in long_[12:]}
        rows.append({'kind': 'minimap2', 'id': 'm%03d' % i, 'ref': ref, 'query': query, 'tstart': int(f[7]), 'tend': int(f[8]), 'qstart': int(f[2]), 'qend': int(f[3]),
                     'cg': tags['cg'], 'md': tags['MD'], 'cs_short': short_tags['cs'], 'cs_long': long_tags['cs']})
    return rows


rows = parasail_rows() + minimap2_rows()
with open(out_path, 'w', newline='\n') as fh:
    for r in rows:
        fh.write(json.dumps(r) + '\n')
print(len(rows), 'rows ->', out_path, '(parasail', sum(r['kind'] == 'parasail' for r in rows), ', minimap2', sum(r['kind'] == 'minimap2' for r in rows), ')')
