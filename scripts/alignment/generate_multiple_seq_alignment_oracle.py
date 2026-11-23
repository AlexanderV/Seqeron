#!/usr/bin/env python3
"""B09 WP67 (F307, R14-10b): oracle of the legacy Bio.Align.MultipleSeqAlignment 1.88 (+ the SeqRecord operations it uses): construction, append / extend, lengths, str(), row / rows / cell / column / sub-alignment indexing,
del, +, sort, substitutions (with sequence weights), the .alignment property and the column annotations (the restricted dictionary).
JSONL.GZ rows: id, records, annotations, column_annotations, ops (a list of {op, ...}), results (one per op: {kind: msa|rec|str|int|matrix|alignment|none|error, ...}, an error is {kind: error, type, message}).
Usage: python3 scripts/alignment/generate_multiple_seq_alignment_oracle.py OUT.jsonl.gz [n=400] [seed]      (needs Biopython 1.88)"""
import gzip, json, random, sys, warnings
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
warnings.simplefilter('ignore')


def mk_record(d):
    return SeqRecord(Seq(d['seq']), id=d['id'], name=d['name'], description=d['description'], annotations=dict(d.get('annotations', {})), letter_annotations=dict(d.get('letter_annotations', {})), dbxrefs=list(d.get('dbxrefs', [])))


def snap_rec(r):
    return {'id': r.id, 'name': r.name, 'description': r.description, 'seq': str(r.seq), 'annotations': dict(r.annotations), 'letter_annotations': {k: (v if isinstance(v, str) else list(v)) for k, v in r.letter_annotations.items()}, 'dbxrefs': list(r.dbxrefs)}


def snap_val(v):
    return v if isinstance(v, str) else list(v)


def snap_msa(m):
    return {'records': [snap_rec(r) for r in m], 'annotations': dict(m.annotations), 'column_annotations': {k: snap_val(v) for k, v in m.column_annotations.items()}}


def sl(x):
    return slice(*x) if isinstance(x, list) else x


def run_op(m, op, rnd_extra):
    k = op['op']
    if k == 'len': return {'kind': 'int', 'value': len(m)}
    if k == 'alignment_length': return {'kind': 'int', 'value': m.get_alignment_length()}
    if k == 'str': return {'kind': 'str', 'value': str(m)}
    if k == 'row': return {'kind': 'rec', 'value': snap_rec(m[op['i']])}
    if k == 'rows': return {'kind': 'msa', 'value': snap_msa(m[sl(op['s'])])}
    if k == 'row_cols': return {'kind': 'rec', 'value': snap_rec(m[op['i'], sl(op['c'])])}
    if k == 'cell': return {'kind': 'str', 'value': m[op['i'], op['j']]}
    if k == 'column': return {'kind': 'str', 'value': m[sl(op['s']), op['j']]}
    if k == 'slice': return {'kind': 'msa', 'value': snap_msa(m[sl(op['s']), sl(op['c'])])}
    if k == 'add':
        other = MultipleSeqAlignment([mk_record(r) for r in op['other']['records']], annotations=dict(op['other']['annotations']), column_annotations=dict(op['other']['column_annotations']))
        return {'kind': 'msa', 'value': snap_msa(m + other)}
    if k == 'sort':
        key = {'id': None, 'len': lambda r: len(r.seq), 'seq': lambda r: str(r.seq), 'description': lambda r: r.description, 'gaps': lambda r: str(r.seq).count('-')}[op['key']]
        m.sort(key=key, reverse=op['reverse'])
        return {'kind': 'msa', 'value': snap_msa(m)}
    if k == 'substitutions':
        a = m.substitutions
        return {'kind': 'matrix', 'letters': ''.join(a.alphabet), 'values': [[float(a[i, j]) for j in range(len(a.alphabet))] for i in range(len(a.alphabet))]}
    if k == 'alignment':
        a = m.alignment
        rows = []
        for r in a.sequences:
            rows.append({'id': r.id, 'description': r.description, 'seq': str(r.seq), 'annotations': dict(r.annotations), 'letter_annotations': {kk: snap_val(v) for kk, v in r.letter_annotations.items()}})
        return {'kind': 'alignment', 'coordinates': a.coordinates.tolist() if len(a.sequences) else [], 'rows': rows, 'annotations': dict(a.annotations), 'column_annotations': {kk: snap_val(v) for kk, v in a.column_annotations.items()}}
    if k == 'append': m.append(mk_record(op['record'])); return {'kind': 'msa', 'value': snap_msa(m)}
    if k == 'extend': m.extend([mk_record(r) for r in op['records']]); return {'kind': 'msa', 'value': snap_msa(m)}
    if k == 'del_row': del m[op['i']]; return {'kind': 'msa', 'value': snap_msa(m)}
    if k == 'del_rows': del m[sl(op['s'])]; return {'kind': 'msa', 'value': snap_msa(m)}
    if k == 'set_column_annotation': m.column_annotations[op['key']] = op['value']; return {'kind': 'msa', 'value': snap_msa(m)}
    if k == 'set_column_annotations': m.column_annotations = dict(op['value']); return {'kind': 'msa', 'value': snap_msa(m)}
    if k == 'column_annotations': return {'kind': 'msa', 'value': snap_msa(m)}
    raise ValueError(k)


def build(rnd):
    ncols = rnd.randint(0, 14)
    n = rnd.choice([0, 1, 2, 3, 3, 4, 5, 6, 21, 22]) if rnd.random() < 0.9 else rnd.randint(1, 4)
    alpha = rnd.choice(['ACGT', 'ACGT-', 'ACGT--', 'ACDEFGHIKLMNPQRSTVWY-', 'ACGTN-acgt'])
    recs = []
    for i in range(n):
        seq = ''.join(rnd.choice(alpha) for _ in range(ncols))
        d = {'id': rnd.choice(['s%d' % i, 'seq_%d' % i, 'Alpha', 'Beta', 'Gamma', 'x']) if rnd.random() < 0.8 else '<unknown id>', 'name': rnd.choice(['n%d' % i, 'same', '<unknown name>']),
             'description': rnd.choice(['d%d' % i, 'same description', '', '<unknown description>']), 'seq': seq}
        if rnd.random() < 0.3: d['annotations'] = {'weight': rnd.choice([0.5, 0.8, 1.0, 2, 1.5]), 'tool': 'demo'}
        elif rnd.random() < 0.2: d['annotations'] = {'tool': rnd.choice(['demo', 'other']), 'molecule_type': 'DNA'}
        if rnd.random() < 0.25: d['letter_annotations'] = {'q': [rnd.randint(0, 40) for _ in range(ncols)], 'ss': ''.join(rnd.choice('HEC') for _ in range(ncols))}
        if rnd.random() < 0.15: d['dbxrefs'] = rnd.choice([['db:1'], ['db:1', 'db:2'], []])
        recs.append(d)
    return recs, ncols


def make_ops(rnd, recs, ncols):
    n = len(recs)

    def rslice():
        return [rnd.choice([None, None, -3, -1, 0, 1, 2, 5]), rnd.choice([None, None, -2, 1, 2, 3, 8]), rnd.choice([None, None, 1, 2, -1, 3])]

    ops = [{'op': 'len'}, {'op': 'alignment_length'}, {'op': 'str'}]
    for _ in range(rnd.randint(2, 7)):
        t = rnd.choice(['row', 'rows', 'row_cols', 'cell', 'column', 'slice', 'add', 'sort', 'substitutions', 'alignment', 'append', 'extend', 'del_row', 'del_rows', 'column_annotations', 'set_column_annotation', 'set_column_annotations', 'str'])
        if t == 'row': ops.append({'op': 'row', 'i': rnd.randint(-n - 1, n + 1)})
        elif t == 'rows': ops.append({'op': 'rows', 's': rslice()})
        elif t == 'row_cols': ops.append({'op': 'row_cols', 'i': rnd.randint(-n - 1, n), 'c': rslice()})
        elif t == 'cell': ops.append({'op': 'cell', 'i': rnd.randint(-n - 1, n), 'j': rnd.randint(-ncols - 1, ncols + 1)})
        elif t == 'column': ops.append({'op': 'column', 's': rslice() if rnd.random() < 0.6 else [None, None, None], 'j': rnd.randint(-ncols - 1, ncols + 1)})
        elif t == 'slice': ops.append({'op': 'slice', 's': rslice() if rnd.random() < 0.5 else [None, None, None], 'c': rslice()})
        elif t == 'add':
            other_recs = []
            extra = rnd.randint(0, 6)
            for i in range(n if rnd.random() < 0.9 else n + 1):
                src = recs[i] if i < n and rnd.random() < 0.7 else {'id': 'z%d' % i, 'name': 'n', 'description': 'd'}
                seq = ''.join(rnd.choice('ACGT-') for _ in range(extra))
                d = {'id': src['id'] if rnd.random() < 0.6 else 'o%d' % i, 'name': src['name'], 'description': src['description'] if rnd.random() < 0.5 else 'other', 'seq': seq}
                if rnd.random() < 0.3 and 'annotations' in src: d['annotations'] = dict(src['annotations'])
                if rnd.random() < 0.3: d['letter_annotations'] = {'q': [rnd.randint(0, 40) for _ in range(extra)]} if rnd.random() < 0.7 else {'q': 'x' * extra}
                if rnd.random() < 0.2: d['dbxrefs'] = ['db:2']
                other_recs.append(d)
            ops.append({'op': 'add', 'other': {'records': other_recs, 'annotations': rnd.choice([{}, {'tool': 'demo'}]), 'column_annotations': {'stats': 'C' * extra} if (other_recs and rnd.random() < 0.5) else {}}})
        elif t == 'sort': ops.append({'op': 'sort', 'key': rnd.choice(['id', 'id', 'len', 'seq', 'description', 'gaps']), 'reverse': rnd.random() < 0.4})
        elif t == 'substitutions': ops.append({'op': 'substitutions'})
        elif t == 'alignment': ops.append({'op': 'alignment'})
        elif t == 'append':
            L = ncols if rnd.random() < 0.8 else ncols + 1
            ops.append({'op': 'append', 'record': {'id': 'new', 'name': 'nn', 'description': 'new row', 'seq': ''.join(rnd.choice('ACGT-') for _ in range(L))}})
        elif t == 'extend':
            L = ncols if rnd.random() < 0.8 else ncols + 2
            ops.append({'op': 'extend', 'records': [{'id': 'e%d' % i, 'name': 'e', 'description': 'e', 'seq': ''.join(rnd.choice('ACGT-') for _ in range(L))} for i in range(rnd.randint(0, 3))]})
        elif t == 'del_row': ops.append({'op': 'del_row', 'i': rnd.randint(-n - 1, n + 1)})
        elif t == 'del_rows': ops.append({'op': 'del_rows', 's': rslice()})
        elif t == 'column_annotations': ops.append({'op': 'column_annotations'})
        elif t == 'set_column_annotation': ops.append({'op': 'set_column_annotation', 'key': rnd.choice(['stats', 'cons']), 'value': rnd.choice([''.join(rnd.choice('CX') for _ in range(ncols)), [rnd.random() for _ in range(ncols)], 'toolong' * 3, 5 if False else 'a'])})
        elif t == 'set_column_annotations': ops.append({'op': 'set_column_annotations', 'value': rnd.choice([{}, {'stats': 'C' * ncols}, {'stats': 'C' * (ncols + 1)}])})
    return ops


def main():
    out = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 400; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20262502
    rnd = random.Random(seed)
    rows = []
    while len(rows) < n:
        recs, ncols = build(rnd)
        if rnd.random() < 0.08 and len(recs) > 1:
            recs[-1]['seq'] += 'A'          # unequal rows: the constructor refuses
        ann = rnd.choice([{}, {}, {'tool': 'demo'}, {'tool': 'demo', 'name': 'start'}])
        cann = {}
        if recs and rnd.random() < 0.4: cann = {'stats': ''.join(rnd.choice('CX') for _ in range(ncols))}
        if recs and rnd.random() < 0.1: cann['bad'] = 'x' * (ncols + 1)
        if not recs and rnd.random() < 0.2: cann = {'stats': 'abc'}
        row = {'id': len(rows), 'records': recs, 'annotations': ann, 'column_annotations': cann}
        try:
            m = MultipleSeqAlignment([mk_record(r) for r in recs], annotations=dict(ann), column_annotations=dict(cann))
        except Exception as e:
            row['ctor_error'] = {'type': type(e).__name__, 'message': str(e)}
            rows.append(row)
            continue
        row['ops'] = make_ops(rnd, recs, ncols)
        row['results'] = []
        for op in row['ops']:
            try:
                row['results'].append(run_op(m, op, rnd))
            except Exception as e:
                row['results'].append({'kind': 'error', 'type': type(e).__name__, 'message': str(e)})
        rows.append(row)
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    errs = sum(1 for r in rows for x in r.get('results', []) if x['kind'] == 'error')
    print(len(rows), 'rows;', sum(len(r.get('ops', [])) for r in rows), 'ops;', errs, 'op errors;', sum(1 for r in rows if 'ctor_error' in r), 'constructor errors')


if __name__ == '__main__':
    main()
