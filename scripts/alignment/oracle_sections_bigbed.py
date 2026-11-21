"""Oracle sections of the Bio.Align bigbed format (B09 A2-27d): rows of biopython_coord_oracle.tsv, called by generate_biopython_coord_oracle.py.

write rows: options = the writer arguments (bedN, compress, itemsPerSlot, blockSize, extraIndex, declaration), data1 = JSON {"alignments": [spec...], "targets": [record spec...]},
data2 = JSON {"file": base64 of the bytes `Bio.Align.write` produces, "zlib": [[base64 input, base64 output], ...] the calls of zlib.compress} (the compressed bytes depend on the zlib build: the test
feeds them back to the port, so the whole file is compared).  read rows: data1 = base64 of the file, data2 = JSON {alignments, targets, declaration} of `Bio.Align.parse`.
"""
import base64
import io
import json
import zlib

from Bio import Align
from Bio.Align import Alignments, bigbed
from Bio.Align.bigbed import AutoSQLTable, Field
from biopython_coord_oracle_lib import alignment_from, alignment_spec, record_from, record_spec


def b64(b):
    return base64.b64encode(b).decode()


class Spy:
    def __init__(self):
        self.calls = []
        self.orig = zlib.compress

    def __enter__(self):
        def compress(data, *a, **k):
            out = self.orig(data, *a, **k)
            self.calls.append([b64(bytes(data)), b64(out)])
            return out
        zlib.compress = compress
        return self

    def __exit__(self, *e):
        zlib.compress = self.orig


def make_declaration(spec, bedN):
    fields = [Field(*f) for f in spec['fields']]
    return AutoSQLTable(spec['name'], spec['comment'], fields)


def write_bigbed(specs, targets, kw):
    aligns = Alignments([alignment_from(s) for s in specs])
    aligns.targets = [record_from(t) for t in targets]
    args = {k: v for k, v in kw.items() if k != 'declaration'}
    if 'declaration' in kw:
        args['declaration'] = make_declaration(kw['declaration'], kw.get('bedN', 12))
    if 'extraIndex' in args:
        args['extraIndex'] = tuple(args['extraIndex'])
    import os
    import tempfile
    fd, path = tempfile.mkstemp(suffix='.bb')
    os.close(fd)
    try:
        with Spy() as spy:
            Align.write(aligns, path, 'bigbed', **args)
        with open(path, 'rb') as fh:
            data = fh.read()
    finally:
        os.unlink(path)
    return data, spy.calls


def read_bigbed(data):
    try:
        it = Align.parse(io.BytesIO(data), 'bigbed')
        aligns = [a for a in it]
        decl = it.declaration
        return {'alignments': [alignment_spec(a) for a in aligns], 'metadata': {}, 'targets': [record_spec(t) for t in it.targets],
                'declaration': {'name': decl.name, 'comment': decl.comment, 'fields': [list(f) for f in decl]}}
    except Exception as e:  # noqa: BLE001
        return {'error': type(e).__name__}


def random_alignments(rng, n_chrom, bedN, custom):
    names = sorted(rng.sample(['chr1', 'chr2', 'chr3', 'chrX', 'chrY', 'ctg_7', 'scaffold_12', 'z'], n_chrom))
    sizes = {nm: rng.randrange(2000, 6000) for nm in names}
    specs = []
    for nm in names:
        starts = sorted(rng.randrange(0, sizes[nm] - 400) for _ in range(rng.choice([1, 2, 3, 5, 9, 14, 30])))
        for st in starts:
            ops = []
            for i in range(rng.randrange(1, 5)):
                if i:
                    ops.append((rng.choice('TQ'), rng.randrange(1, 40)))
                ops.append(('M', rng.randrange(5, 60)))
            ts, qs = [st], [0]
            t = st
            q = 0
            for kind, k in ops:
                if kind == 'M':
                    t += k
                    q += k
                elif kind == 'T':
                    t += k
                else:
                    q += k
                ts.append(t)
                qs.append(q)
            qlen = q + rng.randrange(0, 20)
            if rng.random() < 0.4:
                qs = [qlen - x for x in qs]
            attrs = {}
            if rng.random() < 0.6:
                attrs['score'] = rng.choice([0, 100, 500, 1000, 250.5, 7.25])
            if rng.random() < 0.4:
                attrs['thickStart'] = st + rng.randrange(0, 5)
                attrs['thickEnd'] = ts[-1] - rng.randrange(0, 5)
            if rng.random() < 0.4:
                attrs['itemRgb'] = rng.choice(['0', '255,0,0', '12,34,56'])
            if custom:
                attrs['annotations'] = custom(rng)
            specs.append({'sequences': [{'id': nm, 'length': sizes[nm], 'parts': None}, {'id': rng.choice(['q1', 'read%d' % rng.randrange(1000), 'tx.5', 'a_b']), 'length': qlen, 'parts': None}],
                          'coordinates': [ts, qs], 'attrs': attrs})
    targets = [{'id': nm, 'length': sizes[nm], 'parts': None} for nm in names]
    return specs, targets


def run(ctx):
    rng = ctx.rng
    rows = ctx.rows
    cases = 0
    for i in range(260):
        bedN = rng.choice([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 12, 12])
        kw = {'bedN': bedN}
        if rng.random() < 0.3:
            kw['compress'] = False
        kw['itemsPerSlot'] = rng.choice([512, 512, 1, 2, 3, 5])
        kw['blockSize'] = rng.choice([256, 256, 2, 3, 4, 7])
        custom = None
        if rng.random() < 0.35:
            fields = [list((f.as_type, f.name, f.comment)) for f in bigbed.AutoSQLTable.default[:bedN]]
            fields += [['string', 'label', 'A label'], ['int', 'count', 'A count'], ['float', 'ratio', 'A ratio'], ['int[count]', 'vals', 'Values']]
            kw['declaration'] = {'name': 'bedX', 'comment': 'custom table', 'fields': fields}

            def custom(r):
                n = r.randrange(1, 4)
                return {'label': r.choice(['alpha', 'b', 'gene one']), 'count': n, 'ratio': r.choice([0.5, 2.25, 10.0]), 'vals': [r.randrange(0, 50) for _ in range(n)]}
            if rng.random() < 0.5:
                kw['extraIndex'] = rng.choice([['label'], ['name'], ['chrom'], ['label', 'name']])
        elif rng.random() < 0.2:
            kw['extraIndex'] = rng.choice([['name'], ['chrom']])
        if 'extraIndex' in kw and bedN <= 3 and 'name' in kw['extraIndex']:
            kw['extraIndex'] = ['chrom']
        specs, targets = random_alignments(rng, rng.choice([1, 1, 2, 3]), bedN, custom)
        try:
            data, calls = write_bigbed(specs, targets, kw)
        except Exception as e:  # noqa: BLE001
            rows.append(('bigbed', 'write', json.dumps(kw), json.dumps({'alignments': specs, 'targets': targets}), json.dumps({'error': type(e).__name__})))
            continue
        rows.append(('bigbed', 'write', json.dumps(kw), json.dumps({'alignments': specs, 'targets': targets}), json.dumps({'file': b64(data), 'zlib': calls})))
        rows.append(('bigbed', 'read', '{}', b64(data), json.dumps(read_bigbed(data))))
        cases += 1
    # refusals and odd inputs
    spec, targets = random_alignments(rng, 1, 12, None)
    bad = json.loads(json.dumps(spec))
    bad[0]['attrs']['score'] = 2000
    for kw, sp in (({'bedN': 12}, bad), ({'bedN': 2}, spec), ({'bedN': 13}, spec), ({'bedN': 12, 'extraIndex': ['nofield']}, spec)):
        try:
            data, calls = write_bigbed(sp, targets, kw)
            rows.append(('bigbed', 'write', json.dumps(kw), json.dumps({'alignments': sp, 'targets': targets}), json.dumps({'file': b64(data), 'zlib': calls})))
        except Exception as e:  # noqa: BLE001
            rows.append(('bigbed', 'write', json.dumps(kw), json.dumps({'alignments': sp, 'targets': targets}), json.dumps({'error': type(e).__name__})))
    # a file that is not bigBed, truncated files
    ok = [r for r in rows if r[0] == 'bigbed' and r[1] == 'read']
    for r in ok[:6]:
        raw = base64.b64decode(r[3])
        for cut in (len(raw) // 2, 100):
            rows.append(('bigbed', 'read', '{}', b64(raw[:cut]), json.dumps(read_bigbed(raw[:cut]))))
    rows.append(('bigbed', 'read', '{}', b64(b'not a bigbed file at all'), json.dumps(read_bigbed(b'not a bigbed file at all'))))
