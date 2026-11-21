"""Oracle sections of the Bio.Align bigpsl format (B09 A2-27d): rows of biopython_coord_oracle.tsv, called by generate_biopython_coord_oracle.py.

write rows: options = {"compress", "cds", "fa", "mask", "wildcard"}, data1 = JSON {"alignments": [spec...], "targets": [record spec...]}, data2 = JSON {"file": base64 of the bytes `Bio.Align.write` produces, "zlib": [[input, output], ...]};
read rows: data1 = base64 of the file, data2 = JSON {alignments, targets, declaration} of `Bio.Align.parse`.  The damaged read rows are bigBed files written with the bigPsl declaration whose custom field was damaged.
"""
import base64
import copy
import io
import json
import os
import tempfile

from Bio import Align
from Bio.Align import Alignments, bigbed, bigpsl
from biopython_coord_oracle_lib import alignment_from, alignment_spec, location_from, record_from, record_spec, revcomp
from oracle_sections_bigbed import Spy, b64

AA = 'ACDEFGHIKLMNPQRSTVWY'


def rand_dna(rng, n, mask=False, n_rate=0.0):
    s = ''.join(rng.choice('ACGT') for _ in range(n))
    if n_rate:
        s = ''.join('N' if rng.random() < n_rate else c for c in s)
    if mask:
        out = []
        i = 0
        while i < len(s):
            k = rng.randrange(1, 12)
            seg = s[i:i + k]
            out.append(seg.lower() if rng.random() < 0.35 else seg)
            i += k
        s = ''.join(out)
    return s


def mutate(rng, s, rate):
    out = []
    for c in s:
        if rng.random() < rate:
            out.append(rng.choice('ACGTN') if c.isupper() else rng.choice('acgtn'))
        else:
            out.append(c)
    return ''.join(out)


def random_location(rng, qsize, strand_hint=None):
    def simple(strand):
        a = rng.randrange(0, max(qsize - 4, 1))
        b = min(qsize, a + rng.choice([0, 1, 1, 2, 5, 9, 20]))
        return {'start': a, 'end': b, 'strand': strand}
    kind = rng.random()
    if kind < 0.45:
        return simple(rng.choice([1, -1, None, 1]))
    op = rng.choice(['join', 'join', 'order'])
    n = rng.randrange(2, 4)
    mode = rng.choice(['plus', 'minus', 'mixed', 'none'])
    parts = []
    pos = 0
    for _ in range(n):
        a = pos + rng.randrange(0, 5)
        b = a + rng.randrange(2, 10)
        parts.append([a, b])
        pos = b
    if mode == 'minus':
        parts = parts[::-1]
    strands = {'plus': [1] * n, 'minus': [-1] * n, 'none': [None] * n, 'mixed': [rng.choice([1, -1]) for _ in range(n)]}[mode]
    if mode == 'mixed' and len(set(strands)) == 1:
        strands[0] = -strands[0]
    return {'parts': [{'start': a, 'end': b, 'strand': st} for (a, b), st in zip(parts, strands)], 'operator': op}


def make_alignment(rng, chroms, kind=None):
    name = rng.choice(list(chroms))
    tsize, ttext = chroms[name]
    kind = kind or rng.choice(['nuc', 'nuc', 'nuc', 'prot'])
    strand = rng.choice('+-')
    nblocks = rng.randrange(1, 6)
    ts, qs = [], []
    t = rng.randrange(0, max(tsize // 4, 1))
    q = 0
    t0 = t
    ts.append(t)
    qs.append(q)
    for i in range(nblocks):
        if i:
            gap = rng.choice('TQ')
            g = rng.randrange(1, 9)
            if gap == 'T':
                t += g
            else:
                q += g
            ts.append(t)
            qs.append(q)
        n = rng.randrange(2, 25)
        t += n * (3 if kind == 'prot' else 1)
        q += n
        ts.append(t)
        qs.append(q)
    if t > tsize:
        return None
    qpre, qpost = rng.randrange(0, 6), rng.randrange(0, 6)
    qsize = q + qpre + qpost
    qs = [x + qpre for x in qs]
    if kind == 'prot':
        qtext = ''.join(rng.choice(AA) for _ in range(qsize))
    else:
        qtext = mutate(rng, ttext[t0 - qpre:t0 - qpre + qsize] if t0 >= qpre and len(ttext) >= t0 - qpre + qsize else rand_dna(rng, qsize), 0.1)
        qtext = (qtext + rand_dna(rng, qsize))[:qsize]
    if strand == '-':
        if kind == 'prot':
            ts = [tsize - x for x in ts]
        else:
            qs = [qsize - x for x in qs]
            qtext = revcomp(qtext)
            qtext = qtext  # the stored query sequence is the forward one; the coordinates run backwards
            qtext = revcomp(qtext)
    tparts = {'0': ttext} if ttext is not None and rng.random() < 0.85 else None
    qmode = rng.random()
    qparts = {'0': qtext} if qmode < 0.8 else None
    spec = {'sequences': [{'id': name, 'length': tsize, 'parts': tparts}, {'id': rng.choice(['q1', 'read_7', 'tx.2', 'a b', 'prot%d' % rng.randrange(100)]), 'length': qsize, 'parts': qparts}],
            'coordinates': [ts, qs], 'attrs': {'annotations': {}}}
    query = spec['sequences'][1]
    if rng.random() < 0.7:
        query['annotations'] = {'molecule_type': 'protein' if kind == 'prot' else rng.choice(['DNA', 'DNA', 'RNA', None])}
        if query['annotations']['molecule_type'] is None:
            del query['annotations']
    if rng.random() < 0.5:
        query['features'] = [{'type': rng.choice(['CDS', 'CDS', 'exon', 'gene']), 'location': random_location(rng, qsize), 'qualifiers': {}}]
        if rng.random() < 0.3:
            query['features'].insert(0, {'type': 'gene', 'location': random_location(rng, qsize), 'qualifiers': {}})
    at = spec['attrs']
    if rng.random() < 0.7:
        at['score'] = rng.choice([0, 100, 250, 1000, 777])
    if rng.random() < 0.4:
        at['thickStart'] = min(ts) + rng.randrange(0, 3)
        at['thickEnd'] = max(ts) - rng.randrange(0, 3)
    if rng.random() < 0.3:
        at['itemRgb'] = rng.choice(['0', '255,0,0', '12,34,56'])
    if rng.random() < 0.15:
        at['matches'] = rng.randrange(0, 90)
    if rng.random() < 0.1:
        at['misMatches'] = rng.randrange(0, 9)
    if rng.random() < 0.08:
        at['repMatches'] = rng.randrange(0, 9)
    if rng.random() < 0.08:
        at['nCount'] = rng.randrange(0, 9)
    return spec


def random_file(rng):
    k = rng.choice([1, 1, 2, 3])
    names = sorted(rng.sample(['chr1', 'chr2', 'chrX', 'ctg_9', 'scaf4'], k))
    chroms = {}
    mask = rng.random() < 0.5
    for nm in names:
        size = rng.randrange(300, 900)
        chroms[nm] = (size, rand_dna(rng, size, mask=mask, n_rate=rng.choice([0, 0, 0.03])))
    specs = []
    for _ in range(rng.choice([1, 2, 3, 5, 8, 12])):
        s = make_alignment(rng, chroms)
        if s:
            specs.append(s)
    targets = [{'id': nm, 'length': chroms[nm][0], 'parts': None} for nm in names]
    return specs, targets


def write_bigpsl(specs, targets, kw):
    aligns = Alignments([alignment_from(s) for s in specs])
    aligns.targets = [record_from(t) for t in targets]
    fd, path = tempfile.mkstemp(suffix='.bb')
    os.close(fd)
    captured = []
    orig = bigbed.AlignmentWriter.write

    def spy_write(self, alignments, *a, **k):
        captured.extend(x for x in (alignment_spec(a) for a in alignments) if 'oChromStart' in x['attrs'].get('annotations', {}))
        return orig(self, alignments, *a, **k)
    bigbed.AlignmentWriter.write = spy_write
    try:
        with Spy() as spy:
            Align.write(aligns, path, 'bigpsl', **kw)
        with open(path, 'rb') as fh:
            data = fh.read()
    finally:
        bigbed.AlignmentWriter.write = orig
        os.unlink(path)
    return data, spy.calls, captured, aligns.targets


def read_bigpsl(data):
    try:
        it = Align.parse(io.BytesIO(data), 'bigpsl')
        aligns = [a for a in it]
        decl = it.declaration
        return {'alignments': [alignment_spec(a) for a in aligns], 'targets': [record_spec(t) for t in it.targets],
                'declaration': {'name': decl.name, 'comment': decl.comment, 'fields': [list(f) for f in decl]}}
    except Exception as e:  # noqa: BLE001
        return {'error': type(e).__name__}


def write_raw(specs, targets):
    fd, path = tempfile.mkstemp(suffix='.bb')
    os.close(fd)
    try:
        aligns = Alignments([alignment_from(x) for x in specs])
        aligns.targets = targets
        Align.write(aligns, path, 'bigbed', bedN=12, declaration=bigpsl.declaration)
        with open(path, 'rb') as fh:
            return fh.read()
    finally:
        os.unlink(path)


DAMAGE = [
    ('oChromStart', lambda v: str(int(v) + 1)), ('oChromEnd', lambda v: str(int(v) + 2)), ('oStrand', lambda v: '-' if v == '+' else '+'),
    ('oChromSize', lambda v: str(int(v) + 5)), ('oChromStarts', lambda v: v + ',3'), ('oChromStarts', lambda v: ''), ('oSequence', lambda v: 'ACGT'),
    ('oCDS', lambda v: 'join(1..5,9..'), ('oCDS', lambda v: '3..1'), ('oCDS', lambda v: 'complement(2..5)'), ('oCDS', lambda v: '5'), ('oCDS', lambda v: '6^7'), ('oCDS', lambda v: '<1..>20'),
    ('oCDS', lambda v: 'order(1..3,5..9)'), ('oCDS', lambda v: 'nonsense'), ('chromSize', lambda v: str(int(v) + 1)), ('seqType', lambda v: '3'), ('seqType', lambda v: '2'), ('seqType', lambda v: '1'),
    ('match', lambda v: 'x'), ('misMatch', lambda v: '-1'), ('nCount', lambda v: '2.5'),
]


def run(ctx):
    rng = ctx.rng
    rows = ctx.rows
    for _ in range(330):
        specs, targets = random_file(rng)
        kw = {}
        if rng.random() < 0.3:
            kw['compress'] = False
        if rng.random() < 0.4:
            kw['cds'] = True
        if rng.random() < 0.3 and all(sq['parts'] is not None for sp in specs for sq in sp['sequences'][1:]):
            kw['fa'] = True
        kw['mask'] = rng.choice([None, None, 'lower', 'upper'])
        if rng.random() < 0.3:
            kw['wildcard'] = rng.choice([78, 78, 'N', 'n', 110])
        else:
            kw['wildcard'] = 'N'
        try:
            data, calls, captured, tg = write_bigpsl(specs, targets, kw)
        except Exception as e:  # noqa: BLE001
            rows.append(('bigpsl', 'write', json.dumps(kw), json.dumps({'alignments': specs, 'targets': targets}), json.dumps({'error': type(e).__name__})))
            continue
        rows.append(('bigpsl', 'write', json.dumps(kw), json.dumps({'alignments': specs, 'targets': targets}), json.dumps({'file': b64(data), 'zlib': calls})))
        rows.append(('bigpsl', 'read', '{}', b64(data), json.dumps(read_bigpsl(data))))
    # a file with 3 sequences, without annotations, with an impossible block ratio
    specs, targets = random_file(rng)
    for variant in range(6):
        sp = copy.deepcopy(specs)
        if variant == 0:
            del sp[0]['attrs']['annotations']
        elif variant == 1:
            sp[0]['coordinates'][1][1] = sp[0]['coordinates'][1][1] + 1
        elif variant == 2:
            sp[0]['sequences'].append(copy.deepcopy(sp[0]['sequences'][1]))
            sp[0]['coordinates'].append(list(sp[0]['coordinates'][1]))
        elif variant == 3:
            sp[0]['sequences'][1]['parts'] = None
            kw = {'fa': True}
        elif variant == 4:
            sp[0]['coordinates'] = [[sp[0]['coordinates'][0][0], sp[0]['coordinates'][0][0]], [0, 0]]
        else:
            sp[0]['sequences'][0]['id'] = 'unknown_chrom'
        kw = {'fa': True} if variant == 3 else {}
        try:
            data, calls, captured, tg = write_bigpsl(sp, targets, kw)
            res = {'file': b64(data), 'zlib': calls}
        except Exception as e:  # noqa: BLE001
            res = {'error': type(e).__name__}
        rows.append(('bigpsl', 'write', json.dumps(kw), json.dumps({'alignments': sp, 'targets': targets}), json.dumps(res)))
    # damaged custom fields
    for _ in range(14):
        specs, targets = random_file(rng)
        try:
            data, calls, captured, tg = write_bigpsl(specs, targets, {'cds': True, 'fa': rng.random() < 0.5})
        except Exception:  # noqa: BLE001
            continue
        if not captured:
            continue
        for field, fn in DAMAGE:
            fixed = copy.deepcopy(captured)
            a = rng.choice(fixed)
            a['attrs']['annotations'][field] = fn(a['attrs']['annotations'][field])
            try:
                raw = write_raw(fixed, tg)
            except Exception:  # noqa: BLE001
                continue
            rows.append(('bigpsl', 'read', '{}', b64(raw), json.dumps(read_bigpsl(raw))))
    # oCDS location strings of the whole INSDC grammar (A5-8c3e, F239 follow-up): within / one-of positions, references, compound locations, fuzzy bounds, refusals
    odd = ['(8.9)..(20.25)', '8..(20.25)', '(8.9)..25', 'one-of(6,9)..12', '3..one-of(8,12,5)', 'one-of(6,9)..one-of(20,22)', '(8^9)..25', '<(8.9)..25', '12..>20', 'NC_1.2:3..10', 'complement(NC_1.2:3..10)',
           'join(NC_1:3..10,20..30)', 'join(1..5,NC_1.2:20..30)', '(8.9)', 'one-of(6,9)', 'bond(5)', 'bond(<5)', '<5', '10..<5', '5..3', '(9.8)..12', '(8.9)..(3.4)', '0..5', '(0.5)..9', 'one-of(6,9)..2',
           'join(one-of(3,5)..9,12..15)', 'join(1..5,6..8,10..14)', 'complement(join(10..14,6..8,1..5))', 'order(1..3,complement(5..9))', 'join(complement(5..9),1..3)', 'join(1..5)', 'join(1..5,', 'join()', 'join(1..5,,7..9)',
           'join(a..b,2..3)', 'bond(1..3,5..7)', 'complement(complement(2..5))', 'join(complement(1..3),complement(5..9))', 'a:b:1..5', '1..5xyz', '2^3', '3^4', '3^5', '1.5', '..5', '5..', '1,5', ' 1..5', '1..5 ',
           'join(1..5,7..9)extra', 'order(1..5)', 'bond(1..5)', 'join(<1..>5,7..9)', 'complement(<1..>5)', 'join(1..3, 5..9)', 'JOIN(1..3,5..9)', '-3..5', '1..-5', '3..3', '4..4', 'join(3..3,6..6)']
    for k in range(rng.randrange(1, 3)):
        specs, targets = random_file(rng)
        try:
            data, calls, captured, tg = write_bigpsl(specs, targets, {'cds': True})
        except Exception:  # noqa: BLE001
            continue
        if not captured:
            continue
        for text in odd:
            fixed = copy.deepcopy(captured)
            a = rng.choice(fixed)
            a['attrs']['annotations']['oCDS'] = text
            try:
                raw = write_raw(fixed, tg)
            except Exception:  # noqa: BLE001
                continue
            rows.append(('bigpsl', 'read', '{}', b64(raw), json.dumps(read_bigpsl(raw))))
    ok = [r for r in rows if r[0] == 'bigpsl' and r[1] == 'read'][:6]
    for r in ok:
        raw = base64.b64decode(r[3])
        for cut in (len(raw) // 2, 100):
            rows.append(('bigpsl', 'read', '{}', b64(raw[:cut]), json.dumps(read_bigpsl(raw[:cut]))))
    rows.append(('bigpsl', 'read', '{}', b64(b'not a bigpsl file'), json.dumps(read_bigpsl(b'not a bigpsl file'))))
    # a bigBed file with another declaration is refused by the bigPsl reader
    spec = {'sequences': [{'id': 'chr1', 'length': 100, 'parts': None}, {'id': 'q', 'length': 10, 'parts': None}], 'coordinates': [[0, 10], [0, 10]], 'attrs': {}}
    aligns = Alignments([alignment_from(spec)])
    aligns.targets = [record_from({'id': 'chr1', 'length': 100, 'parts': None})]
    fd, path = tempfile.mkstemp(suffix='.bb')
    os.close(fd)
    try:
        Align.write(aligns, path, 'bigbed')
        raw = open(path, 'rb').read()
    finally:
        os.unlink(path)
    rows.append(('bigpsl', 'read', '{}', b64(raw), json.dumps(read_bigpsl(raw))))
