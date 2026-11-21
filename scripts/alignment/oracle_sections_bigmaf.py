"""Oracle sections of the Bio.Align bigmaf format (B09 A2-27d): rows of biopython_coord_oracle.tsv, called by generate_biopython_coord_oracle.py.

write rows: options = {"compress": bool}, data1 = JSON {"alignments": [spec...], "targets": [record spec...]}, data2 = JSON {"file": base64 of the bytes `Bio.Align.write` produces, "zlib": [[base64 input, base64 output], ...]}
(the compressed bytes depend on the zlib build: the test feeds them back to the port).  read rows: data1 = base64 of the file, data2 = JSON {alignments, targets, declaration} of `Bio.Align.parse`.
The corrupted read rows are bigBed files whose `mafBlock` field was written with the bigbed writer from a damaged text of the block.
"""
import base64
import io
import json
import os
import tempfile

from Bio import Align
from Bio.Align import Alignments, bigmaf
from biopython_coord_oracle_lib import alignment_from, alignment_spec, random_msa, record_from, record_spec
from oracle_sections_bigbed import Spy, b64


def write_bigmaf(specs, targets, compress):
    aligns = Alignments([alignment_from(s) for s in specs])
    aligns.targets = [record_from(t) for t in targets]
    fd, path = tempfile.mkstemp(suffix='.bb')
    os.close(fd)
    try:
        with Spy() as spy:
            Align.write(aligns, path, 'bigmaf', compress=compress)
        with open(path, 'rb') as fh:
            data = fh.read()
    finally:
        os.unlink(path)
    return data, spy.calls


def read_bigmaf(data):
    try:
        it = Align.parse(io.BytesIO(data), 'bigmaf')
        aligns = [a for a in it]
        decl = it.declaration
        return {'alignments': [alignment_spec(a) for a in aligns], 'targets': [record_spec(t) for t in it.targets],
                'declaration': {'name': decl.name, 'comment': decl.comment, 'fields': [list(f) for f in decl]}}
    except Exception as e:  # noqa: BLE001
        return {'error': type(e).__name__}


def maf_spec(rng, chromosome):
    spec = random_msa(rng, rng.choice([2, 2, 3, 4, 5]), with_strands=True)
    seqs = spec['sequences']
    seqs[0]['id'] = 'hg38.' + chromosome
    used = {seqs[0]['id']}
    for sq in seqs[1:]:
        while sq['id'] in used:
            sq['id'] += '_'
        used.add(sq['id'])
    if rng.random() < 0.6:
        spec['attrs']['score'] = rng.choice([0, 12.5, 1234.567891, -3, 99])
    ann = {}
    if rng.random() < 0.3:
        ann['pass'] = rng.randrange(1, 9)
    if rng.random() < 0.25:
        ann['empty'] = [{'record': {'id': 'emp%d' % k, 'length': rng.randrange(50, 200), 'parts': None}, 'segment': [rng.randrange(0, 20), rng.randrange(20, 40)], 'status': rng.choice('CIMn')} for k in range(rng.randrange(1, 3))]
    if ann:
        spec['attrs']['annotations'] = ann
    for sq, co in zip(seqs, spec['coordinates']):
        a = {}
        if rng.random() < 0.25:
            a['quality'] = ''.join(rng.choice('0123456789F') for _ in range(abs(co[-1] - co[0])))
        if rng.random() < 0.25:
            a.update({'leftStatus': rng.choice('CINnMT'), 'leftCount': rng.randrange(0, 500), 'rightStatus': rng.choice('CINnMT'), 'rightCount': rng.randrange(0, 500)})
        if a:
            sq['annotations'] = a
    return spec


def random_file(rng):
    k = rng.choice([1, 1, 2, 3])
    chroms = rng.sample(['chr1', 'chr2', 'chr3', 'chrX', 'chrM', 'scaf_9'], k)
    specs = []
    for ch in chroms:
        for _ in range(rng.choice([1, 2, 3, 5, 8, 15])):
            specs.append(maf_spec(rng, ch))
    rng.shuffle(specs)
    lengths = {}
    for s in specs:
        ch = s['sequences'][0]['id'][5:]
        lengths[ch] = max(lengths.get(ch, 0), s['sequences'][0]['length'])
    targets = [{'id': 'hg38.' + chroms[0], 'length': lengths[chroms[0]], 'parts': None}] + [{'id': ch, 'length': lengths[ch], 'parts': None} for ch in chroms[1:]]
    return specs, targets


def damaged_blocks(rng, text):
    """Damaged copies of a MAF block text (`;` for the line ends)."""
    lines = text.split(';')
    out = []
    s_idx = [i for i, l in enumerate(lines) if l.startswith('s ')]
    if s_idx:
        i = rng.choice(s_idx)
        w = lines[i].split()
        for mut in (lambda: w[:3] + [str(int(w[3]) + 1)] + w[4:],            # wrong size
                    lambda: w[:2] + ['-3'] + w[3:],                         # negative start
                    lambda: w[:6] + [w[6].replace('-', '.')],               # dots as gaps
                    lambda: w[:-1],                                         # six fields
                    lambda: w[:4] + ['*'] + w[5:],                          # wrong strand
                    lambda: w[:5] + [str(int(w[5]) - 3)],                   # source too short
                    lambda: w[:6] + [w[6] + 'AC']):                         # longer row
            ws = mut()
            ls = list(lines)
            ls[i] = ' '.join(ws)
            out.append(';'.join(ls))
    out.append(';'.join(lines[1:]))                                        # no a line
    out.append('#comment;' + text)
    out.append(text.replace('a score=', 'a scor=', 1))
    out.append(text.replace('a', 'a pass=0', 1) if text.startswith('a') else text)
    out.append(text + 'x y z;')
    out.append(text.replace(';', ';q ', 1))
    out.append(text.replace(';', ';i x C 0 C 0;', 1))
    out.append(text.replace(';', ';e x 1 2 + 10 C;', 1))
    out.append(text.replace(';', ';e x 1 2 + 10 Q;', 1))
    out.append(text.replace(';', ';i hg38.chr1 C 0 C 0;', 1))
    out.append(text.replace('a', 'a score=1.5e2', 1) if text.startswith('a') else text)
    return out


def run(ctx):
    rng = ctx.rng
    rows = ctx.rows
    cases = 0
    for _ in range(230):
        specs, targets = random_file(rng)
        compress = rng.random() > 0.3
        kw = {'compress': compress}
        try:
            data, calls = write_bigmaf(specs, targets, compress)
        except Exception as e:  # noqa: BLE001
            rows.append(('bigmaf', 'write', json.dumps(kw), json.dumps({'alignments': specs, 'targets': targets}), json.dumps({'error': type(e).__name__})))
            continue
        rows.append(('bigmaf', 'write', json.dumps(kw), json.dumps({'alignments': specs, 'targets': targets}), json.dumps({'file': b64(data), 'zlib': calls})))
        rows.append(('bigmaf', 'read', '{}', b64(data), json.dumps(read_bigmaf(data))))
        cases += 1
    # refusals: no dot in the reference name, a reverse-strand reference row, a single row, no targets
    specs, targets = random_file(rng)
    bad = json.loads(json.dumps(specs))
    bad[0]['sequences'][0]['id'] = 'nodot'
    rows.append(('bigmaf', 'write', json.dumps({'compress': True}), json.dumps({'alignments': bad, 'targets': targets}), json.dumps({'error': 'ValueError'}) if False else None))
    for sp, tg in ((bad, targets),):
        try:
            data, calls = write_bigmaf(sp, tg, True)
            res = {'file': b64(data), 'zlib': calls}
        except Exception as e:  # noqa: BLE001
            res = {'error': type(e).__name__}
        rows[-1] = ('bigmaf', 'write', json.dumps({'compress': True}), json.dumps({'alignments': sp, 'targets': tg}), json.dumps(res))
    # damaged blocks: write the bigBed file with the text as the mafBlock and read it as bigmaf
    from Bio.Align import bigbed
    good = [r for r in rows if r[0] == 'bigmaf' and r[1] == 'read']
    for r in rng.sample(good, 40):
        raw = base64.b64decode(r[3])
        parsed = Align.parse(io.BytesIO(raw), 'bigmaf')
        first = next(iter(parsed))
        block = format(first, 'maf')[:-1].replace('\n', ';')
        chrom = first.sequences[0].id.split('.', 1)[1]
        for text in damaged_blocks(rng, block):
            length = len(first.sequences[0])
            pair = {'sequences': [{'id': chrom, 'length': length, 'parts': None}, {'id': 'x', 'length': 5, 'parts': None}], 'coordinates': [[0, 4], [0, 4]], 'attrs': {'annotations': {'mafBlock': text}}}
            fd, path = tempfile.mkstemp(suffix='.bb')
            os.close(fd)
            try:
                aligns = Alignments([alignment_from(pair)])
                aligns.targets = [record_from({'id': chrom, 'length': length, 'parts': None})]
                Align.write(aligns, path, 'bigbed', bedN=3, declaration=bigmaf.declaration)
                with open(path, 'rb') as fh:
                    data = fh.read()
            finally:
                os.unlink(path)
            rows.append(('bigmaf', 'read', '{}', b64(data), json.dumps(read_bigmaf(data))))
    ok = [r for r in rows if r[0] == 'bigmaf' and r[1] == 'read'][:6]
    for r in ok:
        raw = base64.b64decode(r[3])
        for cut in (len(raw) // 2, 100):
            rows.append(('bigmaf', 'read', '{}', b64(raw[:cut]), json.dumps(read_bigmaf(raw[:cut]))))
    rows.append(('bigmaf', 'read', '{}', b64(b'not a bigmaf file'), json.dumps(read_bigmaf(b'not a bigmaf file'))))
