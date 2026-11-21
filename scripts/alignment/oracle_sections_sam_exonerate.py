"""Oracle sections of the Bio.Align formats sam and exonerate (B09 A2-27d): rows of biopython_coord_oracle.tsv, called by generate_biopython_coord_oracle.py."""
import copy

import numpy as np
from Bio import Align
from biopython_coord_oracle_lib import random_blocks, random_dna, random_pair

BASES = 'ACGT'


def rec(i, length, rng, desc=None, ann=None):
    d = {'id': i, 'length': length, 'parts': None}
    if desc is not None:
        d['description'] = desc
    if ann:
        d['annotations'] = ann
    return d


def sam_ops(rng, spec, with_ops):
    """Optionally attach `operations` to a random pair (M, I, D and N for a target gap)."""
    co = spec['coordinates']
    ops = ''
    for k in range(len(co[0]) - 1):
        dt = co[0][k + 1] - co[0][k]
        dq = abs(co[1][k + 1] - co[1][k])
        if dt and dq:
            ops += 'M'
        elif dt:
            ops += rng.choice('DDN')
        else:
            ops += 'I'
    return ops


def random_tags(rng):
    ann = {}
    for _ in range(rng.randrange(0, 4)):
        kind = rng.choice('ifAZHiiAZZ')
        tag = rng.choice(['NM', 'XS', 'ZZ', 'RG', 'XY', 'UQ'])
        if kind == 'i':
            ann[tag] = rng.randrange(-5, 300)
        elif kind == 'f':
            ann[tag] = rng.choice([0.5, 1.25, -3.0, 1e-5, 100.0, 12345.678])
        elif kind == 'A':
            ann[tag] = rng.choice('ACGT!')
        elif kind == 'Z':
            ann[tag] = rng.choice(['grp1', 'hello world', 'x:y', ''])
        elif kind == 'H':
            ann[tag] = {'__bytes__': [rng.randrange(0, 100) for _ in range(rng.randrange(1, 4))]}
        elif kind == 'B':
            ann[tag] = {'__array__': [rng.randrange(0, 9) for _ in range(rng.randrange(1, 4))], 'dtype': 'int'}
    return ann


def sam_pair(rng, strand, defined=True):
    spec = random_pair(rng, defined, strand)
    q = spec['sequences'][1]
    t = spec['sequences'][0]
    t['id'] = rng.choice(['chr1', 'chr2', 'chrX', 'ref'])
    q['id'] = rng.choice(['read1', 'r/2', 'seq_7'])
    if defined and rng.random() < 0.6:
        # keep the query unmasked, uppercase ACGT only (the writer reverse-complements with the DNA table)
        q['parts'] = {'0': ''.join(rng.choice(BASES) for _ in range(q['length']))}
        t['parts'] = {'0': ''.join(rng.choice(BASES) for _ in range(t['length']))}
    return spec


def run(ctx):
    rng = ctx.rng
    add_write, add_read, run_write = ctx.add_write, ctx.add_read, ctx.run_write
    meta_choices = [
        None,
        {'HD': {'VN': '1.6', 'SO': 'coordinate'}},
        {'HD': {'VN': '1.5'}, 'PG': [{'ID': 'bwa', 'PN': 'bwa', 'VN': '0.7'}], 'CO': [{'XX': 'y'}]},
        {'RG': [{'ID': 'g1', 'SM': 'a'}, {'ID': 'g2', 'SM': 'b'}]},
    ]
    target_choices = [
        None,
        [rec('chr1', 1000, rng), rec('chr2', 500, rng)],
        [rec('chr1', 1000, rng, desc='<unknown description>', ann={'assembly': 'hg38', 'species': 'human', 'MD5': 'abc', 'names': ['a', 'b'], 'topology': 'linear', 'URI': 'file:x', 'alternate_locus': 'chr1_alt'})],
        [rec('ctg', 77, rng, desc='some desc', ann={'XX': 'val', 'topology': 'circular'})],
        [rec('chrM', 16569, rng, desc='')],
    ]
    # ---- writers
    for i in range(130):
        strand = rng.choice('+-')
        spec = sam_pair(rng, strand, defined=(i % 4 != 0))
        a = spec['attrs']
        if rng.random() < 0.5:
            a['score'] = rng.choice([0, 5, 77, 12.5, 100.5, 3.5, 2.5, 1e6])
        if rng.random() < 0.4:
            a['flag'] = rng.choice([0, 1, 2, 256, 1024, 64, 128, 2048])
        if rng.random() < 0.5:
            a['mapq'] = rng.choice([0, 30, 60, 255, 1])
        if rng.random() < 0.3:
            a['rnext'] = rng.choice(['chr1', 'chr2', 'other', spec['sequences'][0]['id']])
            a['pnext'] = rng.randrange(0, 500)
        if rng.random() < 0.4:
            a['annotations'] = random_tags(rng)
        q = spec['sequences'][1]
        if rng.random() < 0.3:
            q['annotations'] = {'hard_clip_left': rng.randrange(1, 20)}
            if rng.random() < 0.5:
                q['annotations']['hard_clip_right'] = rng.randrange(1, 20)
        elif rng.random() < 0.15:
            q['annotations'] = {'hard_clip_right': rng.randrange(1, 9)}
        if rng.random() < 0.3:
            q['letter_annotations'] = {'phred_quality': [rng.randrange(0, 60) for _ in range(q['length'])]}
        if i % 3 == 0:
            a['operations'] = sam_ops(rng, spec, True)
        if i % 17 == 0:
            a['operations'] = sam_ops(rng, spec, True).replace('M', rng.choice('=X'), 1)
        kw = {}
        md = (i % 4 != 0) and rng.random() < 0.7
        if md:
            kw['md'] = True
        m = rng.choice(meta_choices)
        if m is not None:
            kw['metadata'] = m
        tg = rng.choice(target_choices)
        if tg is not None:
            kw['targets'] = tg
        specs = [spec] if i % 5 else [spec, sam_pair(rng, '+', True)]
        add_write('sam', specs, **kw)
    # unmapped and odd writers
    add_write('sam', [{'sequences': [None, {'id': 'r1', 'length': 4, 'parts': {'0': 'ACGT'}}], 'coordinates': None, 'attrs': {'flag': 4}}])
    spec = sam_pair(rng, '+', True)
    spec['sequences'][1]['id'] = None
    add_write('sam', [spec])
    spec = sam_pair(rng, '+', True)
    spec['coordinates'] = [[spec['coordinates'][0][0]], [spec['coordinates'][1][0]]]
    add_write('sam', [spec])
    spec = sam_pair(rng, '+', True)
    spec['attrs']['annotations'] = {'AB': {'__array__': [1, 2], 'dtype': 'int'}}
    add_write('sam', [spec])
    spec = sam_pair(rng, '+', True)
    spec['attrs']['annotations'] = {'AB': {'k': 1}}
    add_write('sam', [spec])
    add_write('sam', [], metadata={'HD': {'VN': '1.6'}}, targets=[rec('c', 5, rng)])
    add_write('sam', [], metadata={'HD': {'SO': 'x'}})
    add_write('sam', [sam_pair(rng, '+', False)], targets=[rec('c', 5, rng, ann={'topology': 'bogus'})])

    # ---- readers: texts of the writer
    for i in range(70):
        spec = sam_pair(rng, rng.choice('+-'), defined=True)
        a = spec['attrs']
        if rng.random() < 0.5:
            a['score'] = rng.choice([0, 5, 77, 100])
        if rng.random() < 0.5:
            a['mapq'] = rng.choice([0, 30, 60])
        q = spec['sequences'][1]
        if rng.random() < 0.3:
            q['annotations'] = {'hard_clip_left': rng.randrange(1, 20)}
        if rng.random() < 0.3:
            q['letter_annotations'] = {'phred_quality': [rng.randrange(0, 60) for _ in range(q['length'])]}
        if i % 3 == 0:
            a['operations'] = sam_ops(rng, spec, True)
        t = spec['sequences'][0]
        m = rng.choice(meta_choices)
        tg = [rec(n, 100000, rng) for n in ('chr1', 'chr2', 'chrX', 'ref')] + [rec('other', 40, rng)]
        kw = {'targets': tg}
        if m is not None:
            kw['metadata'] = m
        if rng.random() < 0.6:
            kw['md'] = True
        try:
            text = run_write('sam', [spec, sam_pair(rng, '+', True)] if i % 6 == 0 else [spec], **kw)
        except Exception:  # noqa: BLE001
            continue
        add_read('sam', text)
    # hand-made lines
    def cigar_line(rng, with_md, flag=None, extra=''):
        n = rng.randrange(1, 6)
        cig = ''
        qlen = 0
        tlen = 0
        ref_seg = ''
        md = ''
        run_match = 0
        seq = ''
        ops = []
        if rng.random() < 0.4:
            k = rng.randrange(1, 6)
            cig += '%dS' % k
            seq += ''.join(rng.choice(BASES) for _ in range(k))
        for j in range(n):
            if j:
                kind = rng.choice('IDN')
            else:
                kind = 'M'
            if kind == 'M':
                k = rng.randrange(1, 12)
                letter = rng.choice(['M', 'M', 'M', '=', 'X']) if rng.random() < 0.4 else 'M'
                cig += '%d%s' % (k, letter)
                for _ in range(k):
                    r = rng.choice(BASES)
                    if rng.random() < 0.15:
                        qb = rng.choice([b for b in BASES if b != r])
                        md += str(run_match) + r
                        run_match = 0
                    else:
                        qb = r
                        run_match += 1
                    seq += qb
                    ref_seg += r
            elif kind == 'I':
                k = rng.randrange(1, 4)
                cig += '%dI' % k
                seq += ''.join(rng.choice(BASES) for _ in range(k))
            elif kind == 'D':
                k = rng.randrange(1, 4)
                cig += '%dD' % k
                dseq = ''.join(rng.choice(BASES) for _ in range(k))
                md += str(run_match) + '^' + dseq
                run_match = 0
                ref_seg += dseq
            else:
                k = rng.randrange(5, 30)
                cig += '%dN' % k
        md += str(run_match)
        if rng.random() < 0.4:
            k = rng.randrange(1, 6)
            cig += '%dS' % k
            seq += ''.join(rng.choice(BASES) for _ in range(k))
        if rng.random() < 0.3:
            cig = '%dH' % rng.randrange(1, 9) + cig
        if rng.random() < 0.3:
            cig += '%dH' % rng.randrange(1, 9)
        if flag is None:
            flag = rng.choice([0, 0, 16, 16, 256, 4, 99])
        fields = ['read%d' % rng.randrange(100), str(flag), rng.choice(['chr1', 'chr2']), str(rng.randrange(1, 200)), str(rng.choice([0, 30, 255, 60])),
                  cig, rng.choice(['*', '=', 'chr2']), str(rng.choice([0, 17])), str(rng.choice([0, 100, -50])), seq if rng.random() > 0.1 else '*',
                  ''.join(chr(33 + rng.randrange(0, 40)) for _ in range(len(seq))) if rng.random() < 0.4 else '*']
        if with_md:
            fields.append('MD:Z:' + md)
        if rng.random() < 0.4:
            fields.append('AS:i:%d' % rng.randrange(0, 100))
        for _ in range(rng.randrange(0, 3)):
            fields.append(rng.choice(['NM:i:3', 'XS:A:+', 'XT:Z:hello', 'XF:f:1.5', 'XH:H:1234', 'XB:B:i1,2,3', 'XB:B:f1.5,2', 'XC:Z:a:b']))
        return '\t'.join(fields) + extra + '\n'

    header = '@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:5000\tM5:xyz\tTP:linear\n@SQ\tSN:chr2\tLN:3000\tDS:second\tAH:*\n@PG\tID:p\tPN:q\n@CO\tXX:free text\n'
    for i in range(60):
        lines = ''.join(cigar_line(rng, rng.random() < 0.5) for _ in range(rng.randrange(1, 4)))
        add_read('sam', (header if i % 4 else '') + lines)
    for i in range(30):
        add_read('sam', header + cigar_line(rng, True, flag=rng.choice([0, 16])))
    add_read('sam', '')
    add_read('sam', header)
    add_read('sam', '@SQ\tSN:chr1\n')
    add_read('sam', '@SQ\tSN:chr1\tLN:5\tTP:x\n')
    add_read('sam', '@SQ\tSN:chr1\tLN:5\tZZ:q\n@RG\tID:1\n')
    add_read('sam', '@HD\tVN:1\nread\t0\tchr1\t1\n')
    add_read('sam', 'r\t4\t*\t0\t255\t*\t*\t0\t0\tACGT\t*\n')
    add_read('sam', 'r\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n')
    add_read('sam', 'r\t0\tchr9\t1\t255\t4M\t*\t0\t0\tACGT\t*\n')
    add_read('sam', '@SQ\tSN:chr1\tLN:50\nr\t0\tchr9\t1\t255\t4M\t*\t0\t0\tACGT\t*\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\t*\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t4P\t*\t0\t0\tACGT\t*\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t3M\t*\t0\t0\tACGT\t*\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t*\t*\t0\t0\t*\t*\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\t*\tNM:i\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\t*\tXB:B:x1,2\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\t*\tXB:B:i,1\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\t*\tAS:f:1\n')
    add_read('sam', '\n')
    add_read('sam', 'r\t0\tchr1\t1\t255\t4M\t*\t0\t0\tACGT\t*\tXY:i:abc\n')


# ---------------------------------------------------------------------- exonerate

def exo_pair(rng, strands, protein=None):
    """Alignment with the molecule types protein on `protein` in {None, 't', 'q'}; steps in residues (x3 for DNA)."""
    n = rng.randrange(1, 7)
    ts, qs = [0], [0]
    ops = ''
    for k in range(n):
        kind = rng.choice('MMMDI') if k else 'M'
        step = rng.randrange(1, 15)
        if protein is None:
            if kind == 'M':
                dt = dq = step
            elif kind == 'D':
                dt, dq = step, 0
            else:
                dt, dq = 0, step
        elif protein == 'q':
            if kind == 'M':
                dt, dq = 3 * step, step
            elif kind == 'D':
                dt, dq = 3 * step, 0
            else:
                dt, dq = 0, step
        else:
            if kind == 'M':
                dt, dq = step, 3 * step
            elif kind == 'D':
                dt, dq = step, 0
            else:
                dt, dq = 0, 3 * step
        ts.append(ts[-1] + dt)
        qs.append(qs[-1] + dq)
        ops += {'M': 'M', 'D': 'D', 'I': 'I'}[kind]
    tpre, qpre = rng.randrange(0, 5), rng.randrange(0, 5)
    tt, qq = ts[-1], qs[-1]
    tlen, qlen = tpre + tt + rng.randrange(0, 5), qpre + qq + rng.randrange(0, 5)
    ts = [x + tpre for x in ts]
    qs = [x + qpre for x in qs]
    if strands[0] == '-':
        ts = [tlen - x for x in ts]
    if strands[1] == '-':
        qs = [qlen - x for x in qs]
    t = {'id': rng.choice(['gene1', 'chr3', 'T_1']), 'length': tlen, 'parts': None}
    q = {'id': rng.choice(['est7', 'prot', 'q2']), 'length': qlen, 'parts': None}
    if protein == 't' or strands[0] == '.':
        t['annotations'] = {'molecule_type': 'protein'}
    if protein == 'q' or strands[1] == '.':
        q['annotations'] = {'molecule_type': 'protein'}
    return {'sequences': [t, q], 'coordinates': [ts, qs], 'attrs': {'score': rng.choice([0, 12, 100.5, 250, 3.25, 1e6])}}, ops


def run_exonerate(ctx):
    rng = ctx.rng
    add_write, add_read, run_write = ctx.add_write, ctx.add_read, ctx.run_write
    for i in range(160):
        mode = i % 5
        if mode == 0:
            spec, ops = exo_pair(rng, (rng.choice('+-'), rng.choice('+-')))
        elif mode == 1:
            spec, ops = exo_pair(rng, ('+', '+'), protein='q')
        elif mode == 2:
            spec, ops = exo_pair(rng, ('+', '+'), protein='t')
        elif mode == 3:
            spec, ops = exo_pair(rng, ('+', '+'), protein='q')
            spec['sequences'][0]['annotations'] = {'molecule_type': 'DNA'}
        else:
            spec, ops = exo_pair(rng, (rng.choice('+-'), rng.choice('+-')))
        kw = {'exfmt': rng.choice(['vulgar', 'cigar'])}
        if rng.random() < 0.5:
            kw['metadata'] = {'Command line': rng.choice(['exonerate -m est2genome a.fa b.fa', 'x', '']), 'Hostname': rng.choice(['host1', ''])}
        if i % 3 == 0:
            # operations (vulgar-style): M, D, I plus N introns, 5 / 3 splice sites, U, S, F, C
            nops = len(spec['coordinates'][0]) - 1
            o = ''
            for k in range(nops):
                dt = spec['coordinates'][0][k + 1] - spec['coordinates'][0][k]
                dq = spec['coordinates'][1][k + 1] - spec['coordinates'][1][k]
                dt, dq = abs(dt), abs(dq)
                if dt and dq:
                    o += rng.choice('MMMC53FS') if dt == dq else 'M'
                elif dt:
                    o += rng.choice('DDNU')
                else:
                    o += rng.choice('IIU')
            spec['attrs']['operations'] = o
        specs = [spec] if i % 4 else [spec, exo_pair(rng, ('+', '+'))[0]]
        add_write('exonerate', specs, **kw)
    for i in range(40):
        steps = []
        ops = ''
        for k in range(rng.randrange(1, 4)):
            a = rng.randrange(1, 9)
            steps.append((a, a))
            ops += 'M'
            kind = rng.choice(['U', 'N', '5', '3', 'D', 'I', 'F', 'S', 'C'])
            if kind == 'U':
                steps.append((rng.randrange(1, 9), 0))
                steps.append((0, rng.randrange(1, 9)))
                ops += 'UU'
            elif kind in 'N53DF':
                steps.append((2 if kind in '53' else rng.randrange(1, 40), 0))
                ops += kind
            elif kind in 'IS':
                steps.append((0, rng.randrange(1, 5)))
                ops += kind
            else:
                steps.append((3, 3))
                ops += 'C'
        ts, qs = [2], [1]
        for dt, dq in steps:
            ts.append(ts[-1] + dt)
            qs.append(qs[-1] + dq)
        spec = {'sequences': [{'id': 'g', 'length': ts[-1] + 3, 'parts': None}, {'id': 'q', 'length': qs[-1] + 2, 'parts': None}], 'coordinates': [ts, qs], 'attrs': {'score': rng.randrange(1, 500), 'operations': ops}}
        for fmt in ('vulgar', 'cigar'):
            add_write('exonerate', [spec], exfmt=fmt)
            try:
                text = run_write('exonerate', [spec], exfmt=fmt)
            except Exception:  # noqa: BLE001
                continue
            add_read('exonerate', text)
    add_write('exonerate', [], metadata={'Command line': 'c', 'Hostname': 'h'})
    spec, _ = exo_pair(rng, ('+', '+'))
    del spec['attrs']['score']
    add_write('exonerate', [spec], exfmt='vulgar')
    add_write('exonerate', [spec], exfmt='bogus')

    # readers: the writers' text and hand-made lines
    for i in range(80):
        mode = i % 4
        spec, ops = exo_pair(rng, (rng.choice('+-'), rng.choice('+-')), protein=None if mode < 2 else rng.choice('qt'))
        o = None
        if i % 2:
            nops = len(spec['coordinates'][0]) - 1
            o = ''
            for k in range(nops):
                dt = abs(spec['coordinates'][0][k + 1] - spec['coordinates'][0][k])
                dq = abs(spec['coordinates'][1][k + 1] - spec['coordinates'][1][k])
                o += rng.choice('MMC') if dt and dq and dt == dq else ('M' if dt and dq else ('D' if dt else 'I'))
            spec['attrs']['operations'] = o
        for fmt in ('vulgar', 'cigar'):
            try:
                text = run_write('exonerate', [spec] if i % 5 else [spec, spec], exfmt=fmt, metadata={'Command line': 'exonerate x', 'Hostname': 'h'})
            except Exception:  # noqa: BLE001
                continue
            add_read('exonerate', text)
    hdr = 'Command line: [exonerate -m est2genome a b]\nHostname: [box]\n'
    foot = '-- completed exonerate analysis\n'
    add_read('exonerate', hdr + 'vulgar: q 0 10 + t 100 110 + 50 M 10 10\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 30 + t 100 200 + 50 M 10 10 5 0 2 I 0 50 3 0 2 M 18 18\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 30 + t 100 200 + 50 M 10 10 N 5 6 G 0 3 G 4 0\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 10 . t 100 130 + 50 M 10 30\n' + foot)
    add_read('exonerate', hdr + 'cigar: q 0 10 + t 100 110 + 50 M 10\n' + foot)
    add_read('exonerate', hdr + 'cigar: q 0 10 + t 110 100 - 50 M 5 I 2 D 3 M 3\n' + foot)
    add_read('exonerate', hdr + 'cigar: q 0 10 . t 100 130 + 50 M 30\n' + foot)
    add_read('exonerate', hdr + 'cigar: q 0 10 + t 100 110 . 50 D 30 M 5\n' + foot)
    add_read('exonerate', hdr + 'cigar: q 0 10 + t 100 110 + 50 X 10\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 10 + t 100 110 + 50 Z 10 10\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 10 + t 100 110 + 50 G 3 4\n' + foot)
    add_read('exonerate', hdr + 'something else\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 10 + t 100 110 + 50 M 10 10\n')
    add_read('exonerate', hdr + foot + '\n')
    add_read('exonerate', hdr + foot + 'x\n')
    add_read('exonerate', hdr + foot)
    add_read('exonerate', 'bad\n')
    add_read('exonerate', '')
    add_read('exonerate', 'Command line: [x]\nHostname: h\n' + foot)
    add_read('exonerate', 'Command line: x]\nHostname: [h]\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 10 + t 100 110 x 50 M 10 10\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 10 + t 100 110 + abc M 10 10\n' + foot)
    add_read('exonerate', hdr + 'vulgar: q 0 10 + t 100 110 + 50\n' + foot)
