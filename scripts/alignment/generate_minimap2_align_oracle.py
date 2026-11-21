#!/usr/bin/env python3
"""A3-9c2 (B09 WP45) oracle: the base-level alignment of minimap2 2.26 (`minimap2 -c` / `-a`, `align.c`), the WHOLE standard output of `minimap2 [options] r.fa q.fa` compared as text.

JSON lines: {"set": "...", "args": [...], "files": {"r.fa": "...", "q.fa": "..."}, "stdout": "..."}   (see minimap2_oracle_lib.py)
Groups (default 300 cases each, seed 20261221; the six preset groups share the default as per // 2 = 150 each, 900 cases for the PAF `-c` feature):
  paf_map_ont paf_map_hifi paf_map_pb paf_asm5 paf_asm10 paf_asm20 - `-c` PAF with `--cs` / `--cs=long` / `--MD` / `--eqx`, related queries (1 - 15 % divergence, long indels, N bases, both strands, partial, repeats = secondary hits, chimeras,
                       inversions = Z-drop and inversion rescue)
  paf_options          - `-c` with a random subset of the alignment options (-A -B -O -E -z -s --end-bonus -r -k -w -N -p -M --score-N --max-clip-ratio --no-end-flt --min-dp-len --end-seed-pen --cap-sw-mem --mask-len --hard-mask-level
                       --secondary=no -P --no-long-join) and of the options that must have no effect (-t -K -I -2 -v --bucket-bits)
  sam                  - `-a` (flags 0 / 16 / 256 / 2048, MAPQ, SA tags, header with @PG CL, -Y, -L, --secondary-seq, --sam-hit-only, -R, -y, -Q with FASTQ qualities, --eqx, --cs, --MD, unmapped reads)
  sr                   - `-x sr -c` / `-a` single short reads (the ungapped alignment of the short-read mode)
  inversion            - Z-drop and the inversion rescue (mm_test_zdrop, ksw_ll_i16, mm_align1_inv): an inverted, inserted or deleted stretch in the middle of a related copy, -z / -r / --no-end-flt
  secondary            - repeat families: secondary hits (flag 256, tp:A:S), supplementary hits (2048), -N -p -M --mask-len --hard-mask-level -P --secondary=no
  splice               - `-x splice -c` / `-a` (-u f|r|b|n, -J 0|1, -C, --junc-bonus, --junc-bed BED12 / BED6, --splice-flank, -G, `splice:hq`)
Usage: python3 scripts/alignment/generate_minimap2_align_oracle.py [per_group=300] > minimap2_align_oracle.jsonl      (needs minimap2 2.26 on PATH)
"""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from minimap2_oracle_lib import Gen, fasta, fastq, revcomp, emit

g = Gen(20261221)
rnd = g.rnd


TARGET_POOL = {}


def make_target(hp=False):
    """A reference (1 - 2 sequences); drawn from a pool of ~90 per kind so that reference texts repeat across cases (the oracle stores each once)."""
    key = (hp, rnd.randrange(90))
    if key not in TARGET_POOL:
        ln = rnd.randint(1200, 3500)
        ref = g.rseq(ln, hp, gc=rnd.choice([None, None, 0.7, 0.3]))
        refs = [ref]
        if rnd.random() < 0.3:
            refs.append(g.rseq(rnd.randint(800, 2500), hp))
        TARGET_POOL[key] = refs
    return list(TARGET_POOL[key])


def make_query(refs, hp=False, rate=0.05, short=False):
    """A query drawn from refs[0]: related copy, repeat (second copy), chimera, inversion, N bases."""
    ref = refs[0]
    ln = len(ref)
    kind = rnd.choice(['rel', 'rel', 'rel', 'rel', 'rel', 'rep', 'chim', 'inv', 'inv2', 'gap', 'rel'])
    a = rnd.randint(0, max(0, ln // 2))
    b = rnd.randint(min(ln, a + 400), min(ln, a + 2000))
    seg = ref[a:b]
    if kind == 'rep':
        copy = g.mutate(seg[:len(seg) // 2], rate * 1.5)
        refs[0] = ref + g.rseq(rnd.randint(200, 3000)) + (revcomp(copy) if rnd.random() < 0.4 else copy)
        q = g.mutate(seg, rate)
    elif kind == 'chim':
        c = rnd.randint(0, max(0, ln - 800))
        part2 = ref[c:c + rnd.randint(400, 1500)]
        q = g.mutate(seg[:rnd.randint(300, len(seg))], rate) + g.mutate(revcomp(part2) if rnd.random() < 0.5 else part2, rate)
    elif kind == 'inv':
        k = rnd.randint(150, max(151, len(seg) - 150))
        w = rnd.randint(100, 900)
        q = g.mutate(seg[:k], rate) + g.mutate(revcomp(seg[k:k + w]), rate) + g.mutate(seg[k + w:], rate)
    elif kind == 'inv2':
        q = g.mutate(seg, rate) + g.mutate(revcomp(ref[b:b + rnd.randint(300, 2500)]), rate)
    elif kind == 'gap':
        k = rnd.randint(100, max(101, len(seg) - 100))
        q = g.mutate(seg[:k], rate) + g.rseq(rnd.choice([30, 100, 400, 1200])) + g.mutate(seg[k:], rate)
    else:
        q = g.mutate(seg, rate)
    if rnd.random() < 0.15:
        q = g.with_n(q, rnd.randint(1, 8))
    if rnd.random() < 0.4:
        q = revcomp(q)
    if rnd.random() < 0.2:
        q = q[:rnd.randint(max(1, len(q) // 3), len(q))]
    return q


def valid_gap(preset):
    """-O / -E values the program accepts (E1 > E2, O1 + E1 < O2 + E2, sum <= 127)."""
    for _ in range(50):
        q, e = rnd.randint(1, 30), rnd.randint(1, 6)
        q2, e2 = rnd.randint(q, 80), rnd.randint(0, e - 1)
        if q + e < q2 + e2 and (q + e) + (q2 + e2) <= 127:
            return ['-O', '%d,%d' % (q, q2), '-E', '%d,%d' % (e, e2)]
    return []


ALIGN_OPTS = [
    lambda: ['-A', str(rnd.choice([1, 2, 3, 5]))],
    lambda: ['-B', str(rnd.choice([1, 2, 4, 6, 8, 12]))],
    lambda: valid_gap(None),
    lambda: ['-O', str(rnd.choice([2, 4, 8])) ],
    lambda: ['-E', str(rnd.choice([1, 2, 3]))],
    lambda: ['-z', rnd.choice(['100', '200', '400,200', '50,50', '1000,500'])],
    lambda: ['-s', str(rnd.choice([10, 40, 100, 400, 1000]))],
    lambda: ['--end-bonus', str(rnd.choice([0, 3, 10, 30]))],
    lambda: ['-r', rnd.choice(['100', '300,3000', '1000,50000', '2000'])],
    lambda: ['-k', str(rnd.choice([11, 13, 17, 21]))],
    lambda: ['-w', str(rnd.choice([3, 5, 10, 20]))],
    lambda: ['-N', str(rnd.choice([1, 3, 10, 50]))],
    lambda: ['-p', str(rnd.choice([0.3, 0.6, 0.9, 1.0]))],
    lambda: ['-M', str(rnd.choice([0.2, 0.7]))],
    lambda: ['--score-N', str(rnd.choice([0, 1, 3, 5]))],
    lambda: ['--max-clip-ratio', str(rnd.choice([0.1, 0.4, 0.8]))],
    lambda: ['--no-end-flt'],
    lambda: ['--min-dp-len', str(rnd.choice([20, 100, 500, 2000]))],
    lambda: ['--end-seed-pen', str(rnd.choice([0, 2, 10]))],
    lambda: ['--cap-sw-mem', rnd.choice(['10k', '1m', '100k'])],
    lambda: ['--mask-len', str(rnd.choice([50, 300]))],
    lambda: ['--hard-mask-level'],
    lambda: ['--secondary=no'],
    lambda: ['-P'],
    lambda: ['--no-long-join'],
    lambda: ['-t', str(rnd.choice([1, 2, 5]))],
    lambda: ['-K', rnd.choice(['1', '10k', '100m'])],
    lambda: ['-I', rnd.choice(['1', '10k', '4g'])],
    lambda: ['-2'],
    lambda: ['-v', '0'],
    lambda: ['--bucket-bits', str(rnd.choice([4, 10]))],
]


def paf_case(group, preset, out, extra_pool=None, n_extra=(0, 0)):
    for _ in range(40):
        refs = make_target(hp=preset in ('map-pb',))
        rate = {'map-ont': rnd.choice([0.03, 0.08, 0.15]), 'map-pb': rnd.choice([0.05, 0.12]), 'map-hifi': rnd.choice([0.002, 0.01, 0.03]),
                'asm5': rnd.choice([0.001, 0.004, 0.02]), 'asm10': rnd.choice([0.01, 0.05]), 'asm20': rnd.choice([0.03, 0.1, 0.17])}.get(preset, 0.05)
        q = make_query(refs, preset == 'map-pb', rate)
        files = {'r.fa': fasta(['ref%d' % i for i in range(len(refs))], refs, width=rnd.choice([None, 60, 80])), 'q.fa': fasta(['q%d' % rnd.randint(0, 9999)], [q])}
        args = (['-x', preset] if preset else []) + ['-c']
        if rnd.random() < 0.5:
            args.append(rnd.choice(['--cs', '--cs=long', '--cs=short']))
        if rnd.random() < 0.4:
            args.append('--MD')
        if rnd.random() < 0.25:
            args.append('--eqx')
        if extra_pool:
            for f in rnd.sample(extra_pool, rnd.randint(*n_extra)):
                args += f()
        if emit(group, args + ['r.fa', 'q.fa'], files, out):
            return True
    return False


INV_POOL = {}


def inversion_case(group, out):
    """Z-drop and inversion rescue: a query with an inversion, a long insertion or a long deletion in the middle of a related copy (the DP drops, the chain is split, the inverted part is rescued by ksw_ll)."""
    for _ in range(40):
        key = rnd.randrange(50)
        if key not in INV_POOL:
            INV_POOL[key] = g.rseq(rnd.randint(3500, 7000))
        ref = INV_POOL[key]
        a = rnd.randint(0, 600)
        l1 = rnd.randint(300, 1500)
        kind = rnd.choice(['inv', 'inv', 'ins', 'del', 'invdel'])
        rate = rnd.choice([0.0, 0.01, 0.03, 0.08])
        left = ref[a:a + l1]
        w = rnd.randint(150, 1000)
        mid = ref[a + l1:a + l1 + w]
        right = ref[a + l1 + w:a + l1 + w + rnd.randint(300, 1500)]
        if kind == 'inv':
            q = g.mutate(left, rate, False) + g.mutate(revcomp(mid), rate, False) + g.mutate(right, rate, False)
        elif kind == 'ins':
            q = g.mutate(left, rate, False) + g.rseq(rnd.randint(200, 1500)) + g.mutate(mid + right, rate, False)
        elif kind == 'del':
            q = g.mutate(left, rate, False) + g.mutate(right, rate, False)
        else:
            q = g.mutate(left, rate, False) + g.mutate(revcomp(mid), rate, False) + g.rseq(rnd.randint(200, 600)) + g.mutate(right, rate, False)
        if rnd.random() < 0.4:
            q = revcomp(q)
        files = {'r.fa': fasta(['inv'], [ref]), 'q.fa': fasta(['q%d' % rnd.randint(0, 999)], [q])}
        pre = rnd.choice([None, 'map-ont', 'map-hifi', 'asm5', 'asm10', 'map-pb'])
        args = (['-x', pre] if pre else []) + [rnd.choice(['-c', '-a', '-c'])]
        if rnd.random() < 0.4:
            args.append(rnd.choice(['--cs', '--MD', '--eqx']))
        for f in rnd.sample([['-z', rnd.choice(['50,25', '100', '200,100', '400,200', '1000,500'])], ['-r', rnd.choice(['100', '500,5000', '2000'])], ['--no-end-flt'], ['-s', '30'], ['--min-dp-len', '50'], ['-g', '300'], ['-n', '2'], ['-m', '20']], rnd.randint(0, 2)):
            args += f
        if emit(group, args + ['r.fa', 'q.fa'], files, out):
            return True
    return False


FAMILY_POOL = {}


def secondary_case(group, out):
    """A reference with a family of 2 - 6 diverged copies of a repeat unit (and an inverted copy): the query maps to several places, giving secondary hits and supplementary alignments."""
    for _ in range(40):
        key = rnd.randrange(60)
        if key not in FAMILY_POOL:
            unit = g.rseq(rnd.randint(300, 1100))
            ref = g.rseq(rnd.randint(200, 800))
            for c in range(rnd.randint(2, 6)):
                copy = g.mutate(unit, rnd.choice([0.0, 0.02, 0.06, 0.12]), big_indels=False)
                ref += (revcomp(copy) if rnd.random() < 0.25 else copy) + g.rseq(rnd.randint(50, 400))
            FAMILY_POOL[key] = (unit, ref)
        unit, ref = FAMILY_POOL[key]
        q = g.mutate(unit, rnd.choice([0.0, 0.02, 0.05, 0.1]), big_indels=False)
        if rnd.random() < 0.3:
            q = g.rseq(rnd.randint(20, 150)) + q + g.rseq(rnd.randint(20, 150))
        if rnd.random() < 0.4:
            q = revcomp(q)
        files = {'r.fa': fasta(['fam'], [ref]), 'q.fa': fasta(['q%d' % rnd.randint(0, 999)], [q])}
        pre = rnd.choice([None, 'map-ont', 'map-hifi', 'asm10', 'asm5'])
        args = (['-x', pre] if pre else []) + [rnd.choice(['-c', '-a', '-a', '-c'])]
        if rnd.random() < 0.5:
            args.append(rnd.choice(['--cs', '--MD', '--eqx', '-Y', '--secondary-seq']))
        for f in rnd.sample([['-N', str(rnd.choice([1, 2, 5, 20, 100]))], ['-p', str(rnd.choice([0.1, 0.5, 0.95]))], ['--secondary=no'], ['-M', str(rnd.choice([0.1, 0.3, 0.7]))], ['--mask-len', '200'],
                             ['--hard-mask-level'], ['-P'], ['-z', '100,50'], ['-s', '60'], ['--end-bonus', '10']], rnd.randint(0, 2)):
            args += f
        if emit(group, args + ['r.fa', 'q.fa'], files, out):
            return True
    return False


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    out = sys.stdout
    for preset in ['map-ont', 'map-hifi', 'map-pb', 'asm5', 'asm10', 'asm20']:
        made = 0
        while made < per // 2:
            made += paf_case('paf_' + preset.replace('-', '_'), preset, out)
    made = 0
    while made < per:
        made += paf_case('paf_options', rnd.choice(['', '', 'map-ont', 'map-hifi', 'asm10']), out, ALIGN_OPTS, (1, 4))
    made = 0
    while made < per:
        refs = make_target()
        names = ['ref%d' % i for i in range(len(refs))]
        nreads = rnd.choice([1, 1, 1, 2, 3])
        reads = []
        for k in range(nreads):
            q = make_query(refs, False, rnd.choice([0.02, 0.06, 0.12]))
            reads.append(('read%d_%d' % (made, k), q, g.quals(len(q))))
        fq = rnd.random() < 0.5
        files = {'r.fa': fasta(names, refs), 'q.fa': fastq(reads) if fq else fasta([r[0] for r in reads], [r[1] for r in reads])}
        if rnd.random() < 0.15:
            files['q.fa'] = files['q.fa'].replace('read%d_0' % made, 'read%d_0 BC:Z:ACGT\tXX:i:3' % made, 1)
        if rnd.random() < 0.15:
            files['q.fa'] += ('@' if fq else '>') + 'junk\n' + g.rseq(rnd.randint(60, 200)) + '\n' + ('+\n' + g.quals(60 if False else 0) + '\n' if False else '')
            if fq:
                s = g.rseq(120)
                files['q.fa'] = files['q.fa'].rsplit('@junk', 1)[0] + '@junk\n%s\n+\n%s\n' % (s, g.quals(120))
        args = ['-a'] + ([rnd.choice(['-x map-ont', '-x map-hifi', '-x asm10']).split()] if False else [])
        pre = rnd.choice([None, 'map-ont', 'map-hifi', 'asm10', 'asm20'])
        if pre:
            args = ['-x', pre, '-a']
        for f in [['-Y'], ['--secondary-seq'], ['--sam-hit-only'], ['-y'], ['--eqx'], ['--cs'], ['--cs=long'], ['--MD'], ['-L'], ['--secondary=no'], ['-R', '@RG\\tID:g1\\tSM:s1'], ['-Q'], ['-N', '8']]:
            if rnd.random() < 0.2:
                args += f
        if emit('sam', args + ['r.fa', 'q.fa'], files, out):
            made += 1
    made = 0
    while made < per:
        key = ('sr', rnd.randrange(40))
        if key not in TARGET_POOL:
            ref = g.rseq(rnd.randint(3000, 7000))
            TARGET_POOL[key] = [ref] + ([ref[:rnd.randint(1500, 2500)] + g.rseq(300) + ref[1500:rnd.randint(2500, 3000)]] if rnd.random() < 0.3 else [])
        refs = list(TARGET_POOL[key])
        ref = refs[0]
        ln = len(ref)
        a = rnd.randint(0, ln - 400)
        q = g.mutate(ref[a:a + rnd.randint(80, 300)], rnd.choice([0.0, 0.01, 0.03, 0.06]), big_indels=False)
        if rnd.random() < 0.15:
            q = g.with_n(q, rnd.randint(1, 4))
        if rnd.random() < 0.5:
            q = revcomp(q)
        files = {'r.fa': fasta(['ref%d' % i for i in range(len(refs))], refs), 'q.fa': (fastq([('rd%d' % made, q, g.quals(len(q)))]) if rnd.random() < 0.5 else fasta(['rd%d' % made], [q]))}
        args = ['-x', 'sr', rnd.choice(['-c', '-a', '-c', '--cs'])]
        if rnd.random() < 0.3:
            args += rnd.choice([['--eqx'], ['--MD'], ['-F', '300'], ['-k', '15', '-w', '5'], ['-N', '5'], ['--secondary=yes'], ['-A', '1', '-B', '3', '-O', '4,12', '-E', '2,1']])
        if emit('sr', args + ['r.fa', 'q.fa'], files, out):
            made += 1
    made = 0
    while made < per:
        made += splice_case('splice', out)
    made = 0
    while made < per:
        made += secondary_case('secondary', out)
    made = 0
    while made < per:
        made += inversion_case('inversion', out)


def make_gene(strand_ok=True):
    """A gene: exons (80 - 320 bp) separated by introns (80 - 2500 bp); canonical GT..AG introns (+ strand), CT..AC (- strand, the reverse complement) or random ones."""
    kind = rnd.choice(['gtag', 'gtag', 'ctac', 'gcag', 'atac', 'rand'])
    parts = []
    exons = []
    introns = []
    for i in range(rnd.randint(2, 5)):
        exons.append(g.rseq(rnd.randint(80, 300)))
        if i:
            ln = rnd.randint(80, 1500)
            body = g.rseq(ln)
            if kind == 'gtag':
                intron = 'GT' + body[2:-2] + 'AG'
            elif kind == 'ctac':
                intron = 'CT' + body[2:-2] + 'AC'
            elif kind == 'gcag':
                intron = 'GC' + body[2:-2] + 'AG'
            elif kind == 'atac':
                intron = 'AT' + body[2:-2] + 'AC'
            else:
                intron = body
            introns.append(intron)
    lead = g.rseq(rnd.randint(50, 600))
    trail = g.rseq(rnd.randint(50, 600))
    gene = lead
    spans = []   # BED12: exon starts in the gene
    for i, e in enumerate(exons):
        if i:
            gene += introns[i - 1]
        spans.append((len(gene), len(e)))
        gene += e
    gene += trail
    return gene, exons, spans, kind


GENE_POOL = {}


def splice_case(group, out):
    for _ in range(40):
        key = rnd.randrange(70)
        if key not in GENE_POOL:
            GENE_POOL[key] = make_gene()
        gene, exons, spans, kind = GENE_POOL[key]
        extra_ref = []
        refs = [gene] + extra_ref
        q = g.mutate(''.join(exons), rnd.choice([0.0, 0.01, 0.03, 0.08]), big_indels=False)
        if rnd.random() < 0.15:
            q = g.with_n(q, rnd.randint(1, 3))
        if kind in ('ctac',) and rnd.random() < 0.7:
            q = revcomp(q)
        elif rnd.random() < 0.3:
            q = revcomp(q)
        if rnd.random() < 0.2:
            q = g.rseq(rnd.randint(20, 80)) + q + g.rseq(rnd.randint(0, 60))      # clipped ends
        files = {'r.fa': fasta(['gene%d' % i for i in range(len(refs))], refs), 'q.fa': fasta(['tx%d' % rnd.randint(0, 9999)], [q])}
        pre = rnd.choice(['splice', 'splice', 'splice:hq', 'cdna'])
        args = ['-x', pre, rnd.choice(['-c', '-c', '-a'])]
        if rnd.random() < 0.5:
            args.append(rnd.choice(['--cs', '--cs=long', '--MD', '--eqx']))
        for f in rnd.sample([['-u', rnd.choice('fbrn')], ['-J', rnd.choice('01')], ['-C', str(rnd.choice([0, 3, 9, 20]))], ['--junc-bonus', str(rnd.choice([0, 3, 9, 20]))],
                             ['--splice-flank=' + rnd.choice(['yes', 'no'])], ['-G', rnd.choice(['500', '2k', '10k', '500k'])], ['-z', '300,100'], ['-A', '2', '-B', '4', '-O', '4,24', '-E', '2,1'],
                             ['--end-seed-pen', '2'], ['--no-end-flt'], ['-N', '3'], ['-k', '12', '-w', '4'], ['--score-N', '2']], rnd.randint(0, 3)):
            args += f
        if rnd.random() < 0.45:
            kind_bed = rnd.choice(['bed12', 'bed12', 'bed6'])
            strand = rnd.choice('+-')
            if kind_bed == 'bed12':
                st = spans[0][0]
                en = spans[-1][0] + spans[-1][1]
                bed = 'gene0\t%d\t%d\tg\t0\t%s\t%d\t%d\t0\t%d\t%s\t%s\n' % (st, en, strand, st, en, len(spans), ','.join(str(x[1]) for x in spans) + ',', ','.join(str(x[0] - st) for x in spans) + ',')
            else:
                bed = ''.join('gene0\t%d\t%d\tj\t0\t%s\n' % (spans[i - 1][0] + spans[i - 1][1], spans[i][0], strand) for i in range(1, len(spans)))
            files['j.bed'] = bed.replace('\\t', '\t').replace('\\n', '\n')
            args += ['--junc-bed', 'j.bed']
        if emit(group, args + ['r.fa', 'q.fa'], files, out):
            return True
    return False


if __name__ == '__main__':
    main()
