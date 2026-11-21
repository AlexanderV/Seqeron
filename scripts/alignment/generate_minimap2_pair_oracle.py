#!/usr/bin/env python3
"""A5-9b / A5-9c (B09 WP45) oracle: the multi-segment (paired-read) mode of minimap2 2.26 - the WHOLE standard output of `minimap2 -x sr [options] r.fa q1.fa q2.fa` (two query files) and `minimap2 -x sr [options] r.fa q.fa` (interleaved reads
with the same name, `--frag=yes`), PAF without and with `-c`, SAM `-a`.

JSON lines like generate_minimap2_align_oracle.py. Groups (default 300 cases each, seed 20261222):
  pe_files_paf   - two query files, PAF without -c  (collect_minimizers with segment ids, comput_sc multi-segment rule, mm_select_sub_multi, mm_seg_gen, per-segment mm_set_parent / mapq, pe_ori FR)
  pe_files_pafc  - two query files, PAF -c (+ --cs / --MD / --eqx): mm_pair, the paired mapping qualities, proper pairs, pe_thru
  pe_files_sam   - two query files, SAM (-a: flags 0x1 0x2 0x8 0x20 0x40 0x80, mate fields, TLEN, MAPQ)
  pe_inter_paf   - one interleaved file, PAF -c
  pe_inter_sam   - one interleaved file, SAM (also /1 /2 name suffixes, three reads with one name = FI tags)
  pe_ff          - --frag=yes with the default preset (pe_ori 0 = FF, no flipping), PAF / SAM
  pe_options     - -x sr with -F -N -p -k -w -A -B -O -E -s -z --end-bonus --max-clip-ratio --secondary=yes --heap-sort=no --no-pairing -r -g -n -m -M --mask-level and options without effect (-t -K -I -2)
Usage: python3 scripts/alignment/generate_minimap2_pair_oracle.py [per_group=300] > minimap2_pair_oracle.jsonl      (needs minimap2 2.26 on PATH)
"""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from minimap2_oracle_lib import Gen, fasta, fastq, revcomp, emit

g = Gen(20261222)
rnd = g.rnd


REF_POOL = {}


def make_ref():
    key = rnd.randrange(45)
    if key not in REF_POOL:
        REF_POOL[key] = make_ref_new()
    return list(REF_POOL[key])


def make_ref_new():
    ln = rnd.randint(3000, 14000)
    ref = g.rseq(ln, gc=rnd.choice([None, None, 0.65, 0.35]))
    refs = [ref]
    if rnd.random() < 0.35:                       # a repeat: second copy elsewhere (secondary hits), sometimes diverged
        a = rnd.randint(0, ln - 1500)
        rep = ref[a:a + rnd.randint(300, 1200)]
        refs.append(g.rseq(rnd.randint(500, 3000)) + g.mutate(rep, rnd.choice([0.0, 0.01, 0.05]), big_indels=False) + g.rseq(500))
    if rnd.random() < 0.2:
        refs.append(g.rseq(rnd.randint(1000, 5000)))
    return refs


def fragment(refs):
    """(read1, read2) of one fragment: FR pair with an insert 100 - 700, errors, optional adapter / N / chimeric / unmapped mate."""
    ref = refs[rnd.choice([0, 0, 0] + list(range(len(refs))))]
    rl1, rl2 = rnd.choice([100, 150, 150, 250]), rnd.choice([100, 150, 150, 250])
    insert = rnd.choice([rnd.randint(rl1 // 2, 300), rnd.randint(200, 700), rnd.randint(300, 1500)])
    a = rnd.randint(0, max(0, len(ref) - insert - 1))
    frag = ref[a:a + insert]
    rate = rnd.choice([0.0, 0.005, 0.01, 0.03, 0.06])
    r1 = g.mutate(frag[:rl1], rate, big_indels=False)
    r2 = g.mutate(revcomp(frag)[:rl2], rate, big_indels=False)
    x = rnd.random()
    if x < 0.08:
        r2 = g.rseq(rl2)                                        # an unmapped mate
    elif x < 0.12:
        r1 = g.rseq(rl1)
    elif x < 0.20:
        r2 = revcomp(r2)                                        # a pair with the wrong orientation (FF)
    elif x < 0.26:
        other = refs[rnd.randrange(len(refs))]
        b = rnd.randint(0, max(0, len(other) - rl2 - 1))
        r2 = revcomp(other[b:b + rl2])                          # a mate elsewhere
    elif x < 0.30:
        r1 = revcomp(r1)
        r2 = revcomp(r2)
    if rnd.random() < 0.1:
        r1 = g.with_n(r1, rnd.randint(1, 3))
    if rnd.random() < 0.1:
        r2 = r2 + g.rseq(rnd.randint(5, 30))                    # adapter read-through
    if rnd.random() < 0.05:
        r1 = g.rseq(rnd.randint(5, 25)) + r1
    return r1, r2


def pair_case(group, out, kind='files', sam=False, c=True, preset='sr', extra=None, extra_n=(0, 0)):
    for _ in range(40):
        refs = make_ref()
        nfrag = rnd.choice([1, 1, 2, 3])
        names = ['frag%d' % rnd.randint(0, 99999) for _ in range(nfrag)]
        reads = [fragment(refs) for _ in range(nfrag)]
        fq = rnd.random() < 0.5
        suffix = rnd.random() < 0.3
        n1 = [n + ('/1' if suffix else '') for n in names]
        n2 = [n + ('/2' if suffix else '') for n in names]
        files = {'r.fa': fasta(['ref%d' % i for i in range(len(refs))], refs, width=rnd.choice([None, 70]))}
        if fq:
            mk = lambda ns, rs: fastq([(n, r, g.quals(len(r))) for n, r in zip(ns, rs)])
        else:
            mk = lambda ns, rs: fasta(ns, rs)
        args = ['-x', preset] if preset else []
        if kind == 'files':
            files['q1.fa'] = mk(n1, [r[0] for r in reads])
            files['q2.fa'] = mk(n2, [r[1] for r in reads])
            qfiles = ['q1.fa', 'q2.fa']
        else:
            inter = []
            for i in range(nfrag):
                inter.append((n1[i], reads[i][0]))
                inter.append((n2[i], reads[i][1]))
                if kind == 'inter3' and rnd.random() < 0.5:
                    inter.append((n2[i], g.rseq(100)))
            files['q.fa'] = mk([x[0] for x in inter], [x[1] for x in inter])
            qfiles = ['q.fa']
            if not preset:
                args = ['--frag=yes']
        if sam:
            args.append('-a')
        elif c:
            args.append('-c')
        if not sam and c and rnd.random() < 0.4:
            args.append(rnd.choice(['--cs', '--MD', '--eqx', '--cs=long']))
        if sam and rnd.random() < 0.3:
            args.append(rnd.choice(['--MD', '--cs', '--eqx', '-Y', '--secondary-seq', '--sam-hit-only', '-Q', '-y']))
        if extra:
            for f in rnd.sample(extra, rnd.randint(*extra_n)):
                args += f()
        if emit(group, args + ['r.fa'] + qfiles, files, out):
            return True
    return False


PAIR_OPTS = [
    lambda: ['-F', str(rnd.choice([200, 400, 1000, 2000]))],
    lambda: ['-N', str(rnd.choice([1, 3, 10, 30]))],
    lambda: ['-p', str(rnd.choice([0.1, 0.3, 0.8, 1.0]))],
    lambda: ['-k', str(rnd.choice([15, 17, 25]))],
    lambda: ['-w', str(rnd.choice([5, 8, 15]))],
    lambda: ['-A', str(rnd.choice([1, 2, 3])), '-B', str(rnd.choice([2, 4, 8, 10]))],
    lambda: ['-O', str(rnd.choice(['6,20', '12,24', '4,40'])), '-E', str(rnd.choice(['2,1', '3,1', '1,0']))],
    lambda: ['-s', str(rnd.choice([20, 40, 80, 200]))],
    lambda: ['-z', str(rnd.choice(['50', '100,50', '200']))],
    lambda: ['--end-bonus', str(rnd.choice([0, 5, 20]))],
    lambda: ['--max-clip-ratio', str(rnd.choice([0.1, 0.5]))],
    lambda: ['--secondary=yes'],
    lambda: ['--heap-sort=no'],
    lambda: ['--no-pairing'],
    lambda: ['-r', str(rnd.choice(['50', '200', '300,1000']))],
    lambda: ['-g', str(rnd.choice([50, 200, 1000]))],
    lambda: ['-n', str(rnd.choice([1, 3, 5]))],
    lambda: ['-m', str(rnd.choice([15, 40, 60]))],
    lambda: ['-M', str(rnd.choice([0.2, 0.8]))],
    lambda: ['--mask-level', '0.9'],
    lambda: ['-f', str(rnd.choice(['0.001', '50']))],
    lambda: ['-U', '5,100'],
    lambda: ['-t', '2'],
    lambda: ['-K', '1k'],
    lambda: ['-I', '10k'],
    lambda: ['-2'],
    lambda: ['--frag=no'],
    lambda: ['-P'],
    lambda: ['--no-end-flt'],
]


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    out = sys.stdout
    plan = [
        ('pe_files_paf', dict(kind='files', c=False)),
        ('pe_files_pafc', dict(kind='files', c=True)),
        ('pe_files_sam', dict(kind='files', sam=True)),
        ('pe_inter_paf', dict(kind='inter', c=True)),
        ('pe_inter_sam', dict(kind='inter3', sam=True)),
        ('pe_ff', dict(kind='inter', preset=None, c=True)),
        ('pe_options', dict(kind='files', c=True, extra=PAIR_OPTS, extra_n=(1, 4))),
    ]
    for group, kw in plan:
        made = 0
        while made < per:
            k = dict(kw)
            if group == 'pe_ff' and rnd.random() < 0.5:
                k.update(sam=True)
            if group == 'pe_options':
                k['kind'] = rnd.choice(['files', 'inter'])
                k['sam'] = rnd.random() < 0.4
            made += pair_case(group, out, **k)


if __name__ == '__main__':
    main()
