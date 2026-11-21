#!/usr/bin/env python3
"""A3-9c2 / A5-9b (B09 WP45) oracle: the options of the mapping / alignment stages of minimap2 2.26 beyond the base-level alignment itself, the WHOLE standard output of the program as text (and, for the error group, the exit status and the standard error
at verbosity 2: `-v 2` suppresses the [M::...] log).

JSON lines like generate_minimap2_align_oracle.py (error cases also have "stderr" and "exit", `-o` cases "outfiles").
Groups (default 300 cases each, seed 20261223):
  opt_index_parts - references of 2 - 8 sequences with `-I` (the index is built and queried in parts of ~NUM bases: the output of every part follows the previous one, mid_occ of the first part is used by the others, SAM without @SQ lines) and
                    `--split-prefix` (hits of all parts merged: mm_split_merge / merge_hits), PAF / SAM, `-c`, single reads and pairs
  opt_chain       - -c with the chaining / seeding options not covered elsewhere: --max-chain-iter --chain-gap-scale --chain-skip-scale -g -G -F -f -U -e --q-occ-frac -T -H --for-only --rev-only --max-qlen --seed --no-hash-name -D -X -P
                    --dual= --sr --heap-sort= --rmq= --no-long-join --max-chain-skip -r
  opt_output      - --paf-no-hit --sam-hit-only -y -Q -Y -L --secondary-seq -R --eqx --cs --MD --secondary= -o FILE (the output goes to the file), unmapped reads, FASTA / FASTQ, comments
  opt_alt         - --alt LIST and --alt-drop with reference contigs marked ALT
  opt_qstrand     - --qstrand (-c, --cs, --MD)
  opt_ava         - all-vs-all mappings with -c (-x ava-ont / ava-pb, -X, -D, --dual=no)
  opt_input       - unusual input: lowercase / IUPAC / U bases, empty and very short reads, reads of N, multi-line FASTA / FASTQ, CRLF, blank lines, comments, an empty reference sequence
  opt_long_cigar  - one case of a 610 kb query with more than 65 535 CIGAR operations (-L: the CG tag)
  opt_errors      - invalid or special command lines: option checks of mm_check_opt, unknown / missing option arguments, unknown preset, --splice with --frag, -V, -h, missing input, `-N 0`, deprecated options (exit status, stderr and stdout)
Usage: python3 scripts/alignment/generate_minimap2_options_oracle.py [per_group=300] > minimap2_options_oracle.jsonl      (needs minimap2 2.26 on PATH)
"""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
import generate_minimap2_align_oracle as A
from minimap2_oracle_lib import fasta, fastq, revcomp, emit

A.g.rnd.seed(20261223)
g = A.g
rnd = g.rnd


MULTI_POOL = {}


def multi_ref(n=None):
    key = rnd.randrange(50)
    if key not in MULTI_POOL:
        MULTI_POOL[key] = multi_ref_new(n)
    return list(MULTI_POOL[key])


def multi_ref_new(n=None):
    n = n or rnd.randint(2, 6)
    refs = []
    for _ in range(n):
        refs.append(g.rseq(rnd.choice([0, 50, 300, 1500, 4000, 7000]) if rnd.random() < 0.15 else rnd.randint(500, 4500)))
    refs = [r for r in refs if r] or [g.rseq(3000)]
    return refs


def reads_from(refs, n, rate=0.04, fq=False, prefix='rd'):
    out = []
    for k in range(n):
        ref = refs[rnd.randrange(len(refs))]
        if len(ref) < 400:
            q = g.rseq(rnd.randint(100, 300))
        else:
            a = rnd.randint(0, len(ref) - 300)
            q = g.mutate(ref[a:a + rnd.randint(150, min(len(ref) - a, 4000))], rate, big_indels=False)
            if rnd.random() < 0.4:
                q = revcomp(q)
        out.append(('%s%d' % (prefix, k), q))
    return out


def qfile(reads, fq):
    return fastq([(n, q, g.quals(len(q))) for n, q in reads]) if fq else fasta([n for n, _ in reads], [q for _, q in reads])


def case_index_parts(out):
    refs = multi_ref()
    paired = rnd.random() < 0.35
    fq = rnd.random() < 0.3
    files = {'r.fa': fasta(['ref%d' % i for i in range(len(refs))], refs)}
    args = []
    if paired:
        frs = []
        for k in range(rnd.randint(1, 3)):
            ref = refs[rnd.randrange(len(refs))]
            if len(ref) < 800:
                continue
            a = rnd.randint(0, len(ref) - 600)
            frag = ref[a:a + rnd.randint(250, 550)]
            frs.append(('fr%d' % k, g.mutate(frag[:150], 0.01, False), g.mutate(revcomp(frag)[:150], 0.01, False)))
        if not frs:
            return False
        files['q1.fa'] = qfile([(n + '/1', a) for n, a, b in frs], fq)
        files['q2.fa'] = qfile([(n + '/2', b) for n, a, b in frs], fq)
        args = ['-x', 'sr']
        qf = ['q1.fa', 'q2.fa']
    else:
        files['q.fa'] = qfile(reads_from(refs, rnd.randint(1, 4), fq=fq), fq)
        pre = rnd.choice([None, 'map-ont', 'map-hifi', 'asm10', 'sr'])
        args = ['-x', pre] if pre else []
        qf = ['q.fa']
    args += [rnd.choice(['-c', '-a', '-c', '--cs'])]
    args += ['-I', rnd.choice(['1', '1k', '5k', '12k', '30k', '100k', '4g'])] if rnd.random() < 0.8 else []
    if rnd.random() < 0.5:
        args += ['--split-prefix', 'sp']
    if rnd.random() < 0.15:
        args += [rnd.choice(['--eqx', '-N', '2', '--secondary=no', '-P'])] if False else [rnd.choice(['--eqx', '--secondary=no'])]
    if '--cs' in args and '--split-prefix' in args:
        args.remove('--cs'); args.append('-c')
    return emit('opt_index_parts', args + ['-v', '2', 'r.fa'] + qf, files, out, True)


CHAIN_OPTS = [
    lambda: ['--max-chain-iter', str(rnd.choice([1, 3, 20, 200]))],
    lambda: ['--chain-gap-scale', str(rnd.choice([0, 0.3, 1.5, 3]))],
    lambda: ['--chain-skip-scale', str(rnd.choice([0, 0.2, 1.0]))],
    lambda: ['-g', str(rnd.choice([100, 500, 3000, 20000]))],
    lambda: ['-G', str(rnd.choice(['500', '5k', '100k']))],
    lambda: ['-F', str(rnd.choice([100, 1000, 5000]))],
    lambda: ['-f', str(rnd.choice(['0.0001', '0.01', '20', '5,50']))],
    lambda: ['-U', str(rnd.choice(['2,10', '10,50', '5,1000']))],
    lambda: ['-e', str(rnd.choice([0, 100, 2000]))],
    lambda: ['--q-occ-frac', str(rnd.choice([0, 0.001, 0.2]))],
    lambda: ['-T', str(rnd.choice([10, 20, 40]))],
    lambda: ['-H'],
    lambda: ['--for-only'],
    lambda: ['--rev-only'],
    lambda: ['--max-qlen', str(rnd.choice([100, 1000, 3000, 100000]))],
    lambda: ['--seed', str(rnd.choice([0, 1, 11, 12345]))],
    lambda: ['--no-hash-name'],
    lambda: ['-D'],
    lambda: ['-P'],
    lambda: ['--dual=' + rnd.choice(['yes', 'no'])],
    lambda: ['--sr'],
    lambda: ['--heap-sort=' + rnd.choice(['yes', 'no'])],
    lambda: ['--rmq=' + rnd.choice(['yes', 'no'])],
    lambda: ['--no-long-join'],
    lambda: ['--max-chain-skip', str(rnd.choice([2, 10, 50]))],
    lambda: ['-r', rnd.choice(['100', '200,5000', '3000'])],
    lambda: ['--mask-level', str(rnd.choice([0.1, 0.9]))],
    lambda: ['--lj-min-ratio', '0.5'],
    lambda: ['--print-aln-seq'],
    lambda: ['--print-qname'],
    lambda: ['--print-seeds'],
    lambda: ['--print-chains'],
    lambda: ['--min-occ-floor', str(rnd.choice([1, 5, 30, 100]))],
    lambda: ['--splice'],
    lambda: ['-u', rnd.choice('fbrn')],
    lambda: ['-J', rnd.choice('01')],
    lambda: ['--no-kalloc'],
    lambda: ['--cap-kalloc', '1g'],
    lambda: ['--idx-no-seq'],
]


def case_chain(out):
    refs = A.make_target(hp=False)
    q = A.make_query(refs, False, rnd.choice([0.02, 0.06, 0.12]))
    files = {'r.fa': fasta(['ref%d' % i for i in range(len(refs))], refs), 'q.fa': fasta(['q%d' % rnd.randint(0, 999)], [q])}
    pre = rnd.choice([None, None, 'map-ont', 'map-hifi', 'asm10', 'map-pb'])
    args = (['-x', pre] if pre else []) + [rnd.choice(['-c', '-c', '-a'])]
    for f in rnd.sample(CHAIN_OPTS, rnd.randint(1, 4)):
        args += f()
    return emit('opt_chain', args + ['-v', '2', 'r.fa', 'q.fa'], files, out, True)


def case_output(out):
    refs = A.make_target(hp=False)
    n = rnd.randint(1, 4)
    reads = []
    for k in range(n):
        if rnd.random() < 0.2:
            reads.append(('junk%d' % k, g.rseq(rnd.randint(30, 400))))
        else:
            reads.append(('read%d' % k, A.make_query(refs, False, 0.06)))
    fq = rnd.random() < 0.5
    txt = qfile(reads, fq)
    if rnd.random() < 0.4:
        txt = txt.replace('read0', 'read0 comment:Z:xy\tNM:i:7 more words', 1)
    files = {'r.fa': fasta(['ref%d' % i for i in range(len(refs))], refs), 'q.fa': txt}
    args = rnd.choice([['-a'], ['-c'], ['-a'], ['--cs'], ['-x', 'map-hifi', '-a'], ['-x', 'map-pb', '-c'], ['-a', '-x', 'map-ont']])
    for f in [['--paf-no-hit'], ['--sam-hit-only'], ['-y'], ['-Q'], ['-Y'], ['-L'], ['--secondary-seq'], ['-R', '@RG\\tID:grp\\tSM:smp'], ['--eqx'], ['--cs'], ['--cs=long'], ['--MD'], ['--secondary=no'],
              ['-S'], ['-o', 'out.txt'], ['-o', '-'], ['-t', '4'], ['-K', '500k'], ['-v', '1'], ['-N', '0']]:
        if rnd.random() < 0.18:
            args += f
    if rnd.random() < 0.15:
        args += [rnd.choice(['-R', '@RG\\tSM:x']), ] if False else ['-R', rnd.choice(['@RG\\tSM:x', 'RG\\tID:x', '@RG\\tID:a\\tSM:b\\tPL:c'])]
    if '-v' not in args:
        args = ['-v', '2'] + args
    return emit('opt_output', args + ['r.fa', 'q.fa'], files, out, True, ('out.txt',))


ALT_POOL = {}


def case_alt(out):
    key = rnd.randrange(45)
    if key not in ALT_POOL:
        ref = A.make_target()[0]
        ALT_POOL[key] = (ref, g.mutate(ref, rnd.choice([0.0, 0.002, 0.01, 0.03]), big_indels=False), g.mutate(ref, 0.01, False))
    ref, alt1, alt2 = ALT_POOL[key]
    names = ['chr1', 'chr1_alt']
    seqs = [ref, alt1]
    if rnd.random() < 0.4:
        names.append('chr1_alt2'); seqs.append(alt2)
    a = rnd.randint(0, max(0, len(ref) - 3000))
    q = g.mutate(ref[a:a + rnd.randint(600, 3000)], rnd.choice([0.01, 0.05]), big_indels=False)
    if rnd.random() < 0.4:
        q = revcomp(q)
    files = {'r.fa': fasta(names, seqs), 'q.fa': fasta(['q1'], [q]), 'alt.txt': '\n'.join(n + '\t' + str(len(s)) for n, s in zip(names, seqs) if 'alt' in n) + '\n'}
    args = rnd.choice([['-c'], ['-a'], ['-x', 'map-ont', '-c'], ['-x', 'sr', '-a'], ['-x', 'map-pb', '-a']])
    args += ['--alt', 'alt.txt']
    if rnd.random() < 0.6:
        args += ['--alt-drop', str(rnd.choice([0, 0.05, 0.15, 0.5, 1.0]))]
    if rnd.random() < 0.3:
        args += ['-N', '10', '-p', '0.5']
    return emit('opt_alt', args + ['r.fa', 'q.fa'], files, out)


def case_qstrand(out):
    refs = A.make_target()
    q = A.make_query(refs, False, rnd.choice([0.02, 0.08]))
    files = {'r.fa': fasta(['ref%d' % i for i in range(len(refs))], refs), 'q.fa': fasta(['q1'], [q])}
    args = rnd.choice([[], ['-x', 'map-hifi'], ['-x', 'asm10']]) + ['-c', '--qstrand']
    if rnd.random() < 0.5:
        args.append(rnd.choice(['--cs', '--MD', '--eqx', '--cs=long']))
    if rnd.random() < 0.3:
        args += rnd.choice([['-N', '4'], ['-P'], ['--secondary=no']])
    return emit('opt_qstrand', args + ['r.fa', 'q.fa'], files, out)


AVA_POOL = {}


def case_ava(out):
    key = rnd.randrange(45)
    if key not in AVA_POOL:
        base = g.rseq(rnd.randint(2000, 5000))
        reads = []
        for i in range(rnd.randint(3, 6)):
            a = rnd.randint(0, len(base) // 3)
            b = rnd.randint(2 * len(base) // 3, len(base))
            r = g.mutate(base[a:b], rnd.choice([0.02, 0.08, 0.12]), big_indels=False)
            reads.append(('read%d' % i, revcomp(r) if rnd.random() < 0.4 else r))
        AVA_POOL[key] = reads
    reads = AVA_POOL[key]
    files = {'r.fa': fasta([n for n, _ in reads], [r for _, r in reads])}
    qi = rnd.randrange(len(reads))
    files['q.fa'] = fasta([reads[qi][0]], [reads[qi][1]])
    pre = rnd.choice([['-x', 'ava-ont'], ['-x', 'ava-pb'], ['-X', '-x', 'map-ont'], ['-D', '-x', 'map-ont'], ['-x', 'map-ont', '--dual=no'], ['-X']])
    args = pre + [rnd.choice(['-c', '-c', '-a', '--cs'])]
    if rnd.random() < 0.3:
        args.append(rnd.choice(['--eqx', '--MD']))
    return emit('opt_ava', args + ['r.fa', 'q.fa'], files, out)


IUPAC = 'RYKMSWBDHVN'


def noisy(seq):
    """lowercase stretches, IUPAC codes, U instead of T"""
    c = list(seq)
    for _ in range(rnd.randint(0, 6)):
        i = rnd.randrange(len(c))
        c[i] = rnd.choice(IUPAC)
    if rnd.random() < 0.4:
        a = rnd.randrange(len(c))
        b = min(len(c), a + rnd.randint(5, 200))
        c[a:b] = [x.lower() for x in c[a:b]]
    s = ''.join(c)
    if rnd.random() < 0.15:
        s = s.replace('T', 'U').replace('t', 'u')
    return s


def case_input(out):
    refs = A.make_target()
    names = ['ref%d' % i for i in range(len(refs))]
    seqs = [noisy(r) if rnd.random() < 0.5 else r for r in refs]
    reads = []
    for k in range(rnd.randint(1, 5)):
        x = rnd.random()
        if x < 0.12:
            q = ''
        elif x < 0.22:
            q = g.rseq(rnd.randint(1, 14))
        elif x < 0.28:
            q = 'N' * rnd.randint(20, 200)
        else:
            q = noisy(A.make_query(refs, False, 0.04))
        reads.append(('r%d' % k, q))
    fq = rnd.random() < 0.4
    if fq:
        txt = ''.join('@%s\n%s\n+\n%s\n' % (n, q, g.quals(len(q))) for n, q in reads)
    else:
        txt = ''
        for n, q in reads:
            txt += '>%s\n' % n
            w = rnd.choice([0, 60, 70, 7])
            if w and q:
                txt += ''.join(q[i:i + w] + '\n' for i in range(0, len(q), w))
            else:
                txt += q + '\n'
    if rnd.random() < 0.15:
        txt = txt.replace('\n', '\r\n')
    if rnd.random() < 0.15:
        txt = '\n' + txt + '\n\n'
    ref_txt = fasta(names, seqs)
    if rnd.random() < 0.1:
        ref_txt += '>emptyref\n\n'
    if rnd.random() < 0.1:
        ref_txt = ref_txt.replace('>ref0', '>ref0 some description here', 1)
    files = {'r.fa': ref_txt, 'q.fa': txt}
    args = rnd.choice([['-a'], ['-c'], ['-a', '--cs'], ['-c', '--MD'], ['-x', 'sr', '-a'], ['-x', 'map-pb', '-c'], ['-a', '--paf-no-hit'], ['-c', '--paf-no-hit'], ['-a', '--sam-hit-only']])
    return emit('opt_input', ['-v', '2'] + args + ['r.fa', 'q.fa'], files, out, True)


def case_long_cigar(out, k):
    """A 610 kb query at ~17 % substitutions aligned with --eqx: more than 65 535 CIGAR operations (-L writes them to the CG tag, -Y changes the clipping)."""
    ref = g.rseq(640000)
    q = ''.join(c if rnd.random() > 0.17 else rnd.choice('ACGT') for c in ref[10000:620000])
    files = {'r.fa': fasta(['big'], [ref]), 'q.fa': fasta(['bigq'], [q])}
    args = [['-a', '--eqx', '-L'], ['-a', '-L', '--eqx', '-Y'], ['-c', '--eqx', '-L', '--cs']][k % 3]
    return emit('opt_long_cigar', ['-v', '2'] + args + ['r.fa', 'q.fa'], files, out, True)


ERRORS = [
    ['-k', '0'], ['-w', '0'], ['-k', '-3'], ['-N', '-1'], ['-N', '0'], ['-p', '1.5'], ['-p', '-0.1'], ['-r', '2000,100'], ['--rmq', '--sr'], ['--rmq', '--splice'], ['--rmq=yes', '-x', 'splice'],
    ['--for-only', '--rev-only'], ['-O', '0'], ['-E', '0'], ['-O', '4,3'], ['-E', '2,2'], ['-O', '60,70', '-E', '5,1'], ['-z', '100,300'], ['-X', '--secondary=no'], ['-P', '--secondary=no'], ['--qstrand', '-a'],
    ['--qstrand', '-H'], ['--qstrand', '--splice'], ['--qstrand', '--frag=yes'], ['--splice', '--frag=yes'], ['-x', 'nonexistent'], ['-x'], ['-k'], ['--unknown'], ['-Z'], ['-u', 'q'], ['-u', 'f'], ['--cs=bogus', '-c'],
    ['--cs=none', '-c'], ['--cs=long'], ['--frag=maybe'], ['--secondary=maybe'], ['--dual=maybe'], ['--heap-sort=x'], ['--rmq=x'], ['--splice-flank=x'], ['-V'], ['--version'], ['-h'], ['--help'], ['-x', 'sr', '-h'], ['-x', 'asm5', '-h'],
    ['-x', 'splice', '-h'], ['-S'], ['--lj-min-ratio', '0.3'], ['-a', '-R', '@RG\\tSM:x'], ['-a', '-R', 'bad'], ['-a', '-R', '@RG\\tID:q'], ['--MD', '--split-prefix', 'x'], ['--cs', '--split-prefix', 'x'], ['-A', '300'], ['-B', '200'],
    ['-x', 'sr', '-F', '300'], ['-N', '3', '-N', '0'], ['-N', '0', '-N', '7'], ['-v', '0', '-k', '0'], ['-v', '1', '-N', '0'], ['--sr', '--rmq'], ['-J', '0', '-x', 'splice', '-c'], ['-J', '2'], ['-C', '5', '-x', 'splice', '-c'],
    ['-t', 'abc'], ['-K', '10zz'], ['--max-intron-len', '1000', '-x', 'splice'], ['-G', '1k'], ['--no-pairing', '-x', 'sr'], ['--bucket-bits', '8'], ['-d', 'idx.mmi'], ['-I', '0'], ['-I', '1k', '-a'], ['-c', '-a'],
]


def case_errors(out):
    refs = [g.rseq(rnd.randint(1500, 3000))]
    q = g.mutate(refs[0][100:900], 0.03, False)
    files = {'r.fa': fasta(['ref0'], refs), 'q.fa': fasta(['q1'], [q]), 'q2.fa': fasta(['q2'], [revcomp(q)]), 'q3.fa': fasta(['q3'], [q])}
    args = list(rnd.choice(ERRORS))
    if rnd.random() < 0.15:
        args = ['-x', rnd.choice(['sr', 'map-ont', 'splice'])] + args
    if rnd.random() < 0.1:
        args = args + ['-x', rnd.choice(['sr', 'asm5'])]
    pos = rnd.choice([['r.fa', 'q.fa'], ['r.fa', 'q.fa'], ['r.fa'], [], ['missing.fa', 'q.fa'], ['r.fa', 'missing.fa'], ['r.fa', 'q.fa', 'q2.fa'], ['r.fa', 'q.fa', 'q2.fa', 'q3.fa']])
    if '-v' not in args:
        args = ['-v', '2'] + args
    return emit('opt_errors', args + pos, files, out, True)


def main():
    per = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    out = sys.stdout
    case_long_cigar(out, 0)
    for name, fn in [('opt_index_parts', case_index_parts), ('opt_chain', case_chain), ('opt_output', case_output), ('opt_alt', case_alt), ('opt_qstrand', case_qstrand), ('opt_ava', case_ava), ('opt_input', case_input), ('opt_errors', case_errors)]:
        made = 0
        tries = 0
        while made < per and tries < per * 8:
            tries += 1
            made += bool(fn(out))


if __name__ == '__main__':
    main()
