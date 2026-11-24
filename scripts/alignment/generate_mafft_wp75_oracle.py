#!/usr/bin/env python3
"""Builds the oracle of the fragment options of the script mafft that WP75 of batch B09 registered (F350: the program `addsingle`): every case is one run of the installed script `mafft` (Debian 7.505-1), repeated twice
(a case whose runs differ is dropped: the program is non-deterministic there); the oracle holds return code, standard output, standard error and every file the run wrote. Needs no fasta34 / lastal (the installed `mafft` runs these options).

Groups (>= 200 sets each, scale with MC_SCALE):
  frag-6mer       --addfragments / --addfull / --addlong / --addtotop with --6merpair / --10merpair (addsingle -W n: the k-tuple distances of the program itself), protein and DNA, with --keeplength, --mapout, --anysymbol, --reorder,
                  --treeout --distout, --op --ep, --retree 1, --nofft, --bl 80, --thread 2, --smoothing-free variants and the stops of the script
  frag-6mermulti  the same with --6mermultipair / --10mermultipair (addsingle -d -W n: the distances of the new sequences to the old ones and of the old ones among themselves are computed separately)
  frag-multi      --addfragments ... without a pair option and with --multipair (pairlocalalign -Y + dndpre -M 2 + addsingle -d)
  frag-hybrid     --hybridpair (pairlocalalign -Y, disttbfast -E 1 -T -y, addsingle)
  frag-local      --localpair with fragments (pairlocalalign, addsingle with hat3)
  frag-codon      --codonpos with --6merpair --addfragments (DNA of codons); --codonscore with a score file; the script's stops for the other pair options
  frag-treein     --treein tree.txt (the guide tree of the old sequences: the program loadtop, lengths of -1 are recomputed) with every pair option
  frag-opt        every pair option with one or two further options of the script (--compactmapout, --smoothing, --op / --ep, --kimura, --jtt, --tm, --bl, --retree, --thread, --mixedlinkage, --lop / --lep / --lexp, --weighti, --phylipout, --maxiterate ...)
Input: the file in.fa is an alignment made by the installed mafft from a family, new.fa holds fragments (pieces of family members with mutations, 8 - 70 residues), long members or the full-length members.
Usage: python3 scripts/alignment/generate_mafft_wp75_oracle.py [group ...]   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (TSV, gzip when it ends in .gz), MC_THREADS
"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_wp75_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))
FRAG_KINDS = ['--addfragments', '--addfragments', '--addfull', '--addlong', '--addtotop']
EXTRAS = [[], [], [], ['--keeplength'], ['--mapout'], ['--anysymbol'], ['--reorder'], ['--treeout', '--distout'], ['--op', '2', '--ep', '0.1'], ['--retree', '1'], ['--nofft'], ['--bl', '80'], ['--thread', '2'],
          ['--keeplength', '--mapout'], ['--clustalout'], ['--quiet'], ['--leavegappyregion'], ['--maxambiguous', '0.5'], ['--outnum'], ['--unweight'], ['--minimumweight', '0.1'], ['--nuc'], ['--amino'], ['--fmodel']]


def family(rng, alpha, k, length):
    base = ''.join(rng.choice(alpha) for _ in range(length))
    r = rng.choice([0.03, 0.1, 0.25])
    ind = rng.choice([0.0, 0.02, 0.06])
    return base, [L.mut(rng, base, alpha, r, ind) for _ in range(k)], r, ind


def make_case(rng, kind, dna=None, codon=False):
    """The alignment (in.fa), the fragments / long / full sequences (new.fa) for the option `kind`."""
    dna = rng.random() < 0.4 if dna is None else dna
    alpha = 'ACGT' if dna else 'ACDEFGHIKLMNPQRSTVWY'
    length = rng.choice([30, 60, 120, 240]) if not codon else 3 * rng.choice([12, 25, 40])
    base, fam, r, ind = family(rng, alpha, 1 if codon else rng.randint(2, 9), length)     # --codonpos: the program supports a reference MSA of ONE sequence only ("For now, use a single sequence as reference")
    rc, aln, err, files = L.run_once(['--quiet'], L.fasta(fam, ['a%d' % i for i in range(len(fam))]))
    nnew = rng.randint(1, 6)
    new = []
    for _ in range(nnew):
        if kind in ('--addfragments', '--addtotop') and rng.random() < 0.85:
            a = rng.randrange(max(1, len(base) - 8))
            b = min(len(base), a + rng.randint(8, max(9, min(70, len(base)))))
            if codon:
                a -= a % 3
                b -= b % 3
                b = max(b, a + 6)
            piece = base[a:b]
        else:
            piece = base
        new.append(L.mut(rng, piece, alpha, r, ind) if not codon else piece)
    if rng.random() < 0.1:
        new.append(new[0])
    if rng.random() < 0.08:
        new[0] = new[0].lower()
    if rng.random() < 0.08:
        new[-1] = new[-1][:2] + ('N' if dna else 'X') + new[-1][2:]
    return aln, new


CODONS = [a + b + c for a in 'acgt' for b in 'acgt' for c in 'acgt']
GENETIC_CODE = 'KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF'


def codonscore_file(rng):
    lines = ['# codon scores']
    for i, c1 in enumerate(CODONS):
        for j in range(i, 64):
            c2 = CODONS[j]
            a1, a2 = GENETIC_CODE[i], GENETIC_CODE[j]
            score = 6.0 if i == j else (2.0 if a1 == a2 else -1.0 + 0.25 * (i % 3 == j % 3))
            lines.append('%s %s %s %s %.2f' % (c1, a1, c2, a2, score))
    return '\n'.join(lines) + '\n'


def codonpos_file(aln):
    """One line per column of the first sequence of the alignment: the codon position (1, 2, 3 on the residues of the sequence, 0 on its gaps)."""
    row0 = ''.join(aln.split('>')[1].split('\n')[1:])
    pos, out = 0, []
    for ch in row0:
        if ch == '-':
            out.append('0')
        else:
            out.append(str(pos % 3 + 1))
            pos += 1
    return '\n'.join(out) + '\n'


def case(rng, group, k, pairopts, kinds=FRAG_KINDS, extras=EXTRAS, **kw):
    kind = rng.choice(kinds)
    aln, new = make_case(rng, kind, **kw)
    args = [kind, 'new.fa'] + list(pairopts) + list(rng.choice(extras))
    if rng.random() < 0.03:
        args = [kind, 'new.fa'] + list(rng.choice([['--6merpair', '--add', 'new.fa'], ['--6merpair', '--merge', 'x'], ['--6merpair', '--retree', '0'], ['--6merpair', '--maxiterate', '2'], ['--6merpair', '--treein', 'nofile'], ['--6merpair', '--seed', 'x']]))
    files = {'new.fa': L.fasta(new, ['f%d' % i for i in range(len(new))])}
    if '--codonpos' in args:
        args = [('codonpos.txt' if a == '--codonpos' and False else a) for a in args]
        i = args.index('--codonpos')
        args.insert(i + 1, 'codonpos.txt')
        files['codonpos.txt'] = codonpos_file(aln)
        if rng.random() < 0.4:
            args += ['--codonscore', 'codonscore.txt']
            files['codonscore.txt'] = codonscore_file(rng)
    return ('%s-%04d' % (group, k), group, args, aln, files)


OPT_EXTRAS = [['--compactmapout', '--keeplength'], ['--smoothing'], ['--leavegappyregion'], ['--op', '1.2', '--ep', '0.2'], ['--kimura', '1'], ['--jtt', '100'], ['--tm', '100'], ['--bl', '45'], ['--retree', '2'], ['--thread', '-1'],
              ['--mapout', '--keeplength', '--anysymbol'], ['--keeplength', '--outnum'], ['--mixedlinkage', '0.5'], ['--averagelinkage'], ['--minimumlinkage'], ['--phylipout'], ['--clustalout', '--reorder'], ['--noscore'], ['--lop', '-1.5'],
              ['--lep', '0.2'], ['--lexp', '-0.05'], ['--weighti', '1.5'], ['--nuc', '--thread', '3'], ['--inputorder'], ['--treeout', '--reorder'], ['--memsave'], ['--parttree'], ['--maxiterate', '3'], ['--maxambiguous', '0.2'], ['--fft'],
              ['--nofft', '--thread', '2'], ['--6merpair', '--thread', '2', '--keeplength']]
ALL_PAIRS = [['--6merpair'], ['--10merpair'], ['--6mermultipair'], ['--10mermultipair'], ['--multipair'], [], ['--hybridpair'], ['--localpair']]


def treein_case(rng, k):
    """--treein: a rooted bifurcated tree of the old sequences (one line per node: two representatives, two branch lengths; -1 = recompute)."""
    kind = rng.choice(FRAG_KINDS)
    aln, new = make_case(rng, kind)
    n = aln.count('>')
    reps = list(range(1, n + 1))
    lines = []
    while len(reps) > 1:
        i, j = sorted(rng.sample(reps, 2))
        l1 = '%.5f' % (rng.random() * 0.3 + 0.01)
        l2 = '%.5f' % (rng.random() * 0.3 + 0.01) if rng.random() < 0.7 else '-1.0'
        lines.append('%d %d %s %s' % (i, j, l1, l2))
        reps.remove(j)
    args = [kind, 'new.fa', '--treein', 'tree.txt'] + list(rng.choice(ALL_PAIRS)) + (list(rng.choice(OPT_EXTRAS)) if rng.random() < 0.3 else [])
    files = {'new.fa': L.fasta(new, ['f%d' % i for i in range(len(new))]), 'tree.txt': '\n'.join(lines) + '\n'}
    return ('frag-treein-%04d' % k, 'frag-treein', args, aln, files)


def opt_case(rng, k):
    kind = rng.choice(FRAG_KINDS)
    aln, new = make_case(rng, kind)
    ex = list(rng.choice(OPT_EXTRAS)) if rng.random() < 0.95 else []
    if rng.random() < 0.3:
        ex += list(rng.choice(OPT_EXTRAS))
    args = [kind, 'new.fa'] + list(rng.choice(ALL_PAIRS)) + ex
    return ('frag-opt-%04d' % k, 'frag-opt', args, aln, {'new.fa': L.fasta(new, ['f%d' % i for i in range(len(new))])})


def build(rng, want):
    out = []

    def add(group, pairopts, n, **kw):
        if not want or group in want:
            for k in range(n):
                out.append(case(rng, group, k, rng.choice(pairopts) if isinstance(pairopts[0], list) else pairopts, **kw))

    n = int(220 * SCALE)
    add('frag-6mer', [['--6merpair'], ['--6merpair'], ['--10merpair']], n)
    add('frag-6mermulti', [['--6mermultipair'], ['--10mermultipair']], n)
    add('frag-multi', [[], ['--multipair']], n)
    add('frag-hybrid', [['--hybridpair']], n)
    add('frag-local', [['--localpair']], n)
    if not want or 'frag-codon' in want:
        for k in range(int(210 * SCALE)):
            c = case(rng, 'frag-codon', k, ['--6merpair', '--codonpos'], kinds=['--addfragments', '--addfragments', '--addfull', '--addlong'], dna=True, codon=True)
            if rng.random() < 0.15 and '--codonpos' in c[2]:
                args = [x for x in c[2] if x not in ('--6merpair', 'codonpos.txt', '--codonscore', 'codonscore.txt')]
                i = args.index('--codonpos')
                args[i + 1:i + 1] = ['codonpos.txt', rng.choice(['--multipair', '--localpair', '--hybridpair'])]
                c = (c[0], c[1], args, c[3], c[4])
            out.append(c)
    if not want or 'frag-treein' in want:
        for k in range(int(210 * SCALE)):
            out.append(treein_case(rng, k))
    if not want or 'frag-opt' in want:
        for k in range(int(210 * SCALE)):
            out.append(opt_case(rng, k))
    return out


def run_case(case):
    """C.run_case without the runs that end in a segmentation fault of a program (the shell names a process id in its message)."""
    res, bad = C.run_case(case)
    if res is not None and 'Segmentation fault' in res['err']:
        return None, case[0]
    return res, bad


def main():
    rng = random.Random(7500115)
    want = set(sys.argv[1:])
    allcases = build(rng, want)
    L.run_all(allcases, run_case, OUT, os.environ.get('MC_TSV'), THREADS)


if __name__ == '__main__':
    main()
