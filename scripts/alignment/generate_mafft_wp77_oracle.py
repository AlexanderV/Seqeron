#!/usr/bin/env python3
"""Builds the oracle of `mafft --sp` (WP77 of batch B09: the programs D__align and D__align_variousdist of Dalignmm.c): every case is one run of the installed script `mafft` (Debian 7.505-1), repeated twice (a case whose runs
differ is dropped); the oracle holds return code, standard output, standard error and every file the run wrote.

Groups (>= 200 sets each, scale with MC_SCALE):
  sp          --sp alone (disttbfast -@): DNA and protein families of 2 - 12 sequences of 20 - 400 residues, the output options of the script
  sp-opt      --sp with the scoring options: --op / --ep / --bl / --jtt / --tm / --kimura / --aamatrix / --fmodel / --leavegappyregion ...
  sp-iter     --sp with --maxiterate (dvtditr -@ -B -Z -z 1000: Falign with D__align per segment) and --skipiterate / --core / --nofft
  sp-pair     --sp with a pair library (--localpair / --globalpair / --genafpair / --oldgenafpair, with and without --maxiterate: tbfast and dvtditr with the importance matrix in D__align)
  sp-tree     --sp with the guide tree options: --retree n, --treein, --parttree, --pileup, --randomchain, --large, --memsavetree, --stepadd, --6merpair ..., --thread, --oneiteration
  sp-unalign  --sp with --unalignlevel x / --allowshift / --shiftpenalty x (D__align_variousdist and the shift of D__align)
  sp-add      --sp with --add / --addprofile / --merge / --addfragments / --addfull / --addlong / --addtotop / --seed
Usage: python3 scripts/alignment/generate_mafft_wp77_oracle.py [group ...]   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (TSV, gzip when it ends in .gz), MC_THREADS
"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C
import generate_mafft_compacttree_oracle as T
import generate_mafft_wp74_oracle as W74
import generate_mafft_wp75_oracle as W75

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_wp77_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))
TREE = T.TREE
PROT = 'ACDEFGHIKLMNPQRSTVWY'


def fam(rng, dna=None, nmin=2, nmax=12, lmin=20, lmax=400):
    dna = rng.random() < 0.45 if dna is None else dna
    alpha = 'ACGT' if dna else PROT
    n = rng.choice([2, 3, 4, 5, 6, 8, 10, 12]) if nmax >= 12 else rng.randint(nmin, nmax)
    n = max(nmin, min(n, nmax))
    length = rng.choice([20, 35, 60, 100, 160, 250, 400] if lmax >= 400 else [20, 35, 60, 100, 160])
    base = ''.join(rng.choice(alpha) for _ in range(length))
    r = rng.choice([0.03, 0.1, 0.2, 0.35])
    ind = rng.choice([0.0, 0.02, 0.05, 0.1])
    seqs = [L.mut(rng, base, alpha, r, ind) for _ in range(n)]
    if rng.random() < 0.15:
        # a fragment and a long insertion: gaps of several lengths in the columns (the case D__align is made for)
        i = rng.randrange(len(seqs))
        a = rng.randrange(len(seqs[i]) // 2 + 1)
        seqs[i] = seqs[i][a:a + rng.randint(8, max(9, len(seqs[i]) // 2))]
    if rng.random() < 0.12:
        i = rng.randrange(len(seqs))
        p = rng.randrange(len(seqs[i]) + 1)
        seqs[i] = seqs[i][:p] + ''.join(rng.choice(alpha) for _ in range(rng.randint(3, 25))) + seqs[i][p:]
    return seqs, dna


def make(rng, dna=None, **kw):
    seqs, dna = fam(rng, dna, **kw)
    names = C.names_for(rng, len(seqs))
    return seqs, dna, L.fasta(seqs, names, rng.choice([None, None, None, 30, 61]))


OUTPUT_OPTS = [[], [], [], [], ['--reorder'], ['--clustalout'], ['--phylipout'], ['--treeout'], ['--distout'], ['--treeout', '--distout'], ['--inputorder'], ['--outnum'], ['--quiet'], ['--anysymbol'], ['--nuc'], ['--amino'], ['--scoreout'], ['--noscore'],
               ['--preservecase'], ['--namelength', '12', '--clustalout'], ['--reorder', '--treeout']]
SCORE_OPTS = [['--op', '1.0'], ['--op', '2.5', '--ep', '0.2'], ['--ep', '0.3'], ['--op', '0.5'], ['--op', '3.0', '--ep', '0.05'], ['--bl', '30'], ['--bl', '45'], ['--bl', '62'], ['--bl', '80'], ['--jtt', '20'], ['--jtt', '100'], ['--tm', '20'], ['--tm', '100'],
              ['--kimura', '1'], ['--kimura', '200'], ['--fmodel'], ['--leavegappyregion'], ['--nwildcard'], ['--op', '1.5', '--bl', '45'], ['--op', '2.0', '--ep', '0.1', '--jtt', '50'], ['--nuc', '--op', '2.0'], ['--nuc', '--kimura', '2'],
              ['--amino', '--tm', '60'], ['--lop', '-1.0'], ['--gop', '1.5'], ['--termgappenalty', '0.5'], ['--exp', '0.3'], ['--opdist', '1.0'], ['--amax', '0.2'], ['--amax', '0.5']]
LETTERS = 'ARNDCQEGHILKMFPSTWYV'


def matrix_text(rng):
    order = list(LETTERS)
    if rng.random() < 0.7:
        rng.shuffle(order)
    lines = ['   ' + ' '.join(order)]
    for i in range(20):
        vals = [str(rng.randint(-4, 3) if i != j else rng.randint(4, 11)) for j in range(i + 1)]
        lines.append(order[i] + ' ' + ' '.join(vals))
    return '\n'.join(lines) + '\n'


def sp_case(rng, k, group, extras_fn, pre=None, dna=None):
    seqs, d, text = make(rng, dna)
    args, files = ['--sp'], {}
    args += extras_fn(rng, d, files, seqs)
    if '--anysymbol' not in args and rng.random() < 0.0:
        args.append('--anysymbol')
    return ('%s-%04d' % (group, k), group, args, text, files)


def x_sp(rng, d, files, seqs):
    return list(rng.choice(OUTPUT_OPTS))


def x_opt(rng, d, files, seqs):
    v = list(rng.choice(SCORE_OPTS))
    if rng.random() < 0.15:
        v += rng.choice(OUTPUT_OPTS)
    if rng.random() < 0.12 and not d:
        files['mtx.txt'] = matrix_text(rng)
        v = ['--aamatrix', 'mtx.txt'] + v
    return v


ITER = [['--maxiterate', '1'], ['--maxiterate', '2'], ['--maxiterate', '3'], ['--maxiterate', '5'], ['--maxiterate', '10'], ['--maxiterate', '1000'], ['--maxiterate', '2', '--nofft'], ['--maxiterate', '2', '--fft'], ['--maxiterate', '2', '--retree', '1'],
        ['--maxiterate', '3', '--retree', '2'], ['--maxiterate', '2', '--skipiterate', '0.3'], ['--maxiterate', '2', '--core'], ['--maxiterate', '2', '--core', '--coreext'], ['--maxiterate', '2', '--bunkatsunashi'], ['--maxiterate', '2', '--op', '2', '--ep', '0.2'],
        ['--maxiterate', '2', '--bl', '45'], ['--maxiterate', '2', '--jtt', '100'], ['--maxiterate', '2', '--tm', '100'], ['--maxiterate', '2', '--treeout'], ['--maxiterate', '2', '--reorder'], ['--maxiterate', '2', '--clustalout'],
        ['--maxiterate', '2', '--threadit', '0'], ['--maxiterate', '2', '--thread', '-1', '--threadit', '0'], ['--maxiterate', '2', '--thread', '2'], ['--maxiterate', '2', '--minimumweight', '0.1'], ['--maxiterate', '2', '--unweight'],
        ['--maxiterate', '2', '--6merpair'], ['--maxiterate', '2', '--mixedlinkage', '0.5'], ['--maxiterate', '2', '--averagelinkage'], ['--maxiterate', '2', '--anysymbol'], ['--maxiterate', '2', '--quiet'], ['--maxiterate', '2', '--randomseed', '3']]


def x_iter(rng, d, files, seqs):
    v = list(rng.choice(ITER))
    if rng.random() < 0.2:
        v += rng.choice(OUTPUT_OPTS)
    return v


PAIRS = [['--localpair'], ['--globalpair'], ['--genafpair'], ['--oldgenafpair'], ['--localpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '2'], ['--genafpair', '--maxiterate', '2'], ['--localpair', '--maxiterate', '1000'],
         ['--localpair', '--weighti', '0'], ['--globalpair', '--weighti', '1.5'], ['--localpair', '--lop', '-1.5', '--lep', '0.2'], ['--localpair', '--lexp', '-0.05'], ['--genafpair', '--gop', '-2.0', '--gep', '0.1', '--gexp', '-0.2'],
         ['--localpair', '--op', '2.0', '--ep', '0.2'], ['--localpair', '--retree', '1'], ['--localpair', '--nofft'], ['--globalpair', '--treeout'], ['--localpair', '--reorder'], ['--localpair', '--clustalout'], ['--localpair', '--thread', '2'],
         ['--localpair', '--maxiterate', '2', '--threadit', '0'], ['--localpair', '--bl', '45'], ['--localpair', '--jtt', '100'], ['--localpair', '--minimumweight', '0.1'], ['--localpair', '--unweight'], ['--localpair', '--anysymbol'], ['--localpair', '--quiet']]


def x_pair(rng, d, files, seqs):
    return list(rng.choice(PAIRS))


TREEOPT = [['--retree', '0'], ['--retree', '1'], ['--retree', '2'], ['--retree', '3'], ['--retree', '1', '--nofft'], ['--nofft'], ['--fft'], ['--6merpair'], ['--10merpair'], ['--6mermultipair'], ['--10mermultipair'], ['--parttree'], ['--parttree', '--partsize', '5'],
           ['--dpparttree'], ['--pileup'], ['--randomchain'], ['--randomchain', '--randomseed', '4'], ['--large'], ['--memsavetree'], ['--stepadd'], ['--youngestlinkage'], ['--thread', '2'], ['--thread', '-1'],
           ['--threadtb', '2'], ['--oneiteration'], ['--mixedlinkage', '0.5'], ['--averagelinkage'], ['--minimumlinkage'], ['--sparsepickup', '3'], ['--memsave'], ['--nomemsave'], ['--auto'], ['--treein', 'tree.txt'], ['--treein', 'tree.txt', '--maxiterate', '2'],
           ['--treein', 'tree.txt', '--retree', '1'], ['--treein', 'tree.txt', '--localpair'], ['--treein', 'tree.txt', '--treeout'], ['--pileup', '--maxiterate', '2'], ['--randomchain', '--localpair'], ['--large', '--localpair']]


def x_tree(rng, d, files, seqs):
    v = list(rng.choice(TREEOPT))
    if '--treein' in v:
        files['tree.txt'] = T.tree_file(rng, len(seqs))
    return v


UNALIGN = [['--globalpair', '--unalignlevel', '0.1'], ['--globalpair', '--unalignlevel', '0.3'], ['--globalpair', '--unalignlevel', '0.5'], ['--globalpair', '--unalignlevel', '0.8'], ['--globalpair', '--allowshift'], ['--globalpair', '--shiftpenalty', '2'],
           ['--globalpair', '--shiftpenalty', '5'], ['--globalpair', '--shiftpenalty', '50'], ['--globalpair', '--unalignlevel', '0.4', '--allowshift'], ['--globalpair', '--unalignlevel', '0.4', '--maxiterate', '2'], ['--globalpair', '--allowshift', '--maxiterate', '2'],
           ['--globalpair', '--unalignlevel', '0.2', '--allowshift', '--maxiterate', '3'], ['--globalpair', '--unalignlevel', '0.6', '--retree', '1'], ['--globalpair', '--unalignlevel', '0.3', '--treeout'], ['--globalpair', '--unalignlevel', '0.3', '--weighti', '1.5'],
           ['--globalpair', '--unalignlevel', '0.3', '--op', '2.0'], ['--globalpair', '--allowshift', '--reorder'], ['--localpair', '--unalignlevel', '0.3'], ['--genafpair', '--allowshift'], ['--unalignlevel', '0.3'], ['--allowshift'], ['--shiftpenalty', '3'],
           ['--globalpair', '--unalignlevel', '0.3', '--maxiterate', '2', '--threadit', '0'], ['--globalpair', '--unalignlevel', '0.3', '--addfragments', 'new.fa'], ['--globalpair', '--weighti', '0', '--allowshift'], ['--globalpair', '--amax', '0.3'], ['--amax', '0.3']]


def x_unalign(rng, d, files, seqs):
    v = list(rng.choice(UNALIGN))
    if 'new.fa' in v:
        files['new.fa'] = L.fasta([s[:rng.randint(10, max(11, len(s)))] for s in seqs[:2]], ['f0', 'f1'])
    return v


def add_case(rng, k):
    """--sp with the add family: an existing alignment (made by the installed mafft) plus new sequences / subalignments / seeds."""
    dna = rng.random() < 0.4
    alpha = 'ACGT' if dna else PROT
    length = rng.choice([30, 60, 120, 200])
    base = ''.join(rng.choice(alpha) for _ in range(length))
    r = rng.choice([0.05, 0.15, 0.3])
    ind = rng.choice([0.02, 0.05, 0.1])
    old = [L.mut(rng, base, alpha, r, ind) for _ in range(rng.randint(2, 8))]
    new = [L.mut(rng, base, alpha, r, ind) for _ in range(rng.randint(1, 5))]
    rc, aln, err, files0 = L.run_once(['--quiet'], L.fasta(old, ['a%d' % i for i in range(len(old))]))
    kind = rng.choice(['--add', '--add', '--addprofile', '--merge', '--addfragments', '--addfull', '--addlong', '--addtotop', '--seed', '--merge', '--add'])
    files = {}
    ex = list(rng.choice([[], [], [], ['--retree', '1'], ['--maxiterate', '2'], ['--localpair'], ['--6merpair'], ['--nofft'], ['--keeplength'], ['--reorder'], ['--op', '2', '--ep', '0.2'], ['--treeout'], ['--thread', '2'], ['--anysymbol'],
                                       ['--globalpair', '--unalignlevel', '0.3'], ['--clustalout'], ['--mapout'], ['--bl', '45']]))
    if kind == '--seed' and '--thread' in ex:
        ex = []   # the iterative refinement of --seed with threads: D__align keeps static work arrays, the threads would share them
    text = aln
    if kind in ('--add', '--addfragments', '--addfull', '--addlong', '--addtotop'):
        if kind in ('--addfragments', '--addtotop'):
            new = [s[rng.randrange(max(1, len(s) // 3)):] for s in new]
        files['new.fa'] = L.fasta(new, ['n%d' % i for i in range(len(new))])
        args = [kind, 'new.fa']
    elif kind == '--addprofile':
        rc2, aln2, _, _ = L.run_once(['--quiet'], L.fasta(new + [L.mut(rng, base, alpha, r, ind)], ['p%d' % i for i in range(len(new) + 1)]))
        files['new.fa'] = aln2 if rc2 == 0 else L.fasta(new)
        args = [kind, 'new.fa']
    elif kind == '--merge':
        rc2, aln2, _, _ = L.run_once(['--quiet'], L.fasta(new + [L.mut(rng, base, alpha, r, ind)], ['p%d' % i for i in range(len(new) + 1)]))
        rows1 = [''.join(x.split('\n')[1:]) for x in aln.split('>')[1:]]
        rows2 = [''.join(x.split('\n')[1:]) for x in aln2.split('>')[1:]] if rc2 == 0 else new
        n1 = len(rows1)
        text = L.fasta(rows1 + rows2, ['a%d' % i for i in range(n1)] + ['p%d' % i for i in range(len(rows2))])
        files['table.txt'] = ' '.join(str(i + 1) for i in range(n1)) + '\n' + ' '.join(str(n1 + i + 1) for i in range(len(rows2))) + '\n'
        args = [kind, 'table.txt']
    else:
        # --seed: an alignment file as a seed of the sequences of the input
        files['seed.fa'] = aln
        text = L.fasta(new + [L.mut(rng, base, alpha, r, ind)], ['n%d' % i for i in range(len(new) + 1)])
        args = [kind, 'seed.fa']
    return ('sp-add-%04d' % k, 'sp-add', ['--sp'] + args + ex, text, files)


def build(rng, want):
    out = []

    def add(group, fn, n):
        if not want or group in want:
            for k in range(n):
                out.append(fn(k))

    n = int(210 * SCALE)
    add('sp', lambda k: sp_case(rng, k, 'sp', x_sp), n)
    add('sp-opt', lambda k: sp_case(rng, k, 'sp-opt', x_opt), n)
    add('sp-iter', lambda k: sp_case(rng, k, 'sp-iter', x_iter), int(n * 1.12))      # a few runs end with a segmentation fault of dvtditr (use of a freed block with -S): not recorded
    add('sp-pair', lambda k: sp_case(rng, k, 'sp-pair', x_pair), n)
    add('sp-tree', lambda k: sp_case(rng, k, 'sp-tree', x_tree), n)
    add('sp-unalign', lambda k: sp_case(rng, k, 'sp-unalign', x_unalign), n)
    add('sp-add', lambda k: add_case(rng, k), int(n * 1.4))      # --sp with --add: D__align reads mtx[-1] for the new gap character '=' (segmentation fault of disttbfast): not recorded
    return out


def run_case(case):
    res, bad = C.run_case(case)
    if res is not None and 'Segmentation fault' in res['err']:
        return None, case[0]
    return res, bad


def main():
    rng = random.Random(7700117)
    want = set(sys.argv[1:])
    allcases = build(rng, want)
    L.run_all(allcases, run_case, OUT, os.environ.get('MC_TSV'), THREADS)


if __name__ == '__main__':
    main()
