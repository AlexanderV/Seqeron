#!/usr/bin/env python3
"""Builds the oracle of `mafft --hybridpair` without a fragment option (WP78 of batch B09: the chain pairlocalalign -Y -> disttbfast -E 1 -T -y -> tbfast of mafft.tmpl 2615 - 2624): every case is one run of the installed
script `mafft` (Debian 7.505-1), repeated twice (a case whose runs differ is dropped); the oracle holds return code, standard output, standard error and every file the run wrote.

Groups (>= 200 sets each, scale with MC_SCALE):
  hybrid        --hybridpair alone: DNA and protein families of 2 - 12 sequences of 20 - 400 residues, tiny / odd / gapped sets, the output options of the script
  hybrid-opt    --hybridpair with the scoring options (--op / --ep / --bl / --jtt / --tm / --kimura / --aamatrix / --lop / --lep / --lexp / --weighti ...)
  hybrid-iter   --hybridpair with --maxiterate (dvtditr with the distances of disttbfast and the pair library of pairlocalalign -Y) and --skipiterate / --core / --nofft
  hybrid-tree   --hybridpair with the guide tree options: --retree n, --treein, --pileup, --randomchain, --large, --memsavetree, --stepadd, --parttree, --6merpair (the last pair option wins), --thread, --treeout, --distout, --nodeout
  hybrid-add    --hybridpair with --add / --addprofile / --keeplength / --mapout / --merge / --seed / --anchors / --adjustdirection / --sp / --unalignlevel / --text
Usage: python3 scripts/alignment/generate_mafft_wp78_oracle.py [group ...]   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (TSV, gzip when it ends in .gz), MC_THREADS
"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C
import generate_mafft_compacttree_oracle as T
import generate_mafft_wp77_oracle as W77
import generate_mafft_wp76_oracle as W76

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_wp78_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))
HY = '--hybridpair'


def any_input(rng, dna=None):
    """A family (as of WP77) or one of the odd sets of the command oracle (tiny, gapped, odd letters)."""
    x = rng.random()
    if x < 0.7:
        seqs, d, text = W77.make(rng, dna)
        return seqs, d, text
    gen = rng.choice([L.tiny_set, L.odd_set, L.gapped_set, L.small_set, L.pair_set])
    seqs = gen(rng)
    names = C.names_for(rng, len(seqs))
    text = L.fasta(seqs, names, rng.choice([None, None, 30, 61]))
    dna = all(set(s.upper()) <= set('ACGTN-') for s in seqs)
    return seqs, dna, text


def case(rng, k, group, extras_fn, dna=None):
    seqs, d, text = any_input(rng, dna)
    files = {}
    args = [HY] + extras_fn(rng, d, files, seqs)
    return ('%s-%04d' % (group, k), group, args, text, files)


def x_plain(rng, d, files, seqs):
    v = list(rng.choice(W77.OUTPUT_OPTS))
    if rng.random() < 0.15:
        v += rng.choice([['--namelength', '20', '--clustalout'], ['--phylipout', '--reorder'], ['--op', '2'], ['--quiet', '--anysymbol']])
    return v


def x_opt(rng, d, files, seqs):
    v = list(rng.choice(W77.SCORE_OPTS + [['--lop', '-3', '--lep', '0.2'], ['--lexp', '-0.2'], ['--weighti', '1.5'], ['--weighti', '0'], ['--weighti', '5'], ['--LOP', '-5', '--LEXP', '-0.1'], ['--gop', '-3'],
                                          ['--similaritylevel', '0.1'], ['--similaritylevel', '-0.2'], ['--nzero'], ['--unweight'], ['--tuplesize', '4'], ['--6merpair', '--op', '2'], ['--legacygappenalty']]))
    if rng.random() < 0.15:
        v += rng.choice(W77.OUTPUT_OPTS)
    if rng.random() < 0.12 and not d:
        files['mtx.txt'] = W77.matrix_text(rng)
        v = ['--aamatrix', 'mtx.txt'] + v
    return v


def x_iter(rng, d, files, seqs):
    v = list(rng.choice(W77.ITER))
    if rng.random() < 0.2:
        v += rng.choice(W77.OUTPUT_OPTS)
    return v


def x_tree(rng, d, files, seqs):
    v = list(rng.choice(W77.TREEOPT + [['--nodeout'], ['--retree', '0', '--treeout'], ['--retree', '0', '--distout'], ['--retree', '0', '--nodeout'], ['--retree', '2'], ['--memsavetreex'], ['--thread', '4'], ['--threadtb', '3'],
                                         ['--treein', 'tree.txt', '--distout'], ['--treein', 'tree.txt', '--nodeout'], ['--partsize', '4', '--parttree']]))
    if '--treein' in v:
        files['tree.txt'] = T.tree_file(rng, len(seqs))
    return v


def add_case(rng, k):
    """--hybridpair with the add family (existing alignment + new sequences / subalignments / seeds) and the options that the script allows only with it."""
    c = W77.add_case(rng, k)
    cid, group, args, text, files = c
    args = [HY if a == '--sp' else a for a in args]
    ex = rng.choice([[], [], [], ['--keeplength'], ['--mapout'], ['--compactmapout'], ['--smoothing'], ['--keeplength', '--maxiterate', '2'], ['--sp'], ['--unalignlevel', '0.3'], ['--maxiterate', '2']])
    if any(a in args for a in ('--addfragments', '--addfull', '--addlong', '--addtotop', '--seed', '--merge')):
        ex = [e for e in ex if e not in ('--keeplength', '--mapout', '--compactmapout', '--smoothing')] if rng.random() < 0.7 else ex
    return ('hybrid-add-%04d' % k, 'hybrid-add', args + ex, text, files)


def misc_case(rng, k):
    """The rest of the script's own options that act on the hybrid chain."""
    pick = rng.choice(['sp', 'unalign', 'text', 'anchors', 'direction', 'core', 'multi', 'auto', 'pairs'])
    if pick == 'anchors':
        cid, group, args, text, files = W76.anchor_case(rng, k)
        return ('hybrid-add-%04d' % k, 'hybrid-add', [HY] + args, text, files)
    seqs, d, text = any_input(rng)
    files = {}
    if pick == 'sp':
        v = ['--sp'] + rng.choice([[], ['--maxiterate', '2'], ['--localpair']])
    elif pick == 'unalign':
        v = rng.choice([['--unalignlevel', '0.3'], ['--unalignlevel', '0.6', '--allowshift'], ['--allowshift'], ['--shiftpenalty', '3'], ['--unalignlevel', '0.4', '--maxiterate', '2']])
    elif pick == 'text':
        v = ['--text'] + rng.choice([[], ['--maxiterate', '2']])
    elif pick == 'direction':
        v = rng.choice([['--adjustdirection'], ['--adjustdirectionaccurately']])
    elif pick == 'core':
        v = ['--core'] + rng.choice([[], ['--coreext'], ['--corethr', '0.5']])
    elif pick == 'multi':
        v = rng.choice([['--multipair'], ['--6mermultipair'], ['--hybridpair'], ['--10mermultipair'], ['--6merpair']])
        v = [a for a in v if a != HY]
    elif pick == 'auto':
        v = ['--auto']
    else:
        v = rng.choice([['--localpair'], ['--globalpair'], ['--genafpair'], ['--6merpair'], ['--10merpair'], ['--parttree'], ['--dpparttree']])
    args = [HY] + v if rng.random() < 0.6 else v + [HY]
    return ('hybrid-add-%04d' % k, 'hybrid-add', args, text, files)


def build(rng, want):
    out = []

    def add(group, fn, n):
        if not want or group in want:
            for k in range(n):
                out.append(fn(k))

    n = int(235 * SCALE)      # the runs that end with a segmentation fault of the program (--scoreout, threads without --add) are not recorded
    add('hybrid', lambda k: case(rng, k, 'hybrid', x_plain), n)
    add('hybrid-opt', lambda k: case(rng, k, 'hybrid-opt', x_opt), n)
    add('hybrid-iter', lambda k: case(rng, k, 'hybrid-iter', x_iter), n)
    add('hybrid-tree', lambda k: case(rng, k, 'hybrid-tree', x_tree), n)
    add('hybrid-add', lambda k: add_case(rng, k) if k % 3 else misc_case(rng, k), int(n * 1.3))
    return out


def run_case(case):
    res, bad = C.run_case(case)
    if res is not None and 'Segmentation fault' in res['err']:
        return None, case[0]
    return res, bad


def main():
    rng = random.Random(7800117)
    want = set(sys.argv[1:])
    allcases = build(rng, want)
    L.run_all(allcases, run_case, OUT, os.environ.get('MC_TSV'), THREADS)


if __name__ == '__main__':
    main()
