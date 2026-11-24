#!/usr/bin/env python3
"""Builds the oracle of the options of the script mafft that WP73 of batch B09 ported (F334 ff.): one group per option family, every case is one run of the installed script `mafft` (Debian 7.505-1), repeated twice
(a case whose runs differ is dropped: the program is non-deterministic there, or it stops with a message that names a process id); the oracle holds return code, standard output, standard error and every file the run wrote.

Groups: nan (a guide tree with NaN / infinite branch lengths: the shared trace-back matrix commonIP), core-anysymbol (--core --anysymbol with a forced --nuc / --amino), add-pairlib / add-large (--add / --addprofile with a pair
library: tbfast -I; with --large / --memsavetree / --stepadd: nodepair + tbfast -I), smoothing (--add --smoothing:
restorecommongapssmoothly, insertnewgaps_bothorders).

It reuses the infrastructure of WP70 - WP72 (generate_mafft_command_oracle.py: sets, option runner; generate_mafft_compacttree_oracle.py: guide trees; mafft_command_oracle_lib.py: runner, writer).

Usage: python3 scripts/alignment/generate_mafft_wp73_oracle.py [group ...]   (no group = all)   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (TSV, gzip when it ends in .gz), MC_THREADS
"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C
import generate_mafft_compacttree_oracle as T

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_wp73_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))
TREE = T.TREE


def poison(rng, tree_text):
    """Replaces 1 - 3 branch lengths of a tree file by NaN / infinity (the spellings that strtod reads)."""
    lines = tree_text.split('\n')
    idx = [i for i, l in enumerate(lines) if l.strip()]
    for _ in range(rng.choice([1, 1, 2, 3])):
        i = rng.choice(idx)
        parts = lines[i].split()
        parts[rng.choice([2, 3])] = rng.choice(['nan', 'inf', '-nan', 'NaN', '-inf', 'nan', 'inf', 'Infinity', 'NAN(0x12)'])
        lines[i] = ' '.join(parts)
    return '\n'.join(lines)


def build(rng, want):
    out = []
    n = int(330 * SCALE)
    basic = [('small', 10), ('tiny', 3), ('pair', 3), ('gapped', 3), ('odd', 3), ('mid', 1)]
    if not want or 'nan' in want:
        opts = [[], [], [], ['--localpair'], ['--globalpair'], ['--genafpair'], ['--oldgenafpair'], ['--retree', '1'], ['--retree', '2'], ['--nofft'], ['--retree', '1', '--nofft'], ['--maxiterate', '2'], ['--localpair', '--maxiterate', '2'],
                ['--genafpair', '--maxiterate', '2'], ['--treeout'], ['--nuc'], ['--amino'], ['--weighti', '0'], ['--localpair', '--weighti', '0'], ['--genafpair', '--weighti', '0'], ['--reorder'], ['--thread', '2'],
                ['--op', '2', '--ep', '0.1'], ['--6merpair'], ['--unweight'], ['--localpair', '--treeout'], ['--globalpair', '--maxiterate', '2']]
        for k in range(n):
            s = rng.choices([m[0] for m in basic], [m[1] for m in basic])[0]
            seqs, text = C.make_input(rng, C.SETS[s], rng.choice([None, None, 30]))
            tree = poison(rng, T.tree_file(rng, len(seqs)))
            v = ['--treein', TREE] + list(rng.choice(opts))
            if s == 'odd' and '--anysymbol' not in v:
                v.append('--anysymbol')
            out.append(('nan-%04d' % k, 'nan', v, text, {TREE: tree}))
    if not want or 'core-anysymbol' in want:
        bases = [['--maxiterate', '2'], ['--distout'], ['--localpair', '--maxiterate', '2'], ['--maxiterate', '1000', '--nofft'], ['--genafpair', '--maxiterate', '2'], ['--globalpair', '--distout'], ['--maxiterate', '3', '--retree', '1']]
        extra = [[], [], ['--coreext'], ['--corewin', '5'], ['--corethr', '0.3'], ['--clustalout'], ['--reorder'], ['--outnum']]
        for k in range(int(240 * SCALE)):
            s = rng.choices(['small', 'tiny', 'odd', 'odd', 'pair', 'gapped', 'mid'], [6, 2, 5, 5, 2, 2, 2])[0]
            seqs, text = C.make_input(rng, C.SETS[s], None)
            v = ['--core', '--anysymbol', rng.choice(['--nuc', '--amino'])] + rng.choice(bases) + rng.choice(extra)
            out.append(('core-anysymbol-%04d' % k, 'core-anysymbol', v, text, {}))
    if not want or 'add-pairlib' in want:
        pairs = [['--localpair'], ['--globalpair'], ['--genafpair'], ['--oldgenafpair'], ['--localpair', '--maxiterate', '2'], ['--genafpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '1000'],
                 ['--localpair', '--weighti', '0'], ['--globalpair', '--weighti', '0'], ['--localpair', '--weighti', '0', '--maxiterate', '2'], ['--localpair', '--weighti', '1.5'], ['--genafpair', '--retree', '1']]
        mods = [[], [], [], ['--keeplength'], ['--mapout'], ['--compactmapout'], ['--keeplength', '--anysymbol'], ['--reorder'], ['--anysymbol'], ['--retree', '1'], ['--op', '2', '--ep', '0.1'], ['--clustalout'],
                ['--nofft'], ['--treeout'], ['--distout'], ['--thread', '2'], ['--outnum'], ['--unweight'], ['--bl', '80'], ['--minimumweight', '0.1'], ['--6merpair'], ['--lop', '-3', '--lep', '0.2'], ['--quiet']]
        for k in range(int(330 * SCALE)):
            profile = rng.random() < 0.25
            aligned, new = C.aligned_family(rng, rng.randint(2, 9))
            if profile:
                # the added sequences are an alignment themselves (usually no group of the guide tree: the stop of --addprofile)
                new_text = C.aligned_family(rng, rng.randint(2, 5))[0]
            else:
                new_text = L.fasta(new, ['n%d' % i for i in range(len(new))])
            v = ['--addprofile' if profile else '--add', 'new.fa'] + rng.choice(pairs) + rng.choice(mods)
            if '--anysymbol' in v and rng.random() < 0.5 and not profile:
                new_text = L.fasta([new[0][:3] + 'X' + new[0][3:]] + new[1:], ['n%d' % i for i in range(len(new))])
            out.append(('add-pairlib-%04d' % k, 'add-pairlib', v, aligned, {'new.fa': new_text}))
    if not want or 'add-large' in want:
        pairs = [['--localpair'], ['--globalpair'], ['--genafpair'], ['--localpair', '--maxiterate', '2'], ['--genafpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '1000'], ['--localpair', '--weighti', '0', '--maxiterate', '2'],
                 ['--oldgenafpair'], ['--localpair', '--lhlimit', '3'], ['--globalpair', '--lhlimit', '1']]
        trees = [['--large'], ['--large'], ['--memsavetree'], ['--stepadd'], ['--large', '--retree', '1'], ['--youngestlinkage'], ['--large', '--memsavetreex'], ['--stepadd', '--retree', '0']]
        mods = [[], [], [], ['--keeplength'], ['--mapout'], ['--anysymbol'], ['--op', '2', '--ep', '0.1'], ['--treeout'], ['--reorder'], ['--clustalout']]
        for k in range(int(330 * SCALE)):
            profile = rng.random() < 0.2
            aligned, new = C.aligned_family(rng, rng.randint(2, 9))
            new_text = C.aligned_family(rng, rng.randint(2, 5))[0] if profile else L.fasta(new, ['n%d' % i for i in range(len(new))])
            v = ['--addprofile' if profile else '--add', 'new.fa'] + rng.choice(pairs) + rng.choice(trees) + rng.choice(mods)
            out.append(('add-large-%04d' % k, 'add-large', v, aligned, {'new.fa': new_text}))
    if not want or 'smoothing' in want:
        bases = [[], [], [], ['--localpair'], ['--globalpair'], ['--genafpair'], ['--maxiterate', '2'], ['--retree', '1'], ['--nofft'], ['--keeplength'], ['--mapout'], ['--localpair', '--keeplength'], ['--anysymbol'], ['--reorder'],
                 ['--localpair', '--weighti', '0'], ['--large'], ['--localpair', '--large'], ['--op', '2', '--ep', '0.1'], ['--clustalout'], ['--thread', '2']]
        for k in range(int(330 * SCALE)):
            profile = rng.random() < 0.1
            aligned, new = C.aligned_family(rng, rng.randint(2, 12))
            new_text = C.aligned_family(rng, rng.randint(2, 5))[0] if profile else L.fasta(new, ['n%d' % i for i in range(len(new))])
            v = ['--addprofile' if profile else '--add', 'new.fa', '--smoothing'] + rng.choice(bases)
            out.append(('smoothing-%04d' % k, 'smoothing', v, aligned, {'new.fa': new_text}))
    return out


def main():
    rng = random.Random(7350515)
    want = set(sys.argv[1:])
    allcases = build(rng, want)
    L.run_all(allcases, C.run_case, OUT, os.environ.get('MC_TSV'), THREADS)


if __name__ == '__main__':
    main()
