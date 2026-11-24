#!/usr/bin/env python3
"""Builds the oracle of the MAFFT 7.505 guide-tree family (B09 WP72, F326 ff.): --treein (loadtree), --pileup, --randomchain (createchain), --large / --memsavetree, --memsavetreex, --stepadd,
--youngestlinkage (compacttree), and the combinations the script allows.  Every case is one run of the installed script `mafft` (Debian 7.505-1); the oracle holds the return code, standard output,
standard error and every file the run wrote.  Every run is repeated; a case whose two runs differ is reported on stderr and dropped (the program is non-deterministic there).

It reuses the infrastructure of WP70 (generate_mafft_command_oracle.py: sets, option runner, TSV writer; mafft_command_oracle_lib.py).

Usage: python3 scripts/alignment/generate_mafft_compacttree_oracle.py [group ...]   (no group = all)   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (gzip'd TSV), MC_THREADS
"""
import gzip, hashlib, json, os, random, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_compacttree_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))
TREE = 'tree.txt'


def fmt_len(rng, x):
    k = rng.random()
    if k < 0.45:
        return '%f' % x
    if k < 0.65:
        return '%.3f' % x
    if k < 0.75:
        return '%d' % int(round(x * 10))
    if k < 0.82:
        return '%.2e' % x
    if k < 0.88:
        return ('%.4f' % x).lstrip('0')
    return '%.6g' % x


def tree_file(rng, n, style=None):
    """A tree in the format of newick2mafft.rb: one line per node, 'i j len0 len1', i < j, the ids are the smaller sequence number of each group."""
    active = list(range(n))
    lines = []
    for _ in range(n - 1):
        a, b = sorted(rng.sample(active, 2))
        active.remove(b)
        l0, l1 = rng.random() * 0.5, rng.random() * 0.5
        if rng.random() < 0.1:
            l0 = 0.0
        sep = rng.choice([' ', ' ', ' ', '  ', '\t'])
        lines.append('%d%s%d%s%s%s%s' % (a + 1, sep, b + 1, sep, fmt_len(rng, l0), sep, fmt_len(rng, l1)))
    text = '\n'.join(lines) + rng.choice(['\n', '\n', '\n', ''])
    if style == 'crlf':
        text = text.replace('\n', '\r\n')
    elif style == 'cr':
        text = text.replace('\n', '\r')
    elif style == 'blank':
        text = text.replace('\n', '\n\n')
    elif style == 'lead':
        text = '\n' + text
    return text


def bad_tree(rng, n):
    good = tree_file(rng, n).split('\n')
    good = [g for g in good if g]
    kind = rng.choice(['short', 'noeof', 'twoints', 'order', 'range', 'text', 'neglen', 'zerolen', 'dup', 'junk', 'floatids', 'empty-lens', 'sci', 'plus'])
    if kind == 'short':
        return '\n'.join(good[:max(0, len(good) - rng.randint(1, 3))]) + '\n'
    if kind == 'noeof':
        return '\n'.join(good[:-1]) + '\n'
    if kind == 'twoints':
        i = rng.randrange(len(good))
        good[i] = ' '.join(good[i].split()[:2])
        return '\n'.join(good) + '\n'
    if kind == 'order':
        i = rng.randrange(len(good))
        p = good[i].split()
        good[i] = ' '.join([p[1], p[0]] + p[2:])
        return '\n'.join(good) + '\n'
    if kind == 'range':
        i = rng.randrange(len(good))
        p = good[i].split()
        p[1] = str(n + rng.randint(1, 5))
        good[i] = ' '.join(p)
        return '\n'.join(good) + '\n'
    if kind == 'text':
        i = rng.randrange(len(good))
        good[i] = rng.choice(['a b c d', '(1,2);', '1 2 x y', '# comment'])
        return '\n'.join(good) + '\n'
    if kind == 'neglen':
        i = rng.randrange(len(good))
        p = good[i].split()
        p[2] = '-%s' % p[2]
        good[i] = ' '.join(p)
        return '\n'.join(good) + '\n'
    if kind == 'zerolen':
        i = rng.randrange(len(good))
        p = good[i].split()
        p[3] = '0'
        good[i] = ' '.join(p)
        return '\n'.join(good) + '\n'
    if kind == 'dup':
        return '\n'.join(good + [good[-1]]) + '\n'
    if kind == 'junk':
        return '\n'.join(good) + '\nfoo bar\n'
    if kind == 'floatids':
        i = rng.randrange(len(good))
        p = good[i].split()
        p[0] = p[0] + '.5'
        good[i] = ' '.join(p)
        return '\n'.join(good) + '\n'
    if kind == 'empty-lens':
        i = rng.randrange(len(good))
        good[i] = ' '.join(good[i].split()[:3])
        return '\n'.join(good) + '\n'
    if kind == 'sci':
        i = rng.randrange(len(good))
        p = good[i].split()
        p[2] = rng.choice(['1e-1', '1E+0', '0x1p-3', 'inf', 'nan', '1e', '.e1'])
        good[i] = ' '.join(p)
        return '\n'.join(good) + '\n'
    i = rng.randrange(len(good))
    p = good[i].split()
    p[2] = '+' + p[2]
    good[i] = ' '.join(p)
    return '\n'.join(good) + '\n'


def n_of(text):
    return text.count('>')


def case_with_tree(rng, gid, k, opts, mix, treefn, anysymbol=None):
    s = rng.choices([m[0] for m in mix], [m[1] for m in mix])[0]
    v = list(rng.choice(opts))
    if anysymbol is True or (anysymbol is None and s == 'odd'):
        if '--anysymbol' not in v:
            v = ['--anysymbol'] + v
    wrap = rng.choice([None, None, 30, 61])
    seqs, text = C.make_input(rng, C.SETS[s], wrap)
    tree = treefn(rng, len(seqs))
    return ('%s-%04d' % (gid, k), gid, v[:] if TREE not in v else v, text, {TREE: tree})


def build(rng, want):
    out = []
    basic = [('small', 14), ('tiny', 3), ('pair', 2), ('gapped', 3), ('odd', 4), ('mid', 1)]
    n = int(210 * SCALE)

    def add(gid, opts, mix=basic, count=None, anysymbol=None):
        if want and gid not in want:
            return
        for (cid, g, v, text, files) in C.cases(gid, opts, mix, count or n, rng, anysymbol):
            if TREE in v:
                files = {TREE: tree_file(rng, n_of(text))}
            out.append((cid, g, v, text, files))

    generic = [[], ['--retree', '1'], ['--retree', '2'], ['--retree', '3'], ['--reorder'], ['--treeout'], ['--nodeout'], ['--clustalout'], ['--nuc'], ['--amino'], ['--op', '2', '--ep', '0.1'], ['--6merpair'], ['--10merpair', '--nuc'],
               ['--nofft'], ['--bl', '80'], ['--unweight'], ['--anysymbol'], ['--thread', '2'], ['--quiet'], ['--retree', '0'], ['--retree', '0', '--treeout'], ['--outnum'], ['--inputorder'], ['--phylipout'], ['--mixedlinkage', '0.5'],
               ['--averagelinkage'], ['--minimumlinkage'], ['--sparsepickup', '3'], ['--scoreout'], ['--oneiteration']]
    mix_big = [('small', 10), ('tiny', 3), ('pair', 2), ('gapped', 2), ('odd', 3), ('mid', 4)]
    small = [('small', 8), ('tiny', 3), ('pair', 3), ('gapped', 2), ('odd', 3)]

    for gid, base in (('pileup', ['--pileup']), ('randomchain', ['--randomchain']), ('large', ['--large']), ('memsavetree', ['--memsavetree']), ('stepadd', ['--stepadd']), ('youngestlinkage', ['--youngestlinkage']), ('memsavetreex', ['--memsavetreex'])):
        extra = [['--randomseed', '3'], ['--randomseed', '0'], ['--randomseed', '12345'], ['--randomseed', '-7']] if gid == 'randomchain' else []
        add(gid, [base + g for g in generic] + [base + e for e in extra] + [base + e + ['--retree', '1'] for e in extra], mix_big)
    # the combinations the script refuses or passes on
    add('treeext-misc', [['--pileup', '--localpair'], ['--randomchain', '--globalpair'], ['--large', '--maxiterate', '2'], ['--pileup', '--maxiterate', '1'], ['--stepadd', '--maxiterate', '3'], ['--large', '--parttree'],
                         ['--pileup', '--treein', TREE], ['--treein', TREE, '--pileup'], ['--randomchain', '--large'], ['--large', '--pileup'], ['--youngestlinkage', '--stepadd'], ['--pileup', '--nuc', '--auto'],
                         ['--large', '--auto'], ['--large', '--distout'], ['--pileup', '--distout'], ['--stepadd', '--treeout', '--distout'], ['--youngestlinkage', '--distout'], ['--large', '--retree', '0', '--distout'],
                         ['--pileup', '--retree', '0', '--distout'], ['--randomchain', '--sparsepickup', '4'], ['--pileup', '--sparsepickup', '2'], ['--large', '--sparsepickup', '3'], ['--youngestlinkage', '--sparsepickup', '5']],
        basic)
    # --memsave (MSalignmm, MSalignmm.c): alg M for every alignment
    ms_opts = [[], ['--retree', '1'], ['--retree', '3'], ['--nofft'], ['--fft'], ['--treeout'], ['--reorder'], ['--anysymbol'], ['--op', '2', '--ep', '0.1'], ['--6merpair'], ['--clustalout'], ['--unweight'], ['--bl', '80'],
               ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--maxiterate', '2', '--nofft'], ['--retree', '0', '--treeout'], ['--thread', '2'], ['--auto'], ['--sparsepickup', '3'], ['--oneiteration'], ['--scoreout'], ['--nomemsave'],
               ['--parttree'], ['--dpparttree'], ['--parttree', '--partsize', '5'], ['--fastapair'] if False else ['--retree', '2']]
    pair_ms = [[p] + o for p in ('--localpair', '--globalpair', '--genafpair') for o in ([], ['--maxiterate', '2'], ['--weighti', '0'], ['--retree', '0'], ['--treeout'], ['--maxiterate', '1000', '--weighti', '0'], ['--nofft'], ['--anysymbol'])]
    other_ms = [['--large'], ['--pileup'], ['--randomchain'], ['--stepadd'], ['--youngestlinkage'], ['--large', '--retree', '1'], ['--treein', TREE], ['--treein', TREE, '--maxiterate', '2'], ['--add', 'new.fa'], ['--add', 'new.fa', '--keeplength'],
                ['--add', 'new.fa', '--retree', '1'], ['--add', 'new.fa', '--mapout'], ['--memsave']]
    for gid, opts, mix in (('memsave', ms_opts, mix_big), ('memsave-pair', pair_ms, small), ('memsave-other', other_ms, small)):
        if want and gid not in want:
            continue
        for k in range(n):
            o = rng.choice(opts)
            v = ['--memsave'] + list(o)
            sname = rng.choices([m[0] for m in mix], [m[1] for m in mix])[0]
            seqs, text = C.make_input(rng, C.SETS[sname], rng.choice([None, None, 30, 61]))
            if sname == 'odd' and '--anysymbol' not in v:
                v = ['--anysymbol'] + v
            files = {}
            if TREE in v:
                files[TREE] = tree_file(rng, len(seqs))
            if 'new.fa' in v:
                aligned, new = C.aligned_family(rng, rng.randint(2, 9))
                text = aligned
                files['new.fa'] = L.fasta(new, ['n%d' % i for i in range(len(new))])
            out.append(('%s-%04d' % (gid, k), gid, v, text, files))
    # --large / --memsavetree / --stepadd with a pair library: nodepair
    if not want or 'large-pair' in want:
        extras = [[], [], ['--treeout'], ['--retree', '1'], ['--retree', '0'], ['--retree', '3'], ['--weighti', '0'], ['--weighti', '1.5'], ['--lhlimit', '3'], ['--lhlimit', '1'], ['--lhlimit', '0'], ['--anysymbol'], ['--unweight'], ['--reorder'],
                  ['--op', '2', '--ep', '0.1'], ['--thread', '2'], ['--bl', '80'], ['--nofft'], ['--treein', TREE], ['--nodeout'], ['--distout'], ['--quiet'], ['--clustalout'], ['--gop', '-3', '--gep', '0.2'], ['--lop', '-3', '--lep', '0.2'], ['--LOP', '-5'],
                  ['--outnum'], ['--mixedlinkage', '0.5'], ['--sparsepickup', '3'], ['--add', 'new.fa'], ['--maxiterate', '2'], ['--oldgenafpair'], ['--lexp', '-0.2'], ['--bunkatsunashi'], ['--leavegappyregion'], ['--retree', '0', '--treeout']]
        for k in range(n):
            base = rng.choice(['--large', '--memsavetree', '--stepadd', '--large', '--large'])
            pairopt = rng.choice(['--localpair', '--globalpair', '--genafpair'])
            e = list(rng.choice(extras))
            v = [base, pairopt] + e
            sname = rng.choices(['small', 'tiny', 'pair', 'gapped', 'odd', 'mid'], [10, 3, 3, 2, 3, 2])[0]
            seqs, text = C.make_input(rng, C.SETS[sname], rng.choice([None, None, 30, 61]))
            if sname == 'odd' and '--anysymbol' not in v:
                v = ['--anysymbol'] + v
            files = {}
            if TREE in v:
                files[TREE] = tree_file(rng, len(seqs))
            if 'new.fa' in v:
                aligned, new = C.aligned_family(rng, rng.randint(2, 9))
                text = aligned
                files['new.fa'] = L.fasta(new, ['n%d' % i for i in range(len(new))])
            out.append(('large-pair-%04d' % k, 'large-pair', v, text, files))
    # --treein with a tree of the user
    if not want or 'treein' in want:
        opts = generic[:-1] + [['--retree', '1', '--reorder'], ['--treeout', '--reorder'], ['--treeout', '--retree', '0'], ['--nodeout', '--retree', '1']]
        for k in range(n):
            style = rng.choice([None, None, None, None, 'crlf', 'cr', 'blank', 'lead'])
            out.append(case_with_tree(rng, 'treein', k, [['--treein', TREE] + list(o) for o in opts], basic, lambda r, m, s=style: tree_file(r, m, s)))
    small = [('small', 8), ('tiny', 3), ('pair', 3), ('gapped', 2), ('odd', 3)]
    iter_opts = [['--maxiterate', '1'], ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--maxiterate', '2', '--treeout'], ['--maxiterate', '2', '--distout'], ['--maxiterate', '3', '--retree', '1'], ['--maxiterate', '2', '--nofft'],
                 ['--maxiterate', '2', '--treeout', '--distout'], ['--maxiterate', '2', '--anysymbol'], ['--maxiterate', '2', '--reorder'], ['--maxiterate', '2', '--nodeout'], ['--maxiterate', '2', '--op', '2', '--ep', '0.1'],
                 ['--maxiterate', '2', '--retree', '0', '--treeout'], ['--maxiterate', '2', '--thread', '3', '--threadit', '0'], ['--maxiterate', '2', '--bl', '80'], ['--maxiterate', '2', '--unweight'], ['--maxiterate', '2', '--fix', '0.5']]
    pair_opts = [[p] + o for p in ('--localpair', '--globalpair', '--genafpair') for o in ([], ['--maxiterate', '2'], ['--treeout'], ['--distout'], ['--nodeout'], ['--retree', '0'], ['--retree', '0', '--treeout'], ['--weighti', '0'],
                                                                                          ['--maxiterate', '1000', '--treeout', '--distout'], ['--anysymbol'], ['--reorder'], ['--op', '2', '--ep', '0.1'], ['--unweight'], ['--thread', '2'])]
    for gid, opts, mix in (('treein-iter', iter_opts, small), ('treein-pair', pair_opts, small)):
        if not want or gid in want:
            for k in range(n):
                style = rng.choice([None, None, None, 'crlf', 'cr', 'blank'])
                out.append(case_with_tree(rng, gid, k, [['--treein', TREE] + list(o) for o in opts], mix, lambda r, m, s=style: tree_file(r, m, s)))
    if not want or 'treein-pair-err' in want:
        for k in range(n):
            o = rng.choice([['--localpair'], ['--globalpair'], ['--genafpair'], ['--localpair', '--maxiterate', '2'], ['--globalpair', '--treeout'], ['--maxiterate', '2'], ['--maxiterate', '1000', '--anysymbol']])
            c = case_with_tree(rng, 'treein-pair-err', k, [['--treein', TREE] + o], small + [('mid', 1)], bad_tree)
            if '--genafpair' in c[2]:
                # OPEN (B09.md Leftovers, WP72): a NaN / infinite branch length with --genafpair makes the program's dynamic matrix NaN and its alignment degenerate; not ported yet, so not embedded
                c = (c[0], c[1], c[2], c[3], {TREE: c[4][TREE].replace('nan', '9').replace('inf', '9')})
            out.append(c)
    if not want or 'treein-err' in want:
        for k in range(n):
            out.append(case_with_tree(rng, 'treein-err', k, [['--treein', TREE], ['--treein', TREE, '--treeout'], ['--treein', TREE, '--retree', '1'], ['--treein', TREE, '--anysymbol'], ['--treein', TREE, '--quiet']],
                                      [('small', 12), ('tiny', 3), ('pair', 3), ('gapped', 2), ('odd', 3), ('mid', 2)], bad_tree))
    if not want or 'treein-keyword' in want:
        words = ['pileup\n', 'shuffle 5\n', 'shuffle\n', 'shuffle -3\n', 'shuffle 99999\n', 'very compact\n', 'stepadd\n', 'youngestlinkage\n', 'compact  20GB\n', 'compact 0.1k\n', 'compact 100\n', 'compact 5 k\n', 'compact 1M\n', 'compact 0.0001G\n',
                 'compact 3T\n', 'compact 2.5g\n', 'compact abc\n', 'nodepair\n', 'auto 5 20\n', 'auto\n', 'auto 1 1\n', 'test 1 2\n', 'compact 1k\n', 'compact 0.3k\n', 'compact 0.2K\n', 'compact 0.05m\n', 'pileup extra text\n',
                 'shuffle 2 trailing\n', 'compactx 1k\n', 'very compactish\n']
        for k in range(n):
            w = rng.choice(words)
            s = rng.choices(['small', 'tiny', 'pair', 'mid', 'odd'], [10, 2, 2, 3, 3])[0]
            seqs, text = C.make_input(rng, C.SETS[s], rng.choice([None, 30]))
            o = rng.choice([[], ['--retree', '1'], ['--treeout'], ['--reorder'], ['--retree', '3'], ['--anysymbol'], ['--quiet'], ['--localpair'], ['--globalpair', '--treeout'], ['--genafpair'], ['--maxiterate', '2'], ['--localpair', '--maxiterate', '2']])
            out.append(('treein-keyword-%04d' % k, 'treein-keyword', ['--treein', TREE] + o + (['--anysymbol'] if s == 'odd' and '--anysymbol' not in o else []), text, {TREE: w + tree_file(rng, len(seqs))}))
    # the guide tree file "auto" / "test" (varpairscore): 100 or more sequences
    if not want or 'treein-auto' in want:
        words = ['auto 1 5\n', 'auto 2 10\n', 'auto 3 20\n', 'auto 0 2\n', 'auto 7 30\n', 'auto -4 8\n', 'auto 1 100\n', 'test 1 2\n', 'test 5 10\n', 'test 3 3\n', 'test 9 40\n', 'test 0 12\n', 'auto 11\n', 'auto 2 1\n']
        for k in range(n):
            nn = rng.randint(100, 125)
            alpha = rng.choice([L.P, L.N])
            if rng.random() < 0.5:
                seqs = [''.join(rng.choice(alpha) for _ in range(rng.randint(8, 24))) for _ in range(nn)]
            else:
                seqs = L.family(rng, alpha, rng.choice([10, 20]), nn)
            text = L.fasta(seqs)
            w = rng.choice(words)
            o = rng.choice([[], [], ['--retree', '1'], ['--treeout'], ['--reorder'], ['--quiet']])
            out.append(('treein-auto-%04d' % k, 'treein-auto', ['--treein', TREE] + o, text, {TREE: w + tree_file(rng, nn)}))
    # the switches of the program to MSalignmm (groups of more than 30 000 / 10 000 columns) and to the memory-saving tree (20 000 sequences: --auto): heavy cases, few of them
    if not want or 'longcols' in want:
        for k, (nseq, alpha, args) in enumerate([(2, L.N, []), (3, L.N, ['--retree', '1']), (3, L.P, ['--retree', '1', '--nofft'])]):
            base = ''.join(rng.choice(alpha) for _ in range(30500 + 300 * k))
            seqs = [L.mut(rng, base, alpha, 0.08, 0.02) for _ in range(nseq)]
            out.append(('longcols-%04d' % k, 'longcols', args, L.fasta(seqs), {}))
    if not want or 'longparttree' in want:
        for k, (nseq, alpha, args) in enumerate([(6, L.N, ['--parttree']), (8, L.N, ['--parttree', '--partsize', '3']), (5, L.N, ['--dpparttree']), (6, L.P, ['--parttree', '--retree', '1'])]):
            base = ''.join(rng.choice(alpha) for _ in range(10100 + 50 * k))
            seqs = [L.mut(rng, base, alpha, 0.08, 0.02) for _ in range(nseq)]
            out.append(('longparttree-%04d' % k, 'longparttree', args, L.fasta(seqs), {}))
    if not want or 'auto20k' in want:
        seqs = [''.join(rng.choice(L.P) for _ in range(rng.randint(8, 14))) for _ in range(20000)]
        out.append(('auto20k-0000', 'auto20k', ['--auto'], L.fasta(seqs), {}))
    return out


def main():
    rng = random.Random(7250515)
    want = set(sys.argv[1:])
    allcases = build(rng, want)
    n = len(allcases)
    dropped = []
    done = 0
    tsv = os.environ.get('MC_TSV')
    ft = open(tsv, 'w', encoding='utf-8', newline='\n') if tsv else None
    full = 0
    with ThreadPoolExecutor(THREADS) as ex, open(OUT, 'w') as fo:
        for res, bad in ex.map(C.run_case, allcases):
            done += 1
            if res is None:
                dropped.append(bad)
                continue
            fo.write(json.dumps(res) + '\n')
            if ft:
                sha = '-' if res['nostderr'] else hashlib.sha256(res['err'].encode('utf-8')).hexdigest()
                keep = (not res['nostderr']) and len(res['err']) < 200000 and (full % 25 == 0 or len(res['err']) < 400)
                full += 1
                cols = [res['id'], res['group'], json.dumps(res['args']), res['input'], json.dumps(res['files']), str(res['rc']), res['out'], sha,
                        res['err'] if keep else '', json.dumps(res['outfiles'])]
                ft.write('\t'.join(L.esc(c) for c in cols) + '\n')
            if done % 200 == 0:
                print(done, '/', n, file=sys.stderr)
    if ft:
        ft.close()
    print('cases', n, 'written', n - len(dropped), 'non-deterministic (dropped)', len(dropped), dropped[:10], file=sys.stderr)


if __name__ == '__main__':
    main()
