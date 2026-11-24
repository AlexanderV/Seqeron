#!/usr/bin/env python3
"""Builds the oracle of the MAFFT 7.505 command line (B09 WP70, F315 ff.): every case is one run of the installed script `mafft` (Debian 7.505-1) with an option list on a deterministic input; the oracle holds
the return code, standard output, standard error and every file the run wrote.  Every run is repeated; a case whose two runs differ is reported on stderr and dropped (the program is non-deterministic there).

Usage: python3 scripts/alignment/generate_mafft_command_oracle.py [group ...]    (no group = all)   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_THREADS
Output: JSONL, one case per line: {id, group, args, name, input, files, rc, out, err, outfiles}.  scripts/alignment/embed_mafft_command_oracle.py packs it into the gzip'd TSV the tests read.
"""
import gzip, hashlib, json, os, random, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_command_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))


def names_for(rng, n):
    kinds = [lambda i: 's%d' % i,
             lambda i: 'seq_%d some description %s' % (i, rng.choice(['alpha', 'beta|gamma', 'x=y', 'tab\there'])),
             lambda i: ('N%d_' % i) + ''.join(rng.choice('abcdefghij _-.|:') for _ in range(rng.randint(1, 60))),
             lambda i: ('L%d_' % i) + 'x' * rng.choice([100, 250, 260, 300]),
             lambda i: '  lead%d' % i]
    k = rng.random()
    if k < 0.55:
        return ['s%d' % i for i in range(n)]
    f = rng.choice(kinds)
    return [f(i) for i in range(n)]


def make_input(rng, gen, wrap=None, named=True):
    seqs = gen(rng)
    names = names_for(rng, len(seqs)) if named else None
    return seqs, L.fasta(seqs, names, wrap)


SETS = {
    'small': L.small_set, 'pair': L.pair_set, 'mid': L.mid_set, 'tiny': L.tiny_set, 'odd': L.odd_set, 'gapped': L.gapped_set,
}


def cases(group, opt_variants, mix, count, rng, anysymbol=None):
    """mix: list of (set name, weight)."""
    out = []
    names = [m[0] for m in mix]
    weights = [m[1] for m in mix]
    for k in range(count):
        s = rng.choices(names, weights)[0]
        v = list(rng.choice(opt_variants))
        if anysymbol is True or (anysymbol is None and s == 'odd'):
            if '--anysymbol' not in v:
                v = ['--anysymbol'] + v
        wrap = rng.choice([None, None, 30, 61])
        seqs, text = make_input(rng, SETS[s], wrap)
        out.append(('%s-%04d' % (group, k), group, v, text, {}))
    return out


def groups(rng):
    g = []
    basic = [('small', 14), ('tiny', 3), ('pair', 2), ('gapped', 3), ('odd', 4), ('mid', 1)]
    n = int(210 * SCALE)
    g += cases('format', [[], ['--clustalout'], ['--phylipout'], ['--reorder'], ['--clustalout', '--reorder'], ['--phylipout', '--reorder'], ['--inputorder'],
                         ['--clustalout', '--namelength', '8'], ['--clustalout', '--namelength', '40'], ['--phylipout', '--namelength', '20'], ['--phylipout', '--namelength', '3'],
                         ['--clustalout', '--namelength', '0'], ['--quiet'], ['--quiet', '--clustalout'], ['--clustalout', '--anysymbol'], ['--phylipout', '--anysymbol', '--reorder']], basic, n, rng)
    g += cases('scoring', [['--op', '2.5'], ['--op', '0.5', '--ep', '0.1'], ['--ep', '-0.2'], ['--exp', '0.2'], ['--bl', '80'], ['--bl', '45'], ['--jtt', '100'], ['--tm', '150'], ['--kimura', '300'],
                         ['--kappa', '4'], ['--fmodel'], ['--nwildcard'], ['--nzero'], ['--leavegappyregion'], ['--legacygappenalty'], ['--termgappenalty'], ['--unweight'], ['--nuc'], ['--amino'],
                         ['--6merpair'], ['--10merpair', '--nuc'], ['--mixedlinkage', '0.5'], ['--averagelinkage'], ['--minimumlinkage'], ['--fft'], ['--nofft'], ['--opdist', '2'],
                         ['--similaritylevel', '0.1'], ['--similaritylevel', '-0.2', '--op', '1'], ['--retree', '1'], ['--retree', '3'], ['--retree', '0'], ['--retree', '1', '--treeout']], basic, n, rng)
    g += cases('iterative', [['--maxiterate', '1'], ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--maxiterate', '3', '--retree', '1'], ['--maxiterate', '2', '--nofft'], ['--maxiterate', '2', '--thread', '3', '--threadit', '0']],
               [('small', 14), ('tiny', 3), ('pair', 2), ('gapped', 3), ('odd', 4)], n, rng)
    g += cases('pairlib', [['--localpair'], ['--globalpair'], ['--genafpair'], ['--oldgenafpair'], ['--localpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '1000'], ['--genafpair', '--maxiterate', '3'],
                           ['--localpair', '--weighti', '1.5'], ['--globalpair', '--gop', '-3', '--gep', '0.2'], ['--localpair', '--lop', '-3', '--lep', '0.2', '--lexp', '-0.2'], ['--genafpair', '--LOP', '-5', '--LEXP', '-0.1'],
                           ['--localpair', '--op', '2', '--ep', '0.1'], ['--globalpair', '--bl', '80'], ['--localpair', '--thread', '2']],
               [('small', 8), ('tiny', 3), ('pair', 3), ('gapped', 2), ('odd', 3)], n, rng)
    g += cases('auto', [['--auto'], ['--auto', '--op', '2'], ['--auto', '--clustalout'], ['--auto', '--quiet']], [('small', 14), ('tiny', 3), ('pair', 2), ('odd', 3), ('mid', 4)], n, rng)
    g += cases('parttree', [['--parttree'], ['--dpparttree'], ['--parttree', '--partsize', '5'], ['--dpparttree', '--partsize', '8', '--retree', '1'], ['--parttree', '--groupsize', '3'], ['--parttree', '--treeout'],
                            ['--parttree', '--reorder'], ['--parttree', '--clustalout']], [('mid', 6), ('small', 10)], n, rng)
    g += cases('treeout', [['--treeout'], ['--nodeout'], ['--distout'], ['--treeout', '--distout'], ['--retree', '0', '--treeout'], ['--retree', '0', '--distout'], ['--retree', '0', '--nodeout'], ['--retree', '1', '--treeout', '--clustalout']],
               basic, n, rng)
    g += cases('treeout-iter', [['--treeout', '--maxiterate', '2'], ['--distout', '--maxiterate', '1'], ['--treeout', '--distout', '--maxiterate', '3'], ['--maxiterate', '2', '--treeout', '--retree', '1'],
                                ['--maxiterate', '1000', '--distout', '--nofft'], ['--treeout', '--maxiterate', '2', '--op', '2'], ['--retree', '0', '--maxiterate', '2', '--treeout'], ['--retree', '0', '--maxiterate', '2', '--distout'],
                                ['--maxiterate', '2', '--treeout', '--anysymbol']], basic, n, rng)
    g += cases('treeout-pair', [['--localpair', '--treeout'], ['--globalpair', '--distout'], ['--genafpair', '--treeout', '--distout'], ['--localpair', '--nodeout'], ['--localpair', '--maxiterate', '2', '--treeout'],
                                ['--globalpair', '--maxiterate', '2', '--distout'], ['--genafpair', '--maxiterate', '3', '--treeout', '--distout'], ['--retree', '0', '--localpair'], ['--retree', '0', '--globalpair', '--treeout'],
                                ['--retree', '0', '--genafpair', '--distout'], ['--retree', '0', '--localpair', '--nodeout'], ['--retree', '0', '--localpair', '--maxiterate', '2', '--treeout', '--distout'],
                                ['--localpair', '--treeout', '--anysymbol'], ['--auto', '--treeout'], ['--auto', '--distout'], ['--auto', '--retree', '0', '--treeout']],
               [('small', 8), ('tiny', 3), ('pair', 3), ('gapped', 2), ('odd', 3)], n, rng)
    g += cases('weighti0', [['--localpair', '--weighti', '0'], ['--globalpair', '--weighti', '0'], ['--genafpair', '--weighti', '0'], ['--localpair', '--weighti', '0', '--maxiterate', '2'], ['--globalpair', '--weighti', '0.0', '--maxiterate', '1000'],
                            ['--genafpair', '--weighti', '0', '--maxiterate', '3'], ['--localpair', '--weighti', '0', '--treeout', '--distout'], ['--localpair', '--weighti', '0', '--maxiterate', '2', '--treeout', '--distout'],
                            ['--weighti', '0'], ['--weighti', '0.5', '--localpair'], ['--weighti', '0', '--auto'], ['--localpair', '--weighti', '0', '--nofft'], ['--localpair', '--weighti', '0', '--fft'], ['--globalpair', '--weighti', '0', '--op', '2', '--ep', '0.1'],
                            ['--localpair', '--weighti', '0', '--bl', '80'], ['--localpair', '--weighti', '0', '--unweight'], ['--localpair', '--weighti', '0', '--anysymbol']],
               [('small', 8), ('tiny', 3), ('pair', 3), ('gapped', 2), ('odd', 3), ('mid', 2)], n, rng)
    opts = [['--minimumweight', '0.01'], ['--minimumweight', '0.5'], ['--randomseed', '5'], ['--bunkatsunashi'], ['--noscore'], ['--adhoc0'], ['--adhoc1'], ['--adhoc2'], ['--scoreout'], ['--outnum'], ['--nomemsave'], ['--alga'],
            ['--simplehillclimbing'], ['--oneiteration'], ['--sparsepickup', '3']]
    bases = [[], ['--maxiterate', '2'], ['--localpair'], ['--localpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '3'], ['--retree', '1', '--maxiterate', '1000']]
    g += cases('passthrough', [b + o for b in bases for o in opts], [('small', 10), ('tiny', 2), ('pair', 3), ('gapped', 2), ('odd', 3)], int(240 * SCALE), rng)
    g += cases('threads', [['--thread', '2'], ['--thread', '4'], ['--thread', '3', '--threadtb', '2'], ['--threadtb', '2'], ['--thread', '2', '--threadit', '0', '--maxiterate', '2'], ['--thread', '2', '--localpair'],
                           ['--thread', '2', '--localpair', '--maxiterate', '2', '--threadit', '0']], basic, n, rng)
    return g


def edge_text(rng):
    """FASTA texts with the awkward things of real files: empty / gap-only sequences, blank lines, CR / CRLF, no final newline, odd titles, bytes above 127, lower case, tabs."""
    a = pick_alpha_for(rng)
    seqs = L.family(rng, a, rng.choice([12, 30, 50]), rng.randint(2, 7))
    feats = rng.sample(['empty', 'gaponly', 'blank', 'crlf', 'cr', 'nofinalnl', 'hi', 'lower', 'tab', 'wrap', 'emptytitle', 'dots', 'star', 'two-gt'], rng.randint(1, 3))
    names = ['s%d' % i for i in range(len(seqs))]
    if 'empty' in feats:
        seqs[rng.randrange(len(seqs))] = ''
    if 'gaponly' in feats:
        seqs[rng.randrange(len(seqs))] = '-' * rng.randint(1, 5)
    if 'hi' in feats:
        i = rng.randrange(len(seqs))
        seqs[i] = seqs[i][:3] + '\u00e9' + seqs[i][3:]
    if 'lower' in feats:
        seqs = [x.lower() for x in seqs]
    if 'dots' in feats:
        i = rng.randrange(len(seqs))
        seqs[i] = seqs[i][:4] + '..' + seqs[i][4:]
    if 'star' in feats:
        i = rng.randrange(len(seqs))
        seqs[i] = seqs[i][:4] + '*' + seqs[i][4:]
    if 'emptytitle' in feats:
        names[rng.randrange(len(names))] = ''
    body = []
    for n, x in zip(names, seqs):
        if 'tab' in feats:
            x = x[:5] + '\t' + x[5:]
        lines = [x[i:i + 20] for i in range(0, len(x), 20)] if 'wrap' in feats else [x]
        body.append('>' + n + '\n' + '\n'.join(lines) + '\n')
        if 'blank' in feats:
            body.append('\n')
    text = ''.join(body)
    if 'nofinalnl' in feats:
        text = text.rstrip('\n')
    if 'crlf' in feats:
        text = text.replace('\n', '\r\n')
    if 'cr' in feats:
        text = text.replace('\n', '\r')
    return text


def pick_alpha_for(rng):
    return L.pick_alpha(rng)


def edge_cases(rng):
    variants = [[], ['--anysymbol'], ['--localpair'], ['--anysymbol', '--localpair', '--maxiterate', '2'], ['--maxiterate', '2'], ['--parttree'], ['--retree', '1', '--anysymbol'], ['--auto'], ['--globalpair', '--maxiterate', '2', '--anysymbol'],
                ['--treeout', '--distout'], ['--clustalout', '--anysymbol'], ['--phylipout'], ['--retree', '0'], ['--genafpair'], ['--localpair', '--weighti', '0']]
    out = []
    for k in range(int(210 * SCALE)):
        out.append(('edge-%04d' % k, 'edge', list(rng.choice(variants)), edge_text(rng), {}))
    return out


def aligned_family(rng, k, a=None):
    """k related sequences aligned by the installed mafft (the existing alignment of --add), and a few more family members (the sequences to add)."""
    a = a or L.pick_alpha(rng)
    base = ''.join(rng.choice(a) for _ in range(rng.choice([20, 40, 80])))
    r = rng.choice([0.05, 0.15, 0.3])
    ind = rng.choice([0.02, 0.05, 0.1])
    fam = [L.mut(rng, base, a, r, ind) for _ in range(k)]
    new = [L.mut(rng, base, a, r, ind) for _ in range(rng.randint(1, 5))]
    rc, out, err, files = L.run_once(['--quiet'], L.fasta(fam))
    return out, new


def add_cases(rng):
    """--add / --addprofile / --keeplength / --mapout (the sequences of new.fa are aligned into the alignment in.fa)."""
    variants = [['--add'], ['--add', '--retree', '1'], ['--add', '--6merpair'], ['--add', '--nofft'], ['--add', '--op', '2', '--ep', '0.1'], ['--add', '--anysymbol'], ['--add', '--reorder'], ['--add', '--clustalout'],
                ['--add', '--treeout', '--distout'], ['--add', '--keeplength'], ['--add', '--mapout'], ['--add', '--compactmapout'], ['--add', '--keeplength', '--anysymbol'], ['--add', '--retree', '3'], ['--add', '--bl', '80'],
                ['--add', '--minimumweight', '0.1'], ['--add', '--thread', '2', '--threadit', '0'], ['--add', '--unweight'], ['--add', '--maxambiguous', '0.5']]
    out = []
    for k in range(int(210 * SCALE)):
        v = list(rng.choice(variants))
        aligned, new = aligned_family(rng, rng.randint(2, 9))
        if '--anysymbol' in v and rng.random() < 0.5:
            new[0] = new[0][:3] + 'X' + new[0][3:]
        fn = v[0]
        text = aligned
        name = 'new.fa'
        out.append(('add-%04d' % k, 'add', v[:1] + [name] + v[1:], text, {name: L.fasta(new, ['n%d' % i for i in range(len(new))])}))
    return out


def pending_cases(rng):
    """Options whose binaries are not ported yet (not embedded)."""
    pend = [['--thread', '2', '--maxiterate', '2', '--adhoc0'], ['--hybridpair'], ['--sp'], ['--seedtable', 'nofile'], ['--seed', 'nofile'], ['--allowshift', '--globalpair'], ['--unalignlevel', '0.5', '--globalpair'], ['--threadtb', '2']]
    opts = [['--skipiterate', '0.5'], ['--skipiterate', '2'], ['--memsave']]
    bases = [[], ['--maxiterate', '2'], ['--localpair', '--maxiterate', '2']]
    return cases('pending', [b + o for b in bases for o in opts], [('small', 10), ('tiny', 2), ('pair', 3), ('odd', 3)], int(120 * SCALE), rng)


def contract_cases(rng):
    """Option / input errors (no alignment): the messages of the script."""
    out = []
    ok = L.fasta(['MKVLAAGIVGLLLA', 'MKVLAGIVGLLA', 'MKVIAAGLVGLLLA'])
    bad_opts = [['--foo'], ['--op'], ['--op', 'x'], ['--ep', 'x'], ['--exp', 'x'], ['--maxiterate', 'abc'], ['--maxiterate', '-1'], ['--retree', 'x'], ['--thread', 'x'], ['--threadit', 'x'], ['--namelength', 'x'],
                ['--groupsize', 'x'], ['--partsize', 'x'], ['--bl', 'x'], ['--kappa', 'x'], ['--weighti', 'x'], ['--opdist', 'x'], ['--mixedlinkage', '2'], ['--mixedlinkage', '-1'], ['--similaritylevel', '2'], ['--similaritylevel', '0.1', '--unalignlevel', '0.1'],
                ['--unalignlevel', '2'], ['--allowshift'], ['--shiftpenalty', '0.5'], ['--skipiterate', '-1'], ['--progress', 'x'], ['--maxambiguous', '2'], ['--algq'], ['--globalgenafpair'],
                ['--topin', 'x'], ['--pdbidlist', 'x'], ['--pdbfilelist', 'x'], ['--merge', 'nofile'], ['--aamatrix', 'nofile'], ['--treein', 'nofile'], ['--add', 'nofile'], ['--addfragments', 'nofile'], ['--keeplength'], ['--mapout'],
                ['--lastpair'], ['--lastmultipair', '--add', 'x'], ['--larapair'], ['--foldalignpair'], ['--scarnapair'], ['--dafspair'], ['--mccaskill'], ['--contrafold'], ['--dafs'], ['--rnaalifold'], ['--blastpair'], ['--fastapair'],
                ['--multipair'], ['--6mermultipair'], ['--parttree', '--maxiterate', '2'], ['--parttree', '--seed', 'x'], ['--text', '--anysymbol'], ['--retree', '0', '--nodeout', '--maxiterate', '2'], ['--help'], ['--info'], ['--version'],
                ['--parttree', '--merge', 'nofile'], ['--bestfirst', '--threadit', '0'], ['--dash'], ['--anchors', 'nofile']]
    for i, o in enumerate(bad_opts):
        out.append(('contract-%03d' % i, 'contract', o, ok, {}))
    # input problems
    out.append(('contract-empty', 'contract', [], '', {}))
    out.append(('contract-blank-led', 'contract', [], ' >s0\nACGT\n>s1\nACGA\n>s2\nACGG\n', {}))
    out.append(('contract-one', 'contract', [], '>a\nACGT\n', {}))
    out.append(('contract-noseq', 'contract', [], '>a\n>b\n>c\n', {}))
    out.append(('contract-emptyseq', 'contract', [], '>a\nACGTACGT\n>b\n\n>c\nACGTAC\n', {}))
    out.append(('contract-nofasta', 'contract', [], 'ACGT\nACGA\n', {}))
    out.append(('contract-crlf', 'contract', [], '>a\r\nACGTACGT\r\n>b\r\nACGTTCGT\r\n>c\r\nACGAACGT\r\n', {}))
    out.append(('contract-cr', 'contract', ['--anysymbol'], '>a\rACGTACGT\r>b\rACGTTCGT\r>c\rACGAACGT\r', {}))
    out.append(('contract-gt-inside', 'contract', ['--anysymbol'], '>a\nACGT>ACGT\n>b\nACGTTCGT\n>c\nACGAACGT\n', {}))
    out.append(('contract-eq-inside', 'contract', ['--anysymbol'], '>a\nACGT=ACGT\n>b\nACGTTCGT\n>c\nACGAACGT\n', {}))
    out.append(('contract-two', 'contract', [], '>a\nACGTACGTAA\n>b\nACGTTCGTAA\n', {}))
    out.append(('contract-hi', 'contract', ['--anysymbol'], '>a\nMKVLAAG\xe9IVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n', {}))
    out.append(('contract-hi-plain', 'contract', [], '>a\nMKVLAAG\xe9IVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n', {}))
    out.append(('contract-digits', 'contract', ['--anysymbol'], '>a\nMK12VLAAGIVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n', {}))
    out.append(('contract-digits-plain', 'contract', [], '>a\nMK12VLAAGIVGLLLA\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n', {}))
    out.append(('contract-lower', 'contract', ['--clustalout'], '>a\nmkvlaagivgllla\n>b\nmkvlagivglla\n>c\nmkviaaglvgllla\n', {}))
    out.append(('contract-dots', 'contract', [], '>a\nMKVL..AAGIVG\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n', {}))
    out.append(('contract-star', 'contract', [], '>a\nMKVL*AAGIVG\n>b\nMKVLAGIVGLLA\n>c\nMKVIAAGLVGLLLA\n', {}))
    out.append(('contract-out', 'contract', ['--out', 'result.txt'], ok, {}))
    out.append(('contract-out-quiet', 'contract', ['--out', 'result.txt', '--quiet'], ok, {}))
    return out


def threaded(args):
    for i, a in enumerate(args):
        if a in ('--thread', '--threadtb') and i + 1 < len(args) and args[i + 1] not in ('0', '1'):
            return True
        if a == '--threadit' and i + 1 < len(args) and args[i + 1] not in ('0',):
            return True
    return False


def run_case(case):
    cid, group, args, text, extra = case
    (a, det) = L.run_checked(args, text, extra)
    rc, out, err, files = a
    nostderr = threaded(args)
    if not det and nostderr:
        # the progress text of the threads is racy: the rest must still be reproducible
        b = L.run_once(args, text, extra)
        det = (a[0], a[1], a[3]) == (b[0], b[1], b[3])
    if not det:
        return None, cid
    # the temporary directory of the run is named in some messages: refuse to record those
    return dict(id=cid, group=group, args=args, name=L.INPUT_NAME, input=text, files=extra, rc=rc, out=out, err=err, outfiles=files, nostderr=nostderr), None


def main():
    rng = random.Random(7050515)
    want = set(sys.argv[1:])
    allcases = groups(rng) + edge_cases(rng) + (add_cases(rng) if (not want or 'add' in want) else []) + contract_cases(rng) + (pending_cases(rng) if 'pending' in want else [])
    if want:
        allcases = [c for c in allcases if c[1] in want]
    n = len(allcases)
    dropped = []
    done = 0
    tsv = os.environ.get('MC_TSV')
    ft = gzip.open(tsv, 'wt', encoding='utf-8', newline='\n') if tsv else None
    full = 0
    with ThreadPoolExecutor(THREADS) as ex, open(OUT, 'w') as fo:
        for res, bad in ex.map(run_case, allcases):
            done += 1
            if res is None:
                dropped.append(bad)
                continue
            fo.write(json.dumps(res) + '\n')
            if ft and res['group'] != 'pending':
                sha = '-' if res['nostderr'] else hashlib.sha256(res['err'].encode('utf-8')).hexdigest()
                keep = (not res['nostderr']) and (full % 25 == 0 or len(res['err']) < 400)
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
