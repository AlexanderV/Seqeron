#!/usr/bin/env python3
"""Builds the oracle of the options of the script mafft that WP74 of batch B09 ported (F340 ff.): one group per option family, every case is one run of the installed script `mafft` (Debian 7.505-1), repeated twice
(a case whose runs differ is dropped: the program is non-deterministic there); the oracle holds return code, standard output, standard error and every file the run wrote.

Groups: merge (--merge table: subalignments that the programs keep, the constrained UPGMA tree, disttbfast), merge-pair (--merge with --localpair / --globalpair / --genafpair: tbfast), merge-iter (--merge with --maxiterate: dvtditr),
merge-err (tables and groups the programs refuse), merge-tree (--merge with --treein); text (--text: an alphabet of 256 characters, default progressive method), text-pair / text-iter (with a pair library / --maxiterate),
text-matrix (--textmatrix: the entries of the file _aamtx), text-err (stops: '=' '<' '>' in a text, a malformed matrix file), text-format (clustal / phylip output, --reorder, --treeout ...);
seed (--seed file: multi2hat3s, disttbfast -E 1, tbfast -J, dvtditr -l with the structure weights), seed-multi (two or three seed files), seedtable (--seedtable: the hat3.seed file of the user), seed-err (stops),
shift (--shiftpenalty x with --globalpair: trywarp of A__align, G__align11, partA__align for x below 10, the stops of the script for the other methods), seed-large (--seed with --large / --memsavetree / --stepadd and a pair library: nodepair and the seeds of tbfast -U), seed-opt (--seed with the other options of the script: --auto, --retree 0, --large, --treein, --merge, --add, ...).

It reuses the infrastructure of WP70 - WP73 (generate_mafft_command_oracle.py: sets, option runner; mafft_command_oracle_lib.py: runner, writer).

Usage: python3 scripts/alignment/generate_mafft_wp74_oracle.py [group ...]   (no group = all)   env MC_SCALE (default 1), MC_OUT (JSONL path), MC_TSV (TSV, gzip when it ends in .gz), MC_THREADS
"""
import os, random, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mafft_command_oracle_lib as L
import generate_mafft_command_oracle as C
import generate_mafft_compacttree_oracle as T

SCALE = float(os.environ.get('MC_SCALE', '1'))
OUT = os.environ.get('MC_OUT', '/tmp/mafft_wp74_oracle.jsonl')
THREADS = int(os.environ.get('MC_THREADS', '8'))
TABLE = 'table.txt'
TREE = T.TREE


def families(rng, alpha, nfam, size_range, length=None):
    """nfam families (a random base each, so the families are distant), each of 2.. members; returns the list of families."""
    fams = []
    for _ in range(nfam):
        base = ''.join(rng.choice(alpha) for _ in range(length or rng.choice([20, 40, 70])))
        r = rng.choice([0.03, 0.08, 0.15])
        ind = rng.choice([0.02, 0.05])
        k = rng.randint(*size_range)
        fams.append([L.mut(rng, base, alpha, r, ind) for _ in range(k)])
    return fams


def aligned_rows(seqs):
    """The rows of the alignment of `seqs` made by the installed mafft."""
    rc, out, err, files = L.run_once(['--quiet'], L.fasta(seqs, ['t%d' % i for i in range(len(seqs))]))
    rows = []
    for rec in out.split('>')[1:]:
        rows.append(''.join(rec.split('\n')[1:]))
    return rows


def make_case(rng, alpha=None, nfam=None, merge_fams=None, size_range=(2, 4), singles=None):
    """The input file (aligned families that are subalignments + unaligned singletons, shuffled), the names, and the 1-based groups."""
    alpha = alpha or L.pick_alpha(rng)
    nfam = nfam or rng.choice([1, 2, 2, 3])
    fams = families(rng, alpha, nfam, size_range)
    entries = []  # (seq, group index or None)
    groups = []
    for f in fams:
        if rng.random() < 0.85 and len(f) >= 2:
            rows = aligned_rows(f)
            gi = len(groups)
            groups.append([])
            for r in rows:
                entries.append((r, gi))
        else:
            for s in f:
                entries.append((s, None))
    for _ in range(rng.choice([0, 1, 1, 2]) if singles is None else singles):
        entries.append((''.join(rng.choice(alpha) for _ in range(rng.choice([15, 40, 70]))), None))
    if len(entries) < 3:
        entries.append((''.join(rng.choice(alpha) for _ in range(30)), None))
    order = list(range(len(entries)))
    if rng.random() < 0.7:
        rng.shuffle(order)
    seqs = []
    for pos, ei in enumerate(order):
        seq, gi = entries[ei]
        seqs.append(seq)
        if gi is not None:
            groups[gi].append(pos + 1)
    names = C.names_for(rng, len(seqs))
    return seqs, names, groups


def table_text(rng, groups, negative=0.0, junk=True):
    lines = []
    if junk and rng.random() < 0.4:
        lines.append('# the subalignments')
    for g in groups:
        nums = [str(x) for x in g]
        if negative and rng.random() < negative:
            k = rng.randrange(len(nums))
            nums[k] = '-' + nums[k]
        sep = rng.choice([' ', ' ', ' ', '  ', '\t'])
        line = sep.join(nums)
        if junk and rng.random() < 0.2:
            line += ' # comment 3 4'
        if junk and rng.random() < 0.1:
            line = '  ' + line
        lines.append(line)
        if junk and rng.random() < 0.15:
            lines.append('')
    if junk and rng.random() < 0.2:
        lines.append('# end')
    text = '\n'.join(lines) + ('\n' if rng.random() < 0.8 else '')
    if junk and rng.random() < 0.1:
        text = text.replace('\n', '\r\n')
    return text


def consistent_tree(rng, n, groups):
    """A guide tree file (newick2mafft.rb format) in which every group is a monophyletic cluster."""
    active = list(range(n))
    lines = []

    def join(a, b):
        a, b = sorted((a, b))
        active.remove(b)
        l0, l1 = rng.random() * 0.5, rng.random() * 0.5
        lines.append('%d %d %s %s' % (a + 1, b + 1, T.fmt_len(rng, l0), T.fmt_len(rng, l1)))

    for g in groups:
        members = sorted(abs(x) - 1 for x in g)
        reps = list(members)
        while len(reps) > 1:
            a, b = rng.sample(reps, 2)
            lo, hi = min(a, b), max(a, b)
            join(lo, hi)
            reps.remove(hi)
    while len(active) > 1:
        a, b = rng.sample(active, 2)
        join(a, b)
    return '\n'.join(lines) + '\n'


BASES = [[], [], [], ['--retree', '1'], ['--retree', '2'], ['--retree', '3'], ['--nofft'], ['--retree', '1', '--nofft'], ['--op', '2', '--ep', '0.1'], ['--6merpair'], ['--treeout'], ['--distout'], ['--treeout', '--distout'],
         ['--clustalout'], ['--reorder'], ['--anysymbol'], ['--unweight'], ['--minimumweight', '0.1'], ['--bl', '80'], ['--thread', '2'], ['--mixedlinkage', '0.5'], ['--averagelinkage'], ['--minimumlinkage'], ['--retree', '0'],
         ['--retree', '0', '--treeout'], ['--outnum'], ['--clustalout', '--reorder'], ['--amino'], ['--nuc'], ['--retree', '0', '--distout']]


def build(rng, want):
    out = []
    n = int(330 * SCALE)
    if not want or 'merge' in want:
        for k in range(n):
            seqs, names, groups = make_case(rng)
            extra = {TABLE: table_text(rng, groups, negative=0.15)}
            v = ['--merge', TABLE] + list(rng.choice(BASES))
            if '--anysymbol' in v:
                i = rng.randrange(len(seqs))
            out.append(('merge-%04d' % k, 'merge', v, L.fasta(seqs, names, rng.choice([None, None, 30])), extra))
    if not want or 'merge-pair' in want:
        pairs = [['--localpair'], ['--globalpair'], ['--genafpair'], ['--oldgenafpair'], ['--localpair', '--weighti', '0'], ['--globalpair', '--weighti', '0'], ['--localpair', '--weighti', '1.5'], ['--genafpair', '--retree', '1'],
                 ['--localpair', '--treeout'], ['--globalpair', '--distout'], ['--localpair', '--op', '2', '--ep', '0.1'], ['--localpair', '--anysymbol'], ['--localpair', '--reorder'], ['--localpair', '--clustalout']]
        for k in range(n):
            seqs, names, groups = make_case(rng, size_range=(2, 3))
            extra = {TABLE: table_text(rng, groups, negative=0.1)}
            v = ['--merge', TABLE] + list(rng.choice(pairs))
            out.append(('merge-pair-%04d' % k, 'merge-pair', v, L.fasta(seqs, names), extra))
    if not want or 'merge-iter' in want:
        its = [['--maxiterate', '1'], ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--maxiterate', '2', '--retree', '1'], ['--maxiterate', '2', '--nofft'], ['--maxiterate', '3', '--treeout'], ['--maxiterate', '2', '--distout'],
               ['--maxiterate', '2', '--thread', '3', '--threadit', '0'], ['--localpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '2'], ['--genafpair', '--maxiterate', '2'], ['--maxiterate', '2', '--skipiterate', '0.3'],
               ['--maxiterate', '2', '--anysymbol'], ['--maxiterate', '2', '--bunkatsunashi'], ['--maxiterate', '2', '--clustalout'], ['--maxiterate', '2', '--reorder'], ['--maxiterate', '2', '--op', '2'], ['--maxiterate', '2', '--core']]
        for k in range(n):
            seqs, names, groups = make_case(rng)
            extra = {TABLE: table_text(rng, groups, negative=0.1)}
            v = ['--merge', TABLE] + list(rng.choice(its))
            out.append(('merge-iter-%04d' % k, 'merge-iter', v, L.fasta(seqs, names), extra))
    if not want or 'merge-err' in want:
        for k in range(n):
            seqs, names, groups = make_case(rng, singles=rng.choice([1, 2, 3]))
            nseq = len(seqs)
            kind = rng.choice(['unaligned', 'nonmono', 'dup', 'range', 'format', 'empty', 'comment', 'single', 'overlap', 'neg', 'zero', 'far', 'nolead', 'two-lines', 'offbyone', 'allseq', 'pairs'])
            g = [list(x) for x in groups]
            if kind == 'unaligned':
                pick = rng.sample(range(1, nseq + 1), min(nseq, rng.randint(2, 3)))
                g.append(pick)
            elif kind == 'nonmono':
                pick = rng.sample(range(1, nseq + 1), min(nseq, rng.randint(2, 3)))
                g = [pick]
            elif kind == 'dup':
                g = g + [[rng.randint(1, nseq), rng.randint(1, nseq)]]
                if g[-1][0] != g[-1][1]:
                    g[-1][1] = g[-1][0]
            elif kind == 'range':
                g = g + [[rng.randint(nseq + 1, nseq + 2), 1]]
            elif kind == 'offbyone':
                g = [[x + 1 for x in gg] for gg in g]
            elif kind == 'neg':
                g = [[-x if rng.random() < 0.5 else x for x in gg] for gg in g]
            elif kind == 'single':
                g = g + [[rng.randint(1, nseq)]]
            elif kind == 'overlap':
                if g:
                    g = g + [[g[0][0]]]
            elif kind == 'allseq':
                g = [list(range(1, nseq + 1))]
            elif kind == 'pairs':
                g = [[i, i + 1] for i in range(1, nseq, 2)]
            text = table_text(rng, g, junk=True)
            if kind == 'format':
                text = text + rng.choice(['1 x 3\n', 'a b\n', '1 2x\n', '5 q\n', '2 3 +4\n', '1 2 3 4abc\n'])
            elif kind == 'empty':
                text = ''
            elif kind == 'comment':
                text = '# nothing\n# here\n'
            elif kind == 'zero':
                text = text + '0 1\n'
            elif kind == 'far':
                text = text + '%d 1\n' % rng.choice([nseq + 5, nseq + 20, 1000])
            elif kind == 'nolead':
                text = 'x\n' + text
            elif kind == 'two-lines':
                text = text + text
            base = rng.choice([[], [], ['--maxiterate', '2'], ['--localpair'], ['--retree', '1'], ['--clustalout'], ['--quiet'], ['--treein', TREE]])
            extra = {TABLE: text}
            if '--treein' in base:
                extra[TREE] = T.tree_file(rng, nseq)
            out.append(('merge-err-%04d' % k, 'merge-err', ['--merge', TABLE] + base, L.fasta(seqs, names), extra))
    if not want or 'merge-tree' in want:
        for k in range(n):
            seqs, names, groups = make_case(rng, size_range=(2, 3))
            nseq = len(seqs)
            extra = {TABLE: table_text(rng, groups), TREE: consistent_tree(rng, nseq, groups) if rng.random() < 0.7 else T.tree_file(rng, nseq)}
            base = rng.choice([[], ['--retree', '1'], ['--maxiterate', '2'], ['--localpair'], ['--nofft'], ['--treeout'], ['--clustalout'], ['--maxiterate', '2', '--localpair']])
            out.append(('merge-tree-%04d' % k, 'merge-tree', ['--merge', TABLE, '--treein', TREE] + base, L.fasta(seqs, names), extra))
    build_text(rng, want, int(210 * SCALE), out)
    build_seed(rng, want, int(210 * SCALE), out)
    return out


def table_overruns(text, nseq):
    """True when the table names a sequence number that lands behind the padding of the program's block tab01 (calloc( nseq, sizeof( int ) )): the program then reads heap that it does not own (valgrind: invalid read / write
    of size 4, 'bytes after a block'), what it finds there depends on the heap layout: the port refuses these tables."""
    usable = max(32, (nseq * 4 + 8 + 15) // 16 * 16) - 8
    slots = usable // 4
    for line in text.replace('\r', '\n').split('\n'):
        if not line or line[0] == '#':
            continue
        for tok in line.replace('\t', ' ').split(' '):
            if tok.startswith('#'):
                break
            digits = ''.join(ch for ch in tok.lstrip('-+') if ch.isdigit())
            if digits and int(digits) - 1 >= slots:
                return True
    return False



TEXT_ALPHABETS = [
    'abcdefghijklmnopqrstuvwxyz', 'ABCDEFGHIJKLMNOPQRSTUVWXYZ', 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ', '0123456789', 'abcxyz0123', 'ACGT', 'ACDEFGHIKLMNPQRSTVWY', 'aeiou',
    '!"#$%&\'()*+,./:;?@[\\]^_`{|}~', 'abc012!?', 'ab',
]
TEXT_NASTY = ['\u00e9', '\u00fc\u00df', '\t', '\x01', '\x7f', '\u00ff', '*', '.']


def text_set(rng, kind=None):
    """Families of text sequences over an alphabet of the 256-letter alphabet (mixed case, digits, punctuation, now and then a byte above 127 or a control character)."""
    alpha = rng.choice(TEXT_ALPHABETS) if kind is None else kind
    if rng.random() < 0.25:
        alpha = alpha + rng.choice(TEXT_ALPHABETS)
    seqs = L.family(rng, alpha, rng.choice([10, 30, 60, 100]), rng.randint(3, 14))
    out = []
    for s in seqs:
        s = list(s)
        if rng.random() < 0.3:
            for _ in range(rng.randrange(1, 3)):
                s.insert(rng.randrange(len(s) + 1), rng.choice(TEXT_NASTY))
        if rng.random() < 0.2:
            s.insert(rng.randrange(len(s) + 1), ' ')
        if rng.random() < 0.15:
            s.insert(rng.randrange(len(s) + 1), '-')
        out.append(''.join(s))
    return out


def text_input(rng, wrap=None):
    seqs = text_set(rng)
    return seqs, L.fasta(seqs, C.names_for(rng, len(seqs)), wrap)


def matrix_file(rng, seqs):
    """A file for --textmatrix: '0xAA 0xBB score' lines for letters of the data (match scores and a few mismatches), comments."""
    letters = sorted(set(ch for s in seqs for ch in s if ord(ch) < 256 and ch not in ' -=<>\r\n\x00'))
    lines = []
    if rng.random() < 0.4:
        lines.append('# the scores of this text')
    for _ in range(rng.randint(1, 8)):
        if not letters:
            break
        a = rng.choice(letters)
        b = a if rng.random() < 0.5 else rng.choice(letters)
        score = rng.choice(['5', '2.5', '-3', '1', '10', '0', '0.5', '-1.25', '7e0'])
        fmt = rng.choice(['0x%x 0x%x %s', '0x%X 0x%X %s', '%x %x %s', '0x%02x  0x%02x\t%s'])
        line = fmt % (ord(a), ord(b), score)
        if rng.random() < 0.3:
            line += ' # ' + rng.choice(['match', 'x y z', '0x41 0x42 3'])
        lines.append(line)
        if rng.random() < 0.1:
            lines.append('')
    text = '\n'.join(lines) + ('\n' if rng.random() < 0.85 else '')
    if rng.random() < 0.08:
        text = text.replace('\n', '\r\n')
    return text


TEXT_BASES = [[], [], [], ['--retree', '1'], ['--retree', '2'], ['--retree', '3'], ['--retree', '0'], ['--op', '2', '--ep', '0.1'], ['--treeout'], ['--distout'], ['--reorder'], ['--unweight'], ['--minimumweight', '0.1'],
              ['--thread', '2'], ['--mixedlinkage', '0.5'], ['--averagelinkage'], ['--minimumlinkage'], ['--outnum'], ['--6merpair'], ['--nofft'], ['--op', '0.5'], ['--exp', '0.1'], ['--leavegappyregion'], ['--legacygappenalty'],
              ['--termgappenalty'], ['--quiet'], ['--fmodel'], ['--sparsepickup', '3'], ['--oneiteration'], ['--scoreout'], ['--bunkatsunashi']]


def build_text(rng, want, n, out):
    if not want or 'text' in want:
        for k in range(n):
            seqs, text = text_input(rng, rng.choice([None, None, 30]))
            out.append(('text-%04d' % k, 'text', ['--text'] + list(rng.choice(TEXT_BASES)), text, {}))
    if not want or 'text-format' in want:
        fmts = [['--clustalout'], ['--phylipout'], ['--clustalout', '--reorder'], ['--phylipout', '--reorder'], ['--clustalout', '--namelength', '8'], ['--phylipout', '--namelength', '20'], ['--clustalout', '--retree', '1'],
                ['--clustalout', '--maxiterate', '2'], ['--clustalout', '--localpair'], ['--phylipout', '--maxiterate', '2'], ['--clustalout', '--treeout'], ['--out', 'result.txt'], ['--out', 'result.txt', '--clustalout']]
        for k in range(n):
            seqs, text = text_input(rng, rng.choice([None, 30]))
            out.append(('text-format-%04d' % k, 'text-format', ['--text'] + list(rng.choice(fmts)), text, {}))
    if not want or 'text-pair' in want:
        pairs = [['--localpair'], ['--globalpair'], ['--genafpair'], ['--oldgenafpair'], ['--localpair', '--weighti', '0'], ['--globalpair', '--weighti', '0'], ['--localpair', '--weighti', '1.5'], ['--genafpair', '--retree', '1'],
                 ['--localpair', '--treeout'], ['--globalpair', '--distout'], ['--localpair', '--op', '2', '--ep', '0.1'], ['--localpair', '--reorder'], ['--localpair', '--clustalout'], ['--localpair', '--lop', '-3', '--lep', '0.2']]
        for k in range(n):
            seqs, text = text_input(rng)
            seqs = seqs[:rng.randint(2, 7)]
            out.append(('text-pair-%04d' % k, 'text-pair', ['--text'] + list(rng.choice(pairs)), L.fasta(seqs, ['s%d' % i for i in range(len(seqs))]), {}))
    if not want or 'text-iter' in want:
        its = [['--maxiterate', '1'], ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--maxiterate', '2', '--retree', '1'], ['--maxiterate', '3', '--treeout'], ['--maxiterate', '2', '--distout'], ['--maxiterate', '2', '--thread', '3', '--threadit', '0'],
               ['--localpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '2'], ['--genafpair', '--maxiterate', '2'], ['--maxiterate', '2', '--skipiterate', '0.3'], ['--maxiterate', '2', '--bunkatsunashi'],
               ['--maxiterate', '2', '--core'], ['--maxiterate', '2', '--op', '2'], ['--auto']]
        for k in range(n):
            seqs, text = text_input(rng)
            seqs = seqs[:rng.randint(3, 9)]
            out.append(('text-iter-%04d' % k, 'text-iter', ['--text'] + list(rng.choice(its)), L.fasta(seqs, ['s%d' % i for i in range(len(seqs))]), {}))
    if not want or 'text-matrix' in want:
        extra_opts = [[], [], ['--retree', '1'], ['--maxiterate', '2'], ['--localpair'], ['--clustalout'], ['--op', '2'], ['--treeout'], ['--globalpair', '--maxiterate', '2'], ['--reorder'], ['--genafpair']]
        for k in range(n):
            seqs, text = text_input(rng)
            mtx = matrix_file(rng, seqs)
            out.append(('text-matrix-%04d' % k, 'text-matrix', ['--textmatrix', 'mtx.txt'] + list(rng.choice(extra_opts)), text, {'mtx.txt': mtx}))
    if not want or 'text-err' in want:
        for k in range(n):
            seqs, text = text_input(rng)
            kind = rng.choice(['eq', 'lt', 'gt', 'anysymbol', 'nomatrix', 'badline', 'twovalues', 'badchar', 'big', 'dash', 'space', 'noletters', 'negmatch', 'zeromatrix', 'preservecase', 'hexonly', 'empty', 'crlf'])
            extra = {}
            args = ['--text']
            if kind in ('eq', 'lt', 'gt'):
                ch = {'eq': '=', 'lt': '<', 'gt': '>'}[kind]
                i = rng.randrange(len(seqs))
                seqs[i] = seqs[i][:3] + ch + seqs[i][3:]
                text = L.fasta(seqs, ['s%d' % j for j in range(len(seqs))])
                args += rng.choice([[], ['--localpair'], ['--maxiterate', '2'], ['--add', 'x']]) if kind != 'gt' else []
            elif kind in ('anysymbol', 'preservecase'):
                args += ['--' + kind]
            elif kind == 'nomatrix':
                args = ['--textmatrix', 'missing.txt']
            else:
                lines = {'badline': '0x41 0x42\n', 'twovalues': '0x41 0x42 3 4\nx y z\n', 'badchar': '0x00 0x41 3\n', 'big': '0x141 0x41 3\n', 'dash': '0x2d 0x41 3\n', 'space': '0x20 0x41 3\n',
                         'noletters': 'zz zz zz\n', 'negmatch': '0x41 0x41 -50\n0x61 0x61 -50\n0x62 0x62 -50\n', 'zeromatrix': '0x41 0x41 0\n', 'hexonly': '41 42 3.5\n', 'empty': '', 'crlf': '0x41 0x41 4\r\n'}[kind]
                extra['mtx.txt'] = lines
                args = ['--textmatrix', 'mtx.txt']
            out.append(('text-err-%04d' % k, 'text-err', args, text, extra))
    return out



SEED_BASES = [[], [], [], ['--retree', '1'], ['--retree', '3'], ['--maxiterate', '3'], ['--maxiterate', '1000'], ['--nofft'], ['--op', '2', '--ep', '0.1'], ['--treeout'], ['--distout'], ['--clustalout'], ['--reorder'], ['--quiet'], ['--6merpair'],
              ['--unweight'], ['--thread', '2', '--threadit', '0'], ['--weighti', '1'], ['--weighti', '0'], ['--minimumweight', '0.1'], ['--bl', '80'], ['--mixedlinkage', '0.5'], ['--outnum'], ['--maxiterate', '2', '--bunkatsunashi'],
              ['--amino'], ['--nuc'], ['--core'], ['--skipiterate', '0.3'], ['--localpair'], ['--globalpair'], ['--genafpair'], ['--localpair', '--maxiterate', '2'], ['--globalpair', '--maxiterate', '1000'],
              ['--genafpair', '--maxiterate', '3'], ['--localpair', '--weighti', '0'], ['--localpair', '--treeout', '--distout'], ['--localpair', '--clustalout']]


SEED_OPT_STOP = [['--large'], ['--pileup'], ['--randomchain'], ['--randomchain', '--randomseed', '5'], ['--memsavetree'], ['--stepadd'], ['--youngestlinkage'], ['--dpparttree'], ['--parttree']]
SEED_OPT = [['--auto'], ['--retree', '0'], ['--retree', '2'], ['--memsave'], ['--nomemsave'],
            ['--phylipout'], ['--inputorder'], ['--namelength', '10', '--clustalout'], ['--oldgenafpair'], ['--lop', '-1.0', '--lep', '0.2', '--localpair'], ['--gop', '1.0', '--gep', '0.1'], ['--leavegappyregion'], ['--termgappenalty', '0.5'],
            ['--fft'], ['--10merpair'], ['--kimura', '200'], ['--jtt', '100'], ['--tm', '100'], ['--fmodel'], ['--nwildcard'], ['--text'], ['--scoreout'], ['--noscore'], ['--preservecase'], ['--lexp', '-0.2', '--localpair'],
            ['--averagelinkage'], ['--minimumlinkage'], ['--kappa', '3', '--nuc'], ['--opdist', '1.0'], ['--exp', '0.2'], ['--sparsepickup', '3'], ['--smoothing'], ['--thread', '2', '--threadtb', '2', '--threadit', '0'], ['--outnum', '--reorder'],
            ['--maxiterate', '2', '--skipiterate', '0.5'], ['--core', '--coreext'], ['--corethr', '0.5', '--corewin', '3', '--core'], ['--anysymbol', '--treeout'], ['--amino', '--aamatrix', 'x.mtx']]


def seed_case(rng, nseedfiles=1, tag=''):
    """Seed alignments (families aligned by the installed mafft) and the sequences to align (members of the same families, unaligned, and a few unrelated ones)."""
    alpha = L.pick_alpha(rng)
    fams = families(rng, alpha, nseedfiles + rng.choice([0, 1]), (3, 5), None)
    seeds = []
    for f in fams[:nseedfiles]:
        rows = aligned_rows(f)
        names = ['q%d' % i for i in range(len(rows))]
        seeds.append(L.fasta(rows, names))
    inputs = []
    for f in fams:
        inputs += [L.mut(rng, rng.choice(f), alpha, 0.08, 0.03) for _ in range(rng.randint(1, 3))]
    for _ in range(rng.choice([0, 1, 2])):
        inputs.append(''.join(rng.choice(alpha) for _ in range(rng.choice([20, 50]))))
    if len(inputs) < 2:
        inputs.append(''.join(rng.choice(alpha) for _ in range(30)))
    return seeds, L.fasta(inputs, C.names_for(rng, len(inputs)), rng.choice([None, 40]))


def odd_rows(rng, fasta_text):
    """Puts a few unusual letters (U, X, #, digits, a star) into the sequences of a FASTA text (for --anysymbol)."""
    recs = fasta_text.split('>')[1:]
    out = []
    for r in recs:
        name, _, body = r.partition('\n')
        body = list(body.replace('\n', ''))
        for _ in range(rng.randrange(0, 3)):
            if body:
                body[rng.randrange(len(body))] = rng.choice('UXBZ#*7')
        out.append('>' + name + '\n' + ''.join(body) + '\n')
    return ''.join(out)


def build_seed(rng, want, n, out):
    if not want or 'seed' in want:
        for k in range(n):
            seeds, text = seed_case(rng)
            base = list(rng.choice(SEED_BASES))
            files = {'seed1.fa': seeds[0]}
            if rng.random() < 0.2:
                base = ['--anysymbol'] + rng.choice([[], ['--maxiterate', '3'], ['--clustalout'], ['--reorder'], ['--retree', '3'], ['--treeout']])
                files['seed1.fa'] = odd_rows(rng, seeds[0])
                text = odd_rows(rng, text)
            out.append(('seed-%04d' % k, 'seed', ['--seed', 'seed1.fa'] + base, text, files))
    if not want or 'seed-multi' in want:
        for k in range(n):
            ns = rng.choice([2, 2, 3])
            seeds, text = seed_case(rng, ns)
            args = []
            files = {}
            base = list(rng.choice(SEED_BASES))
            anysym = rng.random() < 0.2
            for i, sd in enumerate(seeds):
                args += ['--seed', 'seed%d.fa' % (i + 1)]
                files['seed%d.fa' % (i + 1)] = odd_rows(rng, sd) if anysym else sd
            if anysym:
                base = ['--anysymbol'] + rng.choice([[], ['--maxiterate', '3'], ['--clustalout']])
                text = odd_rows(rng, text)
            out.append(('seed-multi-%04d' % k, 'seed-multi', args + base, text, files))
    if not want or 'seedtable' in want:
        for k in range(n):
            seeds, text = seed_case(rng)
            # a table of the blocks of the seed alignment over the sequences of the input: pairs (i, j) of the input with blocks of the seed
            nin = text.count('>')
            lines = []
            for _ in range(rng.randint(1, 3 * nin)):
                i = rng.randrange(nin - 1)
                j = rng.randrange(i + 1, nin)
                st1 = rng.randrange(0, 15)
                ln = rng.randint(3, 12)
                lines.append('%d %d %d %.3f %d %d %d %d %s' % (i, j, ln, rng.choice([2400.0, 12000.5, 14632.917, 3000.25]), st1, st1 + ln - 1, st1 + rng.randrange(0, 3), st1 + ln - 1, rng.choice('kk h'.replace(' ', ''))))
            tbl = '\n'.join(lines) + '\n'
            out.append(('seedtable-%04d' % k, 'seedtable', ['--seedtable', 'tab.txt'] + list(rng.choice(SEED_BASES)), text, {'tab.txt': tbl}))
    if not want or 'seed-err' in want:
        for k in range(n):
            seeds, text = seed_case(rng)
            kind = rng.choice(['single', 'unaligned', 'both', 'badtable', 'twofiles-one', 'noseq', 'unknownletter', 'anysymbol', 'pair', 'crlf'])
            args = ['--seed', 'seed1.fa']
            files = {'seed1.fa': seeds[0]}
            if kind == 'single':
                files['seed1.fa'] = seeds[0].split('>')[1].join(['>', ''])
            elif kind == 'unaligned':
                rows = [r for r in seeds[0].split('>')[1:]]
                files['seed1.fa'] = '>' + rows[0] + '>' + rows[1].rstrip('\n') + 'ACGT\n' + ''.join('>' + r for r in rows[2:])
            elif kind == 'both':
                args = ['--seed', 'seed1.fa', '--seedtable', 'tab.txt']
                files['tab.txt'] = '0 1 5 100.0 0 4 0 4 k\n'
            elif kind == 'badtable':
                args = ['--seedtable', 'tab.txt']
                files = {'tab.txt': rng.choice(['1 0 5 100.0 0 4 0 4 k\n', '5 1 5 100.0 0 4 0 4 k\n', '2 2 5 100.0 0 4 0 4 k\n', '0 1 5 100.0 0 4 0 4 k\n0 1 5 100.0 6 9 6 9 k\n'])}
            elif kind == 'twofiles-one':
                args = ['--seed', 'seed1.fa', '--seed', 'seed2.fa']
                files['seed2.fa'] = seeds[0]
            elif kind == 'noseq':
                files['seed1.fa'] = '>a\n>b\n'
            elif kind == 'unknownletter':
                files['seed1.fa'] = seeds[0].replace('\n', 'U\n', 1) if '>' in seeds[0] else seeds[0]
            elif kind == 'anysymbol':
                args += ['--anysymbol']
            elif kind == 'pair':
                args += ['--localpair']
            elif kind == 'crlf':
                files['seed1.fa'] = seeds[0].replace('\n', '\r\n')
            out.append(('seed-err-%04d' % k, 'seed-err', args + rng.choice([[], ['--maxiterate', '3'], ['--quiet']]), text, files))
    if not want or 'seed-opt' in want:
        for k in range(int(n * 1.3)):
            seeds, text = seed_case(rng)
            nin = text.count('>')
            nseed = seeds[0].count('>')
            files = {'seed1.fa': seeds[0]}
            kind = rng.choice(['plain'] * 5 + ['treein'] * 2 + ['merge'] * 2 + ['add'])
            args = ['--seed', 'seed1.fa']
            if kind == 'treein':
                files[TREE] = consistent_tree(rng, nin + rng.choice([0, nseed, nseed]), [])
                args += ['--treein', TREE]
            elif kind == 'merge':
                k2 = nin + nseed
                members = list(range(1, k2 + 1))
                rng.shuffle(members)
                files[TABLE] = table_text(rng, [sorted(members[:rng.randint(2, max(2, k2 // 2))])], negative=0.1)
                args += ['--merge', TABLE]
            elif kind == 'add':
                files['add.fa'] = seed_case(rng)[1]
                args += [rng.choice(['--add', '--addprofile']), 'add.fa']
            elif kind == 'addprofile':
                args += ['--keeplength', '--add', 'seed1.fa'] if rng.random() < 0.5 else ['--mapout']
            elif kind == 'seedtable':
                args = ['--seedtable', 'tab.txt', '--seed', 'seed1.fa']
                files['tab.txt'] = '0 1 5 100.0 0 4 0 4 k\n'
            # the tree builders that the iterative refinement (which a seed forces) does not support stop with a message: a few of them
            opt = list(rng.choice(SEED_OPT_STOP if rng.random() < 0.12 else SEED_OPT))
            if rng.random() < 0.4:
                # the nodepair branch of tbfast with seeds (a pair library with --large / --memsavetree / --stepadd / --youngestlinkage) is registered as A2-19c-seed-nodepair
                opt += rng.choice([['--maxiterate', '2'], ['--localpair'], ['--genafpair', '--maxiterate', '2'], ['--globalpair', '--weighti', '0'], ['--retree', '1']])
            out.append(('seed-opt-%04d' % k, 'seed-opt', args + opt, text, files))
    if not want or 'seed-large' in want:
        # the nodepair branch of tbfast with seeds: a seed with --large / --memsavetree / --stepadd and a pair library (no iterative refinement)
        for k in range(int(n * 1.3)):
            ns = rng.choice([1, 1, 2, 3])
            seeds, text = seed_case(rng, ns)
            args = []
            files = {}
            anysym = rng.random() < 0.12
            for i, sd in enumerate(seeds):
                args += ['--seed', 'seed%d.fa' % (i + 1)]
                files['seed%d.fa' % (i + 1)] = odd_rows(rng, sd) if anysym else sd
            if anysym:
                args.append('--anysymbol')
                text = odd_rows(rng, text)
            args += [rng.choice(['--large', '--memsavetree', '--stepadd']), rng.choice(['--globalpair', '--localpair', '--genafpair', '--oldgenafpair'])]
            args += rng.choice([[], [], ['--weighti', '0'], ['--weighti', '1.5'], ['--retree', '1'], ['--treeout'], ['--distout'], ['--clustalout'], ['--reorder'], ['--op', '2', '--ep', '0.1'], ['--thread', '2'], ['--lhlimit', '3'], ['--unweight'], ['--minimumweight', '0.1']])
            out.append(('seed-large-%04d' % k, 'seed-large', args, text, files))
    return out


def run_case(case):
    """C.run_case without the runs that end in a segmentation fault of a program (the shell names a process id in its message, and a crash is not a result that can be compared)."""
    res, bad = C.run_case(case)
    if res is not None and 'Segmentation fault' in res['err']:
        return None, case[0]
    return res, bad


SHIFT_VALUES = ['1.0', '1.5', '2', '2.5', '3', '4.5', '9.9', '10', '20', '50', '99']
SHIFT_BASES = [[], [], [], ['--retree', '1'], ['--retree', '2'], ['--retree', '0'], ['--maxiterate', '1'], ['--maxiterate', '2'], ['--maxiterate', '1000'], ['--weighti', '0'], ['--weighti', '1'], ['--weighti', '3'], ['--treeout'], ['--distout'],
               ['--treeout', '--distout'], ['--reorder'], ['--clustalout'], ['--nofft'], ['--fft'], ['--6merpair'], ['--thread', '2', '--threadit', '0'], ['--amino'], ['--nuc'], ['--leavegappyregion'], ['--op', '2', '--ep', '0.1'], ['--outnum'],
               ['--unweight'], ['--minimumweight', '0.1'], ['--maxiterate', '2', '--skipiterate', '0.3'], ['--maxiterate', '2', '--bunkatsunashi'], ['--maxiterate', '2', '--core'], ['--anysymbol'], ['--maxiterate', '3', '--treeout'],
               ['--large'], ['--memsavetree'], ['--stepadd'], ['--large', '--weighti', '0'], ['--memsave'], ['--nwildcard'], ['--bl', '80'], ['--averagelinkage'], ['--lhlimit', '3', '--large'], ['--phylipout'], ['--inputorder']]


def build_shift(rng, want, n, out):
    """--shiftpenalty x with --globalpair (the only method the script allows): trywarp in A__align, G__align11 and partA__align when x < 10; the script's stops for the other methods."""
    if not want or 'shift' in want:
        for k in range(n):
            seqs, names, groups = make_case(rng, size_range=(2, 4))
            seqs = [x.replace('-', '') for x in seqs]
            x = rng.choice(SHIFT_VALUES)
            args = ['--globalpair', '--shiftpenalty', x] + list(rng.choice(SHIFT_BASES))
            if rng.random() < 0.1:
                args = ['--shiftpenalty', x] + rng.choice([['--localpair'], ['--genafpair'], ['--oldgenafpair'], [], ['--retree', '2'], ['--addfragments', 'frag.fa']])
            if rng.random() < 0.03:
                args = ['--globalpair', '--shiftpenalty', rng.choice(['0.5', '0', '1', '-2'])] + rng.choice([[], ['--retree', '1']])
            files = {'frag.fa': L.fasta(seqs[:2], ['f0', 'f1'])} if 'frag.fa' in args else {}
            out.append(('shift-%04d' % k, 'shift', args, L.fasta(seqs, names, rng.choice([None, None, 30])), files))
    return out


def main():
    rng = random.Random(7450515)
    want = set(sys.argv[1:])
    allcases = []
    refused = 0
    for c in build(rng, want) + build_shift(random.Random(7450516), want, int(260 * SCALE), []):
        nseq = c[3].count('\n>') + (1 if c[3].startswith('>') else 0)
        if TABLE in c[4] and table_overruns(c[4][TABLE], nseq):
            refused += 1
            continue
        allcases.append(c)
    print('refused (table behind the block):', refused, file=sys.stderr)
    L.run_all(allcases, run_case, OUT, os.environ.get('MC_TSV'), THREADS)


if __name__ == '__main__':
    main()
