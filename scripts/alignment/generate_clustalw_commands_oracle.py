#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-14b, B09 F104) oracle: ClustalW 2.1 guide-tree files and the row order.

Seqeron API: SequenceAligner.RunClustalW / RunClustalWProfile / RunClustalWSequences (ClustalWRunOptions: UseTree, UseTree1/2, GuideTreeOnly, OutputOrder).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7):
  align     `-ALIGN -NEWTREE=t.dnd -OUTPUT=FASTA -OUTORDER=ALIGNED|INPUT`   rows in the file order of the output + the guide-tree file
  newtree   `-NEWTREE=t.dnd`                                                  the guide-tree file only (no alignment)
  usetree   `-USETREE=u.dnd -OUTPUT=FASTA -OUTORDER=...`                      alignment along a user tree (valid trees of every shape, no distances, odd spellings, malformed files)
  profile   `-PROFILE1=a -PROFILE2=b -PROFILE [-USETREE1= -USETREE2=]`       rows + the two guide-tree files a.dnd / b.dnd (or the program's refusal)
  sequences `-PROFILE1=a -PROFILE2=b -SEQUENCES [-USETREE=]`                  rows + the guide-tree file b.dnd
with the option variants of the earlier ClustalW oracles (-GAPOPEN -GAPEXT -MAXDIV -NOWEIGHTS -NEGATIVE -QUICKTREE -MATRIX -PWMATRIX -TRANSWEIGHT -NOPGAP -NOHGAP -GAPDIST -PWGAPOPEN -PWGAPEXT).

Output TSV (tests/.../TestData/Alignment/clustalw_commands_oracle.tsv):
  id command type names(;) seqs(;) p1 options(k=v,k=v) outorder tree_in1 tree_in2 out_names(;) out_rows(;) tree_out1 tree_out2 error
Trees are written with \\n / \\t / \\\\ escaped. `error` = the program's refusal (WrongFormat | NotCompatible | NameNotFound | UnrootedWithoutDistances), Hang (killed after the timeout: the
program never terminates), Crash (killed by a signal) or empty.
Usage: python3 scripts/alignment/generate_clustalw_commands_oracle.py OUTFILE [setsPerKind] [seed]      (needs `clustalw` on PATH)
"""
import importlib.util, os, random, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dna = load('generate_clustalw_msa_full_oracle')
prot = load('generate_clustalw_protein_msa_oracle')
FLAGS = dict(prot.FLAGS)
FLAGS.update({'transweight': 'TRANSWEIGHT'})
TIMEOUT = 6
ODD = 0.0


def esc(t):
    return t.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t')


def cw(args, cwd):
    try:
        r = subprocess.run(['clustalw'] + args, capture_output=True, text=True, cwd=cwd, timeout=TIMEOUT)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return 'hang', '', ''


def read_fasta(path):
    names, rows = [], {}
    if not os.path.exists(path): return [], {}
    name = None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'): name = line[1:].strip(); names.append(name); rows[name] = ''
        elif name: rows[name] += line.strip()
    return names, rows


def read_text(path):
    return open(path).read() if os.path.exists(path) and os.path.getsize(path) > 0 else ''


def write_fa(path, names, seqs):
    with open(path, 'w') as f:
        for n, s in zip(names, seqs): f.write(f'>{n}\n{s}\n')


def opt_args(o):
    return [f'-{FLAGS[k]}' + ('' if v is None else f'={v}') for k, v in o.items()]


def ostr(o):
    return ','.join(k if v is None else f'{k}={v}' for k, v in o.items())


def classify(rc, err):
    if rc == 'hang': return 'Hang'
    if isinstance(rc, int) and rc < 0: return 'Crash'
    if isinstance(rc, int) and rc > 1 and rc != 255: return 'Crash'
    if 'Wrong format in tree file' in err: return 'WrongFormat'
    if 'unrooted and has no distances' in err: return 'UnrootedWithoutDistances'
    if 'not found' in err and 'not compatible' in err: return 'NameNotFound'
    if 'ERROR: tree' in err: return 'NotCompatible'
    return ''


# ----------------------------------------------------------------------------------------------------------------------- tree generators
def rand_names(rnd, k):
    pool = []
    style = rnd.random()
    if style < 0.5: return [f's{j:03d}' for j in range(k)]
    out = set()
    res = []
    while len(res) < k:
        n = ''.join(rnd.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-') for _ in range(rnd.randint(1, 14)))
        if n[0] in '-.': n = 'x' + n
        if n.lower() in out: continue
        out.add(n.lower()); res.append(n)
    return res


def fmt_dist(rnd, d):
    r = rnd.random()
    if r < 0.5: return f'{d:.5f}'
    if r < 0.7: return f'{d:.3f}'
    if r < 0.8: return repr(d)
    if r < 0.85: return f'{d:.2e}'
    if r < 0.9: return f'{d:.4f}'.lstrip('0') if d < 1 else f'{d:.4f}'
    return f'{d:.6g}'


def random_tree(rnd, names, rooted, distances, spaced):
    """random topology; rooted binary or unrooted trichotomy; returns Newick text"""
    leaves = list(names)
    rnd.shuffle(leaves)

    def build(ls):
        if len(ls) == 1:
            n = ls[0]
            return (n, True)
        cut = rnd.randint(1, len(ls) - 1)
        return ([build(ls[:cut]), build(ls[cut:])], False)

    def dist():
        return rnd.choice([0.0, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0]) * rnd.random() if rnd.random() < 0.8 else rnd.random() * 2

    def text(node, top=False):
        if node[1]:
            s = node[0]
            if distances: s += ':' + fmt_dist(rnd, dist())
            return s
        kids = node[0]
        inner = (',\n' if spaced and rnd.random() < 0.5 else ',').join(text(k) for k in kids)
        s = '(' + inner + ')'
        if distances and not top and rnd.random() < 0.9: s += ':' + fmt_dist(rnd, dist())
        return s

    if len(leaves) == 2:
        node = [(leaves[0], True), (leaves[1], True)], False
        return '(' + text(node[0][0]) + ',' + text(node[0][1]) + ');'
    if rooted:
        node = build(leaves)
        return text(node, top=True) + ';'
    # unrooted: three subtrees at the top
    c1 = rnd.randint(1, len(leaves) - 2)
    c2 = rnd.randint(c1 + 1, len(leaves) - 1)
    parts = [leaves[:c1], leaves[c1:c2], leaves[c2:]]
    return '(' + ','.join(text(build(p)) for p in parts) + ');'


def mutate_tree(rnd, text, names):
    """malformed / odd variants of a valid tree text"""
    kind = rnd.choice(['nolparen', 'empty', 'extraleaf', 'missingleaf', 'unknown', 'dup', 'poly4', 'single', 'garbage', 'badnum', 'trunc', 'rename', 'spaces', 'comment', 'noparens', 'swapcase', 'noterm', 'prefixspace', 'underscore'])
    if kind == 'nolparen': return kind, text[1:] if text.startswith('(') else 'x' + text
    if kind == 'empty': return kind, ''
    leaf_re = re.compile(r'[^(),:;\s]+')
    toks = leaf_re.findall(text)
    if kind == 'extraleaf':
        return kind, text.replace(')', ',zzextra:0.1)', 1) if ')' in text else text
    if kind == 'missingleaf':
        n = rnd.choice(names)
        return kind, re.sub(r'(?<![\w.\-])' + re.escape(n) + r'(:[^,()]*)?,?', '', text, count=1)
    if kind == 'unknown':
        n = rnd.choice(names)
        return kind, text.replace(n, 'nosuchname', 1)
    if kind == 'dup':
        a, b = rnd.sample(names, 2) if len(names) > 1 else (names[0], names[0])
        return kind, text.replace(b, a, 1)
    if kind == 'poly4':
        a = rnd.choice(names)
        return kind, text.replace(a, a + ',dup4x:0.1,dup5x:0.2', 1)
    if kind == 'single':
        a = rnd.choice(names)
        return kind, text.replace(a, '(' + a + ')', 1)
    if kind == 'garbage': return kind, text + '\ngarbage here'
    if kind == 'badnum':
        return kind, re.sub(r':[0-9.eE+-]+', ':abc', text, count=1)
    if kind == 'trunc':
        return kind, text[:rnd.randint(1, max(1, len(text) - 1))]
    if kind == 'rename':
        a = rnd.choice(names)
        return kind, text.replace(a, a.upper(), 1)
    if kind == 'spaces':
        return kind, re.sub(r',', ' , ', text).replace('(', '( ').replace(':', ' : ')
    if kind == 'comment':
        return kind, text.replace(',', '[c],', 1)
    if kind == 'noparens':
        return kind, ','.join(names) + ';'
    if kind == 'swapcase':
        return kind, ''.join(c.swapcase() if rnd.random() < 0.5 else c for c in text)
    if kind == 'noterm':
        return kind, text.rstrip(';')
    if kind == 'prefixspace':
        return kind, '  \n\t' + text
    if kind == 'underscore':
        return kind, text.replace(names[0], names[0].replace('_', ' '), 1)
    return kind, text


def odd_tree(rnd, names):
    """random arities (1..4 children), leaves from the names (and an occasional empty / odd leaf): the shapes the reader of the program is not meant for"""
    leaves = list(names)
    rnd.shuffle(leaves)

    def build(ls):
        if len(ls) == 1 and rnd.random() < 0.85:
            return ls[0] + (':' + fmt_dist(rnd, rnd.random() * 0.3) if rnd.random() < 0.9 else '')
        if len(ls) == 1:
            ls = ls + ['']
        arity = rnd.choice([1, 2, 2, 2, 3, 3, 4]) if len(ls) > 1 else 1
        arity = min(arity, len(ls))
        cuts = sorted(rnd.sample(range(1, len(ls)), arity - 1)) if arity > 1 else []
        parts = [ls[a:b] for a, b in zip([0] + cuts, cuts + [len(ls)])]
        return '(' + ','.join(build(p) for p in parts) + ')' + (':' + fmt_dist(rnd, rnd.random() * 0.3) if rnd.random() < 0.7 else '')

    return build(leaves).rstrip(',') + ';'


def user_tree(rnd, names, own):
    r = rnd.random()
    k = len(names)
    if rnd.random() < ODD: return odd_tree(rnd, names)
    if r < 0.12 and own: return own
    if r < 0.45:
        t = random_tree(rnd, names, rooted=True, distances=True, spaced=True)
    elif r < 0.65:
        t = random_tree(rnd, names, rooted=False, distances=True, spaced=True) if k >= 3 else random_tree(rnd, names, True, True, False)
    elif r < 0.75:
        t = random_tree(rnd, names, rooted=True, distances=False, spaced=False)
    elif r < 0.8:
        t = random_tree(rnd, names, rooted=False, distances=False, spaced=False) if k >= 3 else random_tree(rnd, names, True, False, False)
    else:
        base = own if (own and rnd.random() < 0.5) else random_tree(rnd, names, rnd.random() < 0.5, True, False)
        _, t = mutate_tree(rnd, base, names)
    return t


# ----------------------------------------------------------------------------------------------------------------------- cases
def make_options(rnd, typ):
    # the option generators of the earlier oracles draw from their own module-level random state: seed it from ours
    m = dna if typ == 'dna' else prot
    m.rnd.seed(rnd.randint(0, 2 ** 31))
    o = m.options()
    if typ == 'protein' and 'novgap' in o: o.pop('novgap')
    if 'endgaps' in o: o.pop('endgaps')
    if 'hgapresidues' in o: o.pop('hgapresidues')
    if typ == 'dna':
        o.pop('ktuple', None); o.pop('window', None); o.pop('topdiags', None); o.pop('pairgap', None)
    else:
        for k in ('ktuple', 'window', 'topdiags', 'pairgap'): o.pop(k, None)
    return o


def make_seqs(rnd, typ):
    maxlen = rnd.choice([40, 60, 90, 140])
    kmax = rnd.choice([4, 6, 8, 10, 14])
    dna.rnd.seed(rnd.randint(0, 2 ** 31)); prot.rnd.seed(rnd.randint(0, 2 ** 31))
    seqs = dna.make_set(maxlen, 'ACGT', 3, kmax) if typ == 'dna' else prot.make_set(maxlen, 3, kmax)
    return seqs


def case_main(rnd, i, command, typ):
    seqs = make_seqs(rnd, typ)
    if rnd.random() < 0.08: seqs = seqs[:2]
    names = rand_names(rnd, len(seqs))
    o = make_options(rnd, typ)
    outorder = rnd.choice(['aligned', 'aligned', 'input'])
    t = typ.upper()
    with tempfile.TemporaryDirectory() as d:
        write_fa(os.path.join(d, 'x.fa'), names, seqs)
        base = [f'-INFILE=x.fa', f'-TYPE={t}'] + opt_args(o)
        tree_in = ''
        if command == 'align':
            rc, out, err = cw(base + ['-ALIGN', '-OUTPUT=FASTA', f'-OUTORDER={outorder.upper()}', '-OUTFILE=o.fa', '-NEWTREE=t.dnd'], d)
        elif command == 'newtree':
            rc, out, err = cw(base + ['-NEWTREE=t.dnd'], d)
        else:
            # usetree: the program's own tree (written by a first run) or a generated / mutated one
            cw(base + ['-NEWTREE=own.dnd'], d)
            own = read_text(os.path.join(d, 'own.dnd'))
            tree_in = user_tree(rnd, names, own if own else None)
            with open(os.path.join(d, 'u.dnd'), 'w', newline='') as f: f.write(tree_in)
            rc, out, err = cw(base + ['-USETREE=u.dnd', '-OUTPUT=FASTA', f'-OUTORDER={outorder.upper()}', '-OUTFILE=o.fa'], d)
        error = classify(rc, err)
        on, rows = read_fasta(os.path.join(d, 'o.fa'))
        tree_out = read_text(os.path.join(d, 't.dnd'))
        if command == 'usetree' and not on: error = error or 'Unclassified'
        if command == 'align' and not on: error = error or 'Unclassified'
    if error in ('Hang', 'Crash') or error == 'Unclassified' or (command != 'newtree' and not on and not error):
        on, rows = [], {}
    if error: on, rows = [], {}
    return [i, command, typ, ';'.join(names), ';'.join(seqs), 0, ostr(o), outorder, esc(tree_in), '', ';'.join(on), ';'.join(rows[n] for n in on), esc(tree_out), '', error]


def make_profiles(rnd, typ):
    seqs = make_seqs(rnd, typ)
    if len(seqs) < 2: seqs = seqs + seqs[:1]
    names = rand_names(rnd, len(seqs))
    cut = rnd.randint(1, len(seqs) - 1)
    return names, seqs, cut


def case_profile(rnd, i, command, typ):
    names, seqs, cut = make_profiles(rnd, typ)
    o = make_options(rnd, typ)
    if command == 'sequences':
        for k in ('pwmatrix', 'pwgapopen', 'pwgapext'):
            pass
    outorder = rnd.choice(['aligned', 'aligned', 'input'])
    t = typ.upper()
    with tempfile.TemporaryDirectory() as d:
        def align_part(part_names, part_seqs, tag):
            write_fa(os.path.join(d, tag + '0.fa'), part_names, part_seqs)
            if len(part_seqs) == 1: return list(part_seqs)
            cw([f'-INFILE={tag}0.fa', f'-TYPE={t}', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={tag}0.aln.fa'], d)
            n, r = read_fasta(os.path.join(d, tag + '0.aln.fa'))
            return [r[x] for x in part_names]
        p1 = align_part(names[:cut], seqs[:cut], 'a')
        if command == 'profile':
            p2 = align_part(names[cut:], seqs[cut:], 'b')
        else:
            p2 = align_part(names[cut:], seqs[cut:], 'b') if rnd.random() < 0.25 else seqs[cut:]
        write_fa(os.path.join(d, 'a.fa'), names[:cut], p1)
        write_fa(os.path.join(d, 'b.fa'), names[cut:], p2)
        for f in ('a.dnd', 'b.dnd'):
            if os.path.exists(os.path.join(d, f)): os.remove(os.path.join(d, f))
        flag = '-PROFILE' if command == 'profile' else '-SEQUENCES'
        args = ['-PROFILE1=a.fa', '-PROFILE2=b.fa', flag, f'-TYPE={t}', '-OUTPUT=FASTA', f'-OUTORDER={outorder.upper()}', '-OUTFILE=o.fa'] + opt_args(o)
        tin1 = tin2 = ''
        mode = rnd.random()
        if command == 'profile':
            n1, n2 = names[:cut], names[cut:]
            if mode < 0.35:   # user trees (one or both, valid or mutated)
                if rnd.random() < 0.8 or len(n1) < 2:
                    if len(n2) >= 1 and (rnd.random() < 0.6 or len(n1) < 2): tin2 = (user_tree(rnd, n2, None) if len(n2) >= 2 else random_tree(rnd, n2 * 2, True, True, False)) or '\n'
                if rnd.random() < 0.6 and len(n1) >= 1:
                    tin1 = (user_tree(rnd, n1, None) if len(n1) >= 2 else random_tree(rnd, n1 * 2, True, True, False)) or '\n'
                if tin1: open(os.path.join(d, 'u1.dnd'), 'w').write(tin1); args.append('-USETREE1=u1.dnd')
                if tin2: open(os.path.join(d, 'u2.dnd'), 'w').write(tin2); args.append('-USETREE2=u2.dnd')
            elif mode < 0.5:
                args.append('-NEWTREE1=x1.dnd' if rnd.random() < 0.5 else '-NEWTREE2=x2.dnd')
        else:
            if mode < 0.3:
                tin1 = user_tree(rnd, names, None) or '\n'
                open(os.path.join(d, 'u1.dnd'), 'w').write(tin1); args.append('-USETREE=u1.dnd')
            elif mode < 0.45:
                args.append('-NEWTREE=x1.dnd')
        rc, out, err = cw(args, d)
        error = classify(rc, err)
        on, rows = read_fasta(os.path.join(d, 'o.fa'))
        t1 = read_text(os.path.join(d, 'a.dnd')) if command == 'profile' else ''
        t2 = read_text(os.path.join(d, 'b.dnd'))
        if command == 'sequences':
            t1 = read_text(os.path.join(d, 'x1.dnd')) or t2
            t2 = ''
        trees_only = any(a.startswith('-NEWTREE') for a in args)
        if not on and not error and not trees_only: error = 'Unclassified'
        if error: on, rows = [], {}
    return [i, command + ('-trees' if trees_only else ''), typ, ';'.join(names), ';'.join(p1 + p2), cut, ostr(o), outorder, esc(tin1), esc(tin2), ';'.join(on), ';'.join(rows[n] for n in on), esc(t1), esc(t2), error]


def main(outfile, per, seed):
    global ODD
    ODD = float(os.environ.get('CW_ODD', '0.12'))
    rnd = random.Random(seed)
    rows = []
    i = 0
    plan = [('align', 'dna'), ('align', 'protein'), ('newtree', 'dna'), ('newtree', 'protein'), ('usetree', 'dna'), ('usetree', 'protein'),
            ('profile', 'dna'), ('profile', 'protein'), ('sequences', 'dna'), ('sequences', 'protein')]
    for command, typ in plan:
        for _ in range(per):
            rows.append((case_main if command in ('align', 'newtree', 'usetree') else case_profile)(rnd, i, command, typ))
            i += 1
    with open(outfile, 'w') as f:
        f.write('id\tcommand\ttype\tnames\tseqs\tp1\toptions\toutorder\ttree_in1\ttree_in2\tout_names\tout_rows\ttree_out1\ttree_out2\terror\n')
        for r in rows: f.write('\t'.join(str(x) for x in r) + '\n')
    from collections import Counter
    print(len(rows), 'sets;', Counter((r[1], r[-1]) for r in rows))


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 330, int(sys.argv[3]) if len(sys.argv) > 3 else 20261121)
