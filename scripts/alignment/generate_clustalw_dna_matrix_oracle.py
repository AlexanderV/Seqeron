#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-10d part 3) oracle: ClustalW 2.1 nucleotide matrices `-DNAMATRIX=IUB|CLUSTALW|file`, `-PWDNAMATRIX=IUB|CLUSTALW|file`, `-TRANSWEIGHT`.

Seqeron API: ClustalWOptions.DnaMatrix / PairwiseDnaMatrix / UserDnaMatrix / PairwiseUserDnaMatrix / UserDnaMatrixUse, ClustalWDnaPairwiseOptions,
SequenceAligner.MultipleAlignClustalW(seqs, ClustalWOptions), SequenceAligner.ClustalWPairwiseAlign(a, b, ClustalWDnaPairwiseOptions).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7), two ways of giving the matrices:
  * mode `cli`  : `clustalw -TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT [-DNAMATRIX=x] [-PWDNAMATRIX=y] [options]`. x / y are IUB, CLUSTALW or a user file. For a FILE the
                  command line reads and validates it (SubMatrix::getUserMatFromFile -> readUserMatrix) and then registers its name under the PROTEIN slot
                  (CommandLineParser: setCurrentNameAndNum(name, 3, Protein, ...)): the DNA stages keep IUB. A malformed file ends the program (exit 1, no output): rows = EXIT.
  * mode `menu` : the interactive menu (`clustalw -INTERACTIVE -TYPE=DNA`: without the explicit type the menu's sequence reader classifies a set with >= 15 % non-ACGTUN letters as
                  PROTEIN) driven through stdin (multiple-alignment parameters -> "6. DNA weight matrix" -> "3. User defined" -> file; pairwise parameters ->
                  "4. DNA weight matrix" -> "3. User defined" -> file; the gap penalties / transition weight / delay-divergent / negative-matrix items; "1. Do complete multiple alignment").
                  The menu registers the file for DNA, so the matrix IS used. (The menu's built-in IUB / CLUSTALW items register an empty name - a second defect - and are not used.)
The matrix files are produced by this script in several layouts: EMBOSS / NCBI (comments, labelled rows, '*' row and column), unlabelled, integer / floating-point rows (x10 per row),
tab separated, junk in the upper triangle, blank / comment / trailer lines; header orders (IUB packing order, ATGC..., alphabetical, shuffled), 13-17 letters (the transition weights
of getProfileAlignMatrix index the matrix with Xref[0], [2], [4], [11], [12], i.e. the letters at those header positions) and short headers (4-12 letters: Xref -1, a write outside the
C++ array - the installed binary crashed with SIGSEGV on 8 of 37 such menu runs; those are not in the oracle, the others are recorded).

Output TSV (tests/.../TestData/Alignment/clustalw_dna_matrix_oracle.tsv):
  id  mode  seqs(;)  options(k=v,...)  files(JSON name -> content)  rows(;|EXIT)  pair_scores(i:j:score;...)  group_scores(g:entries:score | g:D ;...)
Usage: python3 scripts/alignment/generate_clustalw_dna_matrix_oracle.py OUTDIR [nSets] [seed]      (needs `clustalw` on PATH)
"""
import json, os, random, re, subprocess, sys, tempfile, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_clustalw_msa_full_oracle as g

rnd = random.Random(20261129)
g.rnd.seed(20261129)

IUB = "ABCDGHKMNRSTUVWXY"            # the packing order of matrices.h
ORDERS = {
    'iub': IUB,
    'atgc': "ATGCSWRYKMBVHDNU",
    'alpha': "ABCDGHKMNRSTUVWXY",
    'acgt_amb': "ACGTRYKMSWBDHVN",
}
AMBIG = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'U': 'T', 'R': 'AG', 'Y': 'CT', 'K': 'GT', 'M': 'AC', 'S': 'CG', 'W': 'AT', 'B': 'CGT', 'D': 'AGT',
         'H': 'ACT', 'V': 'ACG', 'N': 'ACGT', 'X': 'ACGT'}


def score_pair(a, b, scheme):
    kind = scheme['kind']
    if kind == 'match':
        return scheme['m'] if a == b else scheme['x']
    ma, mb = AMBIG[a], AMBIG[b]
    if kind == 'overlap':     # IUB-like: any overlap of the sets
        return scheme['m'] if set(ma) & set(mb) else scheme['x']
    if kind == 'ednafull':    # fraction-weighted
        if len(ma) == 1 and len(mb) == 1:
            return scheme['m'] if ma == mb else scheme['x']
        ov = len(set(ma) & set(mb))
        return round(scheme['m'] * ov / (len(ma) * len(mb)) * 2 + scheme['x'] * (1 - ov / (len(ma) * len(mb))) * 0.5, 1) if ov else scheme['x'] // 2
    if kind == 'ts':          # transition / transversion
        if ma == mb and len(ma) == 1:
            return scheme['m']
        if len(ma) == 1 and len(mb) == 1:
            return scheme['ts'] if {ma, mb} in ({'A', 'G'}, {'C', 'T'}) else scheme['x']
        return scheme['amb'] if set(ma) & set(mb) else scheme['x']
    return float(rnd.randint(-6, 9))


def random_scheme():
    k = rnd.choice(['match', 'match', 'overlap', 'ednafull', 'ts', 'ts', 'random'])
    m = rnd.choice([1, 2, 4, 5, 10, 15, 20])
    return dict(kind=k, m=m, x=rnd.choice([-1, -2, -3, -4, -5, -9, 0, -m]), ts=rnd.choice([-1, 0, 1, 2, -2]), amb=rnd.choice([0, 1, 2, -1]))


def pick_letters():
    q = rnd.random()
    if q < 0.55:   # at least 13 letters: Xref[0], [2], [4], [11], [12] exist
        base = rnd.choice(list(ORDERS.values()))
        letters = list(base)
        if rnd.random() < 0.25:
            letters = rnd.sample(letters, len(letters))
        if len(letters) < 13:
            letters += [c for c in "ABCDGHKMNRSTUVWXY" if c not in letters][:13 - len(letters)]
        return letters
    if q < 0.80:   # short headers
        n = rnd.choice([4, 5, 6, 8, 10, 12])
        pool = list("ACGT") + rnd.sample(list("NUBDHKMRSVWXY"), 9)
        letters = list("ACGT") if rnd.random() < 0.5 else rnd.sample(pool, n)
        letters = letters[:n] if len(letters) > n else letters + [c for c in pool if c not in letters][:max(0, n - len(letters))]
        if rnd.random() < 0.4:
            letters = rnd.sample(letters, len(letters))
        return letters
    n = rnd.randint(13, 17)   # 13-17 letters of a random subset
    return rnd.sample(list(IUB), n)


def fmt_value(v, mode):
    if mode == 'int':
        return str(int(round(v)))
    if mode == 'f1':
        return '%.1f' % (v / 2.0)
    return '%g' % (v / 2.0)


def render(scheme, letters, layout):
    mode, sep, star = layout['mode'], layout['sep'], layout['star']
    out = []
    if layout['comments']:
        out += ['# DNA matrix made for the user-matrix oracle', '# Entries for the matrix at a scale of 1/2 bit', '']
    header = sep.join(list(letters) + [star])
    out.append((sep + header) if layout['indent'] else header)
    rows = list(letters) + [star]
    cache = {}
    for i, a in enumerate(rows):
        cells = []
        for j, b in enumerate(rows):
            if a == star or b == star:
                v = -4.0
            elif j <= i or not layout['junk']:
                key = (min(a, b), max(a, b))
                if scheme['kind'] == 'random':
                    v = cache.setdefault(key, float(rnd.randint(-6, 9)))
                else:
                    v = float(score_pair(a, b, scheme))
            else:
                v = float(rnd.randint(-9, 9))
            cells.append(fmt_value(v, mode))
        if layout['labels']:
            cells = [a] + cells
        out.append(sep.join(cells))
        if layout['blank_lines'] and rnd.random() < 0.2 and i < len(rows) - 1:
            out.append('')
    if layout['trailer']:
        out.append(rnd.choice(['# end of matrix', '! end', '#']))
        out.append('1 2 3 this is ignored')
    return '\n'.join(out) + '\n'


def random_layout():
    return dict(mode=rnd.choice(['int', 'int', 'int', 'f1', 'f1', 'g']), sep=rnd.choice([' ', '  ', ' ', '\t']), comments=rnd.random() < 0.6,
                labels=rnd.random() < 0.75, indent=rnd.random() < 0.5, junk=rnd.random() < 0.15, blank_lines=rnd.random() < 0.2,
                trailer=rnd.random() < 0.2, star=rnd.choice(['*', '*', '*', '-', '.']))


def make_matrix_file():
    return render(random_scheme(), pick_letters(), random_layout())


def malformed_file():
    kind = rnd.choice(['series', 'norows', 'short', 'empty', 'two'])
    if kind == 'series':
        return 'CLUSTAL_SERIES\nMATRIX 0 100 other.mat\n'
    if kind == 'norows':      # no '*' row
        t = make_matrix_file().rstrip('\n').split('\n')
        return '\n'.join(t[:-1]) + '\n'
    if kind == 'short':       # a row with too few values
        t = make_matrix_file().rstrip('\n').split('\n')
        i = max(2, len(t) // 2)
        t[i] = ' '.join(t[i].split()[:2])
        return '\n'.join(t) + '\n'
    if kind == 'empty':
        return '# nothing here\n\n'
    return 'A C G T *\nA 1 0 0 0 0\nC 0 1 0 0 0\nG 0 0 1 0 0\nT 0 0 0 1 0\n'


def common_options():
    o = {}
    if rnd.random() < 0.35: o['gapopen'] = rnd.choice([5, 10, 15, 20, 25, 7.5])
    if rnd.random() < 0.35: o['gapext'] = rnd.choice([0.5, 1, 3, 6.66, 8])
    if rnd.random() < 0.45: o['transweight'] = rnd.choice([0, 0.25, 0.5, 0.75, 1])
    if rnd.random() < 0.25: o['maxdiv'] = rnd.choice([0, 10, 30, 50, 70, 100])
    if rnd.random() < 0.2: o['negative'] = None
    if rnd.random() < 0.25: o['quick'] = None
    if rnd.random() < 0.25: o['pwgapopen'] = rnd.choice([5, 10, 15, 20])
    if rnd.random() < 0.25: o['pwgapext'] = rnd.choice([1, 3.33, 6.66, 9])
    return o


FLAGS = {'gapopen': 'GAPOPEN', 'gapext': 'GAPEXT', 'transweight': 'TRANSWEIGHT', 'maxdiv': 'MAXDIV', 'noweights': 'NOWEIGHTS', 'negative': 'NEGATIVE',
         'quick': 'QUICKTREE', 'pwgapopen': 'PWGAPOPEN', 'pwgapext': 'PWGAPEXT', 'dnamatrix': 'DNAMATRIX', 'pwdnamatrix': 'PWDNAMATRIX'}


def parse_stdout(text):
    pairs = [(a, b, s) for a, b, s in re.findall(r'Sequences \((\d+):(\d+)\) Aligned\. Score:\s+(-?\d+)', text)]
    groups = []
    for m in re.finditer(r'Group (\d+):\s+(?:Sequences:\s*(\d+)\s+Score:(-?\d+)|\s*Delayed)', text):
        groups.append(f'{m.group(1)}:{m.group(2)}:{m.group(3)}' if m.group(2) else f'{m.group(1)}:D')
    return pairs, groups


def read_fasta(path, n):
    rows, name = {}, None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'): name = line[1:].strip(); rows[name] = ''
        elif name: rows[name] += line.strip()
    return [rows[f's{i}'] for i in range(n)]


def read_aln(path, n):
    rows = {f's{i}': '' for i in range(n)}
    for line in open(path):
        p = line.rstrip('\n').split()
        if len(p) == 2 and p[0] in rows:
            rows[p[0]] += p[1]
    return [rows[f's{i}'] for i in range(n)]


def run_cli(seqs, o, files):
    with tempfile.TemporaryDirectory() as d:
        for name, text in files.items():
            open(os.path.join(d, name), 'w', newline='').write(text)
        fa = os.path.join(d, 'x.fa')
        open(fa, 'w').write(''.join(f'>s{i}\n{s}\n' for i, s in enumerate(seqs)))
        out = os.path.join(d, 'o.fa')
        args = ['clustalw', f'-INFILE={fa}', '-TYPE=DNA', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}']
        for k, v in o.items():
            args.append(f'-{FLAGS[k]}' + ('' if v is None else f'={v}'))
        p = subprocess.run(args, capture_output=True, text=True, cwd=d, timeout=120)
        if not os.path.exists(out):
            return None, p.stdout + p.stderr
        return (read_fasta(out, len(seqs)),) + parse_stdout(p.stdout) + (p.stdout,)


def menu_script(o, matrix, pwmatrix, fasta_name):
    L = ['1', fasta_name, '2']
    # multiple alignment parameters (6), pairwise parameters (5)
    multi = []
    if 'gapopen' in o: multi += ['1', str(o['gapopen'])]
    if 'gapext' in o: multi += ['2', str(o['gapext'])]
    if 'maxdiv' in o: multi += ['3', str(o['maxdiv'])]
    if 'transweight' in o: multi += ['4', str(o['transweight'])]
    if matrix: multi += ['6', '3', matrix, '']
    if 'negative' in o: multi += ['7']
    if multi:
        L += ['6'] + multi + ['']
    pw = []
    if 'pwgapopen' in o: pw += ['1', str(o['pwgapopen'])]
    if 'pwgapext' in o: pw += ['2', str(o['pwgapext'])]
    if pwmatrix: pw += ['4', '3', pwmatrix, '']
    if pw:
        L += ['5'] + pw + ['']
    if 'quick' in o:
        L += ['4']
    L += ['1', '', '', '', '', 'X']
    return '\n'.join(L) + '\n'


def run_menu(seqs, o, files, matrix, pwmatrix):
    with tempfile.TemporaryDirectory() as d:
        for name, text in files.items():
            open(os.path.join(d, name), 'w', newline='').write(text)
        fa = os.path.join(d, 'x.fa')
        open(fa, 'w').write(''.join(f'>s{i}\n{s}\n' for i, s in enumerate(seqs)))
        script = menu_script(o, matrix, pwmatrix, 'x.fa')
        p = subprocess.Popen(['clustalw', '-INTERACTIVE', '-TYPE=DNA'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=d)
        p.stdin.write(script.encode()); p.stdin.close()
        buf = b''
        t0 = time.time()
        while time.time() - t0 < 60:
            chunk = p.stdout.read1(65536) if hasattr(p.stdout, 'read1') else p.stdout.read(65536)
            if not chunk:
                break
            buf += chunk
            if len(buf) > 3_000_000:
                break
        rc = p.poll()
        p.kill(); p.wait()
        text = buf.decode(errors='replace')
        aln = os.path.join(d, 'x.aln')
        if not os.path.exists(aln) or 'Alignment Score' not in text:
            return None, f'rc={rc} script={script!r} files={files!r} tail={text[-400:]!r}'
        return (read_aln(aln, len(seqs)),) + parse_stdout(text) + (text,)


def make_seqs():
    kinds = [(0.6, 60), (0.3, 120), (0.1, 250)]
    r = rnd.random(); acc = 0
    for pr, ml in kinds:
        acc += pr
        if r < acc: break
    seqs = g.make_set(ml, g.ALPHA, 2, rnd.choice([4, 6, 8, 9, 12]))
    if rnd.random() < 0.15: seqs = g.ambiguous(seqs)
    return seqs


def main(outdir, nsets):
    os.makedirs(outdir, exist_ok=True)
    done, tries = 0, 0
    with open(os.path.join(outdir, 'clustalw_dna_matrix_oracle.tsv'), 'w') as f:
        f.write('id\tmode\tseqs\toptions\tfiles\trows\tpair_scores\tgroup_scores\n')
        while done < nsets:
            tries += 1
            seqs = make_seqs()
            o = common_options()
            files = {}
            mode = rnd.choice(['cli', 'cli', 'menu', 'menu', 'menu', 'menu'])
            if mode == 'cli':
                variant = rnd.choice(['builtin', 'builtin', 'file', 'file', 'file', 'mixed', 'bad'])
                if variant in ('builtin', 'mixed'):
                    o['dnamatrix'] = rnd.choice(['IUB', 'CLUSTALW', 'CLUSTALW', 'iub', 'clustalw'])
                    o['pwdnamatrix'] = rnd.choice(['IUB', 'CLUSTALW', 'CLUSTALW', 'iub', 'clustalw'])
                    if rnd.random() < 0.3: del o['pwdnamatrix']
                    if variant == 'mixed':
                        files['d.mat'] = make_matrix_file(); o['pwdnamatrix'] = 'd.mat'
                elif variant == 'file':
                    if rnd.random() < 0.75: files['m.mat'] = make_matrix_file(); o['dnamatrix'] = 'm.mat'
                    if rnd.random() < 0.6 or not files: files['p.mat'] = make_matrix_file(); o['pwdnamatrix'] = 'p.mat'
                else:
                    files['bad.mat'] = malformed_file(); o['dnamatrix' if rnd.random() < 0.5 else 'pwdnamatrix'] = 'bad.mat'
                res = run_cli(seqs, o, files)
            else:
                variant = 'menu'
                m = pw = None
                kind = rnd.choice(['both', 'both', 'profile', 'pairwise'])
                if kind in ('both', 'profile'): files['m.mat'] = make_matrix_file(); m = 'm.mat'
                if kind in ('both', 'pairwise'): files['p.mat'] = make_matrix_file(); pw = 'p.mat'
                if kind == 'both' and rnd.random() < 0.3: files['p.mat'] = files['m.mat']
                o.pop('noweights', None)
                res = run_menu(seqs, o, files, m, pw)
                o_store = dict(o)
                if m: o_store['dnamatrix'] = m
                if pw: o_store['pwdnamatrix'] = pw
                o = o_store
            rows = res[0] if res and res[0] is not None else None
            if rows is None and variant != 'bad':
                print('run failed on try', tries, mode, variant, res[-1] if res else ''); continue
            ostr = ','.join(k if v is None else f'{k}={v}' for k, v in o.items())
            if rows is None:
                f.write(f"{done}\t{mode}\t{';'.join(seqs)}\t{ostr}\t{json.dumps(files)}\tEXIT\t\t\n")
            else:
                assert all(len(x) == len(rows[0]) for x in rows), (tries, o)
                ps = ';'.join(f'{a}:{b}:{s}' for a, b, s in res[1])
                f.write(f"{done}\t{mode}\t{';'.join(seqs)}\t{ostr}\t{json.dumps(files)}\t{';'.join(rows)}\t{ps}\t{';'.join(res[2])}\n")
            done += 1
    print(done, 'sets in', tries, 'tries')


if __name__ == '__main__':
    if len(sys.argv) > 3:
        rnd.seed(int(sys.argv[3])); g.rnd.seed(int(sys.argv[3]))
    main(sys.argv[1] if len(sys.argv) > 1 else '.', int(sys.argv[2]) if len(sys.argv) > 2 else 330)
