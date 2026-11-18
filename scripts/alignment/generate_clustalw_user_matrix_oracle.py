#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A1-10d, file part) oracle: ClustalW 2.1 user-defined matrix files (`-MATRIX=file`, `-PWMATRIX=file`, matrix SERIES files).

Seqeron API: ClustalWUserMatrix.Parse / ClustalWMatrixFile.Parse + ClustalWProteinOptions.UserMatrix / UserMatrixSeries and
ClustalWProteinPairwiseOptions.UserMatrix, SequenceAligner.MultipleAlignClustalW(...) and SequenceAligner.ClustalWPairwiseAlign(...).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7) `-TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT` with matrix files written by this script
(SubMatrix::readUserMatrix / readMatrixSeries). The files are produced from the Biopython 1.88 substitution matrices (BLOSUM45/50/62/80/90, PAM30/70/250, GONNET1992,
JONES, DAYHOFF, ...) and random matrices in several layouts: NCBI / EMBOSS (comment lines, labelled rows, '*' row and column), unlabelled rows, integer rows, floating
point rows (values x10; "%.1f" or "%g" which makes the scale per row), tab separated, header subsets (unlisted residues score 0), shuffled residue order, junk in the
upper triangle (ignored), a trailing '!' / '#' line, extra blank / comment lines; and CLUSTAL_SERIES files with 2-4 "MATRIX llimit ulimit file" entries (with gaps and
overlaps in the limits), and the quirks of readUserMatrix (a header letter ClustalW does not know - J or lower case - only warns and shifts the stored values; matrices with two or three residues). The pairwise "Sequences (i:j) Aligned. Score:" lines are recorded as well (pairwise stage with -PWMATRIX=file).

Output TSV (tests/.../TestData/Alignment/clustalw_user_matrix_oracle.tsv):  id  seqs(;)  rows(;)  options(k=v,k=v)  files(JSON name -> content)  pair_scores(;)
Usage: python3 scripts/alignment/generate_clustalw_user_matrix_oracle.py OUTDIR [nSets] [seed]      (needs `clustalw` on PATH and Biopython)
"""
import json, os, random, re, subprocess, sys, tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import generate_clustalw_protein_msa_oracle as g

rnd = random.Random(20261127)
g.rnd.seed(20261127)

from Bio.Align import substitution_matrices as sm

LIB = ['BLOSUM45', 'BLOSUM50', 'BLOSUM62', 'BLOSUM80', 'BLOSUM90', 'PAM30', 'PAM70', 'PAM250', 'GONNET1992', 'JONES', 'DAYHOFF', 'MDM78', 'RISLER', 'JOHNSON']
AA20 = "ACDEFGHIKLMNPQRSTVWY"
CODES = "ABCDEFGHIKLMNOPQRSTUVWXYZ"


def lib_matrix(name):
    m = sm.load(name)
    letters = [c for c in m.alphabet if c in CODES]
    return {(a, b): float(m[a, b]) for a in letters for b in letters}, letters


def random_matrix():
    n = rnd.choice([8, 12, 20, 20, 20, 23])
    letters = rnd.sample(AA20, min(n, 20)) + (list("BZX") if n == 23 else [])
    d = {}
    for i, a in enumerate(letters):
        for j in range(i + 1):
            b = letters[j]
            v = rnd.randint(3, 9) if i == j else rnd.randint(-5, 3)
            d[(a, b)] = d[(b, a)] = float(v)
    return d, letters


def pick_matrix():
    q = rnd.random()
    if q < 0.10:   # quirks of readUserMatrix: unknown letters (J, lower case) only warn and misread the table; tiny matrices
        d, letters = lib_matrix(rnd.choice(LIB))
        kind = rnd.choice(['j_end', 'j_mid', 'two', 'three', 'lower', 'j_mid'])
        if kind == 'j_end':
            letters = letters + ['J']
        elif kind == 'j_mid':
            letters = list(letters)
            letters.insert(rnd.randint(0, len(letters)), 'J')
        elif kind == 'two':
            letters = rnd.sample(letters, 2)
        elif kind == 'three':
            letters = rnd.sample(letters, 3)
        else:
            d = {(a.lower(), b.lower()): v for (a, b), v in d.items()}
            letters = [c.lower() for c in letters]
        return d, letters
    if rnd.random() < 0.7:
        d, letters = lib_matrix(rnd.choice(LIB))
        if rnd.random() < 0.3:   # header subset
            letters = rnd.sample(letters, rnd.randint(max(3, len(letters) // 2), len(letters)))
        elif rnd.random() < 0.3:  # shuffled residue order
            letters = rnd.sample(letters, len(letters))
        return d, letters
    return random_matrix()


def fmt_value(v, mode):
    if mode == 'int':
        return str(int(round(v)))
    if mode == 'f1':
        return '%.1f' % (v / 2.0)
    return '%g' % (v / 2.0)   # 'g': whole numbers lose the '.', so the x10 scale is decided per row


def render(d, letters, layout):
    """The text of a ClustalW matrix file; layout is a dict of formatting switches."""
    mode, sep, star = layout['mode'], layout['sep'], layout['star']
    out = []
    if layout['comments']:
        out += ['# Matrix made for the user-matrix oracle', '# Entries for the matrix at a scale of 1/2 bit', '']
    header = sep.join(list(letters) + [star])
    out.append((sep + header) if layout['indent'] else header)
    rows = list(letters) + [star]
    for i, a in enumerate(rows):
        cells = []
        for j, b in enumerate(rows):
            if a == star or b == star:
                v = -4.0
            elif j <= i or not layout['junk']:
                v = d.get((a, b), 0.0)
            else:   # upper triangle: junk that must be ignored
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


def make_files(kind):
    """Returns (files dict, options dict)."""
    files, opts = {}, {}
    n = 0

    def add_matrix():
        nonlocal n
        d, letters = pick_matrix()
        name = f'm{n}.mat'
        n += 1
        files[name] = render(d, letters, random_layout())
        return name

    if kind in ('matrix', 'both'):
        opts['matrix'] = add_matrix()
    if kind == 'series':
        k = rnd.randint(2, 4)
        cuts = sorted(rnd.sample(range(5, 96), k - 1))
        bounds = [0] + cuts + [100]
        lines = []
        if rnd.random() < 0.5:
            lines.append('# series of matrices')
        lines.append('CLUSTAL_SERIES')
        style = rnd.choice(['tile', 'tile', 'gaps', 'overlap'])
        for i in range(k):
            lo, hi = bounds[i], bounds[i + 1]
            if style == 'gaps' and 0 < lo < 100 and rnd.random() < 0.5:
                lo = min(lo + rnd.randint(1, 5), hi - 1)
            if style == 'overlap' and i > 0:
                lo = max(0, lo - rnd.randint(1, 10))
            lines.append(f'MATRIX {lo} {hi} {add_matrix()}')
            if rnd.random() < 0.2:
                lines.append('# comment between entries')
        files['series.txt'] = '\n'.join(lines) + '\n'
        opts['matrix'] = 'series.txt'
    if kind in ('pw', 'both'):
        opts['pwmatrix'] = add_matrix()
    return files, opts


def run(seqs, o, files):
    with tempfile.TemporaryDirectory() as d:
        for name, text in files.items():
            with open(os.path.join(d, name), 'w', newline='') as f:
                f.write(text)
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
        out = os.path.join(d, 'o.fa')
        args = ['clustalw', f'-INFILE={fa}', '-TYPE=PROTEIN', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}']
        for k, v in o.items():
            args.append(f'-{g.FLAGS[k]}' + ('' if v is None else f'={v}'))
        p = subprocess.run(args, capture_output=True, text=True, cwd=d)
        if not os.path.exists(out):
            return None, None, p.stdout + p.stderr
        rows, name = {}, None
        for line in open(out):
            line = line.rstrip('\n')
            if line.startswith('>'): name = line[1:].strip(); rows[name] = ''
            elif name: rows[name] += line.strip()
        scores = re.findall(r'Sequences \((\d+):(\d+)\) Aligned\. Score:\s+(-?\d+)', p.stdout)
        return [rows[f's{i}'] for i in range(len(seqs))], scores, p.stdout


def main(outdir, nsets):
    os.makedirs(outdir, exist_ok=True)
    kinds = [(0.5, 70), (0.3, 130), (0.15, 220), (0.05, 320)]
    done = 0
    with open(os.path.join(outdir, 'clustalw_user_matrix_oracle.tsv'), 'w') as f:
        f.write('id\tseqs\trows\toptions\tfiles\tpair_scores\n')
        tries = 0
        while done < nsets:
            r = rnd.random(); acc = 0
            for p, ml in kinds:
                acc += p
                if r < acc: break
            seqs = g.make_set(ml, 2, rnd.choice([4, 6, 8, 8, 12]))
            if rnd.random() < 0.08: seqs = g.ambiguous(seqs)
            kind = rnd.choice(['matrix', 'matrix', 'series', 'series', 'pw', 'both', 'series'])
            files, uo = make_files(kind)
            o = g.options()
            o.pop('matrix', None); o.pop('pwmatrix', None)
            o = {**uo, **o}
            rows, scores, log = run(seqs, o, files)
            tries += 1
            if rows is None:
                print('clustalw failed on try', tries, kind, log[-300:]); continue
            assert all(len(x) == len(rows[0]) for x in rows), (tries, o)
            ostr = ','.join(k if v is None else f'{k}={v}' for k, v in o.items())
            ps = ';'.join(f'{a}:{b}:{s}' for a, b, s in scores)
            f.write(f"{done}\t{';'.join(seqs)}\t{';'.join(rows)}\t{ostr}\t{json.dumps(files)}\t{ps}\n")
            done += 1
    print(done, 'sets in', tries, 'tries')


if __name__ == '__main__':
    if len(sys.argv) > 3:
        rnd.seed(int(sys.argv[3])); g.rnd.seed(int(sys.argv[3]))
    main(sys.argv[1] if len(sys.argv) > 1 else '.', int(sys.argv[2]) if len(sys.argv) > 2 else 240)
