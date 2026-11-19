#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-14, B09 F102) oracle: ClustalW 2.1 `-ITERATION=TREE | ALIGNMENT` / `-NUMITER` (src/multipleAlign/Iteration.cpp, src/alignment/ObjectiveScore.cpp).

Seqeron API: SequenceAligner.MultipleAlignClustalW(seqs, new ClustalWOptions / ClustalWProteinOptions { Iteration = ..., NumIterations = ... }).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7) `-TYPE=DNA|PROTEIN -ALIGN -ITERATION=tree|alignment -NUMITER=n -OUTPUT=FASTA -OUTORDER=INPUT` (the program accepts only the
LOWER-CASE words; "ALIGNMENT" prints "Unknown option for iteration. Setting to NONE") with the option variants of the plain multiple-alignment oracles. The binary computes its own pairwise
distances, guide tree, weights, profile alignments and the remove-first iteration; the C# port must reproduce the FASTA rows exactly (input order) AND the "FINAL score: n" lines (the
ObjectiveScore of the best re-alignment of every Iteration::removeFirstIterate call, in printed order: one for -ITERATION=alignment, one per profile step of more than two sequences for tree).

Output TSV (tests/.../TestData/Alignment/clustalw_iteration_oracle.tsv):  id  type(dna|protein)  seqs(;)  rows(;)  options(k=v,k=v)  finals(,)
Usage: python3 scripts/alignment/generate_clustalw_iteration_oracle.py OUTFILE [nSets] [seed]      (needs `clustalw` on PATH)
Committed file: from `... clustalw_iteration_oracle.tsv 700`.
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

FLAGS = {'gapopen': 'GAPOPEN', 'gapext': 'GAPEXT', 'maxdiv': 'MAXDIV', 'noweights': 'NOWEIGHTS', 'negative': 'NEGATIVE', 'quick': 'QUICKTREE', 'matrix': 'MATRIX',
         'nopgap': 'NOPGAP', 'nohgap': 'NOHGAP', 'novgap': 'NOVGAP', 'gapdist': 'GAPDIST', 'transweight': 'TRANSWEIGHT', 'iteration': 'ITERATION', 'numiter': 'NUMITER'}


def options(rnd, typ):
    o = {}
    r = rnd.random()
    o['iteration'] = 'alignment' if r < 0.5 else 'tree'
    if rnd.random() < 0.4: o['numiter'] = rnd.choice([1, 2, 3, 5, 8])
    if rnd.random() < 0.25: o['gapopen'] = rnd.choice([5, 10, 15, 20, 25]) if typ == 'dna' else rnd.choice([2, 5, 8, 10, 12.5, 15, 20])
    if rnd.random() < 0.25: o['gapext'] = rnd.choice([0.5, 1, 3, 6.66, 8]) if typ == 'dna' else rnd.choice([0.05, 0.1, 0.2, 0.5, 1, 2])
    if rnd.random() < 0.2: o['maxdiv'] = rnd.choice([0, 10, 30, 50, 70, 100])
    if rnd.random() < 0.15: o['noweights'] = None
    if rnd.random() < 0.15: o['negative'] = None
    if rnd.random() < 0.15: o['quick'] = None
    if typ == 'dna':
        if rnd.random() < 0.15: o['transweight'] = rnd.choice([0, 0.25, 0.5, 1])
    else:
        if rnd.random() < 0.3: o['matrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
        if rnd.random() < 0.15: o['nopgap'] = None
        if rnd.random() < 0.15: o['nohgap'] = None
        if rnd.random() < 0.2: o['gapdist'] = rnd.choice([0, 1, 2, 4, 6, 8, 12])
    return o


def run(seqs, typ, o):
    with tempfile.TemporaryDirectory() as d:
        fa = os.path.join(d, 'x.fa')
        with open(fa, 'w') as f:
            for i, s in enumerate(seqs): f.write(f'>s{i:03d}\n{s}\n')
        out = os.path.join(d, 'o.fa')
        args = ['clustalw', f'-INFILE={fa}', f'-TYPE={typ.upper()}', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}']
        for k, v in o.items():
            args.append(f'-{FLAGS[k]}' + ('' if v is None else f'={v}'))
        r = subprocess.run(args, capture_output=True, text=True, cwd=d)
        rows = {}
        name = None
        for line in open(out):
            line = line.rstrip('\n')
            if line.startswith('>'): name = line[1:].strip(); rows[name] = ''
            elif name: rows[name] += line.strip()
        finals = re.findall(r'FINAL score: (-?\d+)', r.stdout)
    return [rows[f's{i:03d}'] for i in range(len(seqs))], finals


def main(outfile, nsets, seed):
    rnd = random.Random(seed)
    dna.rnd.seed(seed + 1)
    prot.rnd.seed(seed + 2)
    changed = 0
    with open(outfile, 'w') as f:
        f.write('id\ttype\tseqs\trows\toptions\tfinals\n')
        for i in range(nsets):
            typ = 'dna' if i % 2 == 0 else 'protein'
            maxlen = rnd.choice([40, 60, 90, 140])
            kmax = rnd.choice([4, 6, 8, 10, 14])
            if typ == 'dna':
                seqs = dna.make_set(maxlen, "ACGT", 3, kmax)
            else:
                seqs = prot.make_set(maxlen, 3, kmax)
            o = options(rnd, typ)
            rows, finals = run(seqs, typ, o)
            base, _ = run(seqs, typ, {k: v for k, v in o.items() if k not in ('iteration', 'numiter')})
            changed += rows != base
            ostr = ','.join(k if v is None else f'{k}={v}' for k, v in o.items())
            f.write(f"{i}\t{typ}\t{';'.join(seqs)}\t{';'.join(rows)}\t{ostr}\t{','.join(finals)}\n")
    print(nsets, 'sets;', changed, 'with rows different from the non-iterated alignment')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 700, int(sys.argv[3]) if len(sys.argv) > 3 else 20261119)
