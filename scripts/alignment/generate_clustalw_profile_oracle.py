#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A2-14, B09 F103) oracle: ClustalW 2.1 profile alignment commands `-PROFILE1=a -PROFILE2=b -PROFILE` (two alignments merged by profile-profile alignment,
`Clustal::profileAlign`) and `-PROFILE1=a -PROFILE2=b -SEQUENCES` (the sequences of file b added one by one to the alignment of file a, `Clustal::sequencesAlignToProfile`).

Seqeron API: SequenceAligner.AlignClustalWProfiles(profile1, profile2, options) / SequenceAligner.AlignClustalWSequencesToProfile(profile, sequences, options).

Reference (run here): the installed `clustalw` 2.1 (Debian 2.1+lgpl-7) `-TYPE=DNA|PROTEIN -PROFILE1=a.fa -PROFILE2=b.fa -PROFILE|-SEQUENCES -OUTPUT=FASTA -OUTORDER=INPUT` with the option variants
-GAPOPEN -GAPEXT -MATRIX -NOWEIGHTS -NEGATIVE -TRANSWEIGHT -NOPGAP -NOHGAP -GAPDIST -MAXDIV (a profile is an alignment produced by `clustalw -ALIGN` itself, with gaps; profile 2 of `-SEQUENCES` is a
list of unaligned sequences, or in 25 % of the sets an alignment) and -PWMATRIX / -PWGAPOPEN / -PWGAPEXT / -QUICKTREE for the pairwise stage of -SEQUENCES. Names are unique (s000, s001, ...).

Output TSV (tests/.../TestData/Alignment/clustalw_profile_oracle.tsv):  id  mode(profile|sequences)  type(dna|protein)  profile1(;)  profile2(;)  rows(;)  options(k=v,k=v)
Usage: python3 scripts/alignment/generate_clustalw_profile_oracle.py OUTFILE [nSets] [seed]      (needs `clustalw` on PATH)
"""
import importlib.util, os, random, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, name + '.py'))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


dna = load('generate_clustalw_msa_full_oracle')
prot = load('generate_clustalw_protein_msa_oracle')

FLAGS = {'gapopen': 'GAPOPEN', 'gapext': 'GAPEXT', 'maxdiv': 'MAXDIV', 'noweights': 'NOWEIGHTS', 'negative': 'NEGATIVE', 'quick': 'QUICKTREE', 'matrix': 'MATRIX', 'pwmatrix': 'PWMATRIX',
         'nopgap': 'NOPGAP', 'nohgap': 'NOHGAP', 'gapdist': 'GAPDIST', 'transweight': 'TRANSWEIGHT', 'pwgapopen': 'PWGAPOPEN', 'pwgapext': 'PWGAPEXT', 'iteration': 'ITERATION'}


def options(rnd, typ, mode):
    o = {}
    if rnd.random() < 0.3: return o
    if rnd.random() < 0.3: o['gapopen'] = rnd.choice([5, 10, 15, 20, 25]) if typ == 'dna' else rnd.choice([2, 5, 8, 10, 12.5, 15, 20])
    if rnd.random() < 0.3: o['gapext'] = rnd.choice([0.5, 1, 3, 6.66, 8]) if typ == 'dna' else rnd.choice([0.05, 0.1, 0.2, 0.5, 1, 2])
    if rnd.random() < 0.2: o['maxdiv'] = rnd.choice([0, 10, 30, 50, 100])
    if rnd.random() < 0.15: o['noweights'] = None
    if rnd.random() < 0.15: o['negative'] = None
    if mode == 'sequences':
        if rnd.random() < 0.15: o['quick'] = None
        if rnd.random() < 0.15: o['pwgapopen'] = rnd.choice([5, 10, 15, 20])
        if rnd.random() < 0.15: o['pwgapext'] = rnd.choice([0.5, 1, 3.33, 6.66])
    if typ == 'dna':
        if rnd.random() < 0.15: o['transweight'] = rnd.choice([0, 0.25, 0.5, 1])
    else:
        if rnd.random() < 0.3: o['matrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
        if mode == 'sequences' and rnd.random() < 0.15: o['pwmatrix'] = rnd.choice(['BLOSUM', 'PAM', 'ID', 'GONNET'])
        if rnd.random() < 0.15: o['nopgap'] = None
        if rnd.random() < 0.15: o['nohgap'] = None
        if rnd.random() < 0.2: o['gapdist'] = rnd.choice([0, 1, 2, 4, 6, 8, 12])
    return o


def write(path, names_rows):
    with open(path, 'w') as f:
        for n, r in names_rows: f.write(f'>{n}\n{r}\n')


def read(path):
    rows = {}
    name = None
    for line in open(path):
        line = line.rstrip('\n')
        if line.startswith('>'): name = line[1:].strip(); rows[name] = ''
        elif name: rows[name] += line.strip()
    return rows


def align(d, seqs, names, typ, tag):
    fa = os.path.join(d, tag + '.fa'); out = os.path.join(d, tag + '.aln.fa')
    write(fa, list(zip(names, seqs)))
    if len(seqs) == 1:
        return [seqs[0]]
    subprocess.run(['clustalw', f'-INFILE={fa}', f'-TYPE={typ.upper()}', '-ALIGN', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}'], capture_output=True, text=True, cwd=d)
    r = read(out)
    return [r[n] for n in names]


def main(outfile, nsets, seed):
    rnd = random.Random(seed)
    dna.rnd.seed(seed + 1)
    prot.rnd.seed(seed + 2)
    with open(outfile, 'w') as f:
        f.write('id\tmode\ttype\tprofile1\tprofile2\trows\toptions\n')
        for i in range(nsets):
            typ = 'dna' if i % 2 == 0 else 'protein'
            mode = 'profile' if (i // 2) % 2 == 0 else 'sequences'
            maxlen = rnd.choice([40, 60, 90, 140])
            kmax = rnd.choice([4, 6, 8, 10, 14])
            seqs = dna.make_set(maxlen, "ACGT", 2, kmax) if typ == 'dna' else prot.make_set(maxlen, 2, kmax)
            if len(seqs) < 2: seqs = seqs + seqs[:1]
            cut = rnd.randint(1, len(seqs) - 1)
            names = [f's{j:03d}' for j in range(len(seqs))]
            with tempfile.TemporaryDirectory() as d:
                p1 = align(d, seqs[:cut], names[:cut], typ, 'a')
                if mode == 'profile':
                    p2 = align(d, seqs[cut:], names[cut:], typ, 'b')
                else:
                    # -SEQUENCES: profile 2 is normally a list of sequences, but the program also accepts an aligned file (the gaps stay)
                    p2 = align(d, seqs[cut:], names[cut:], typ, 'b') if rnd.random() < 0.25 else seqs[cut:]
                o = options(rnd, typ, mode)
                write(os.path.join(d, 'a.fa'), list(zip(names[:cut], p1)))
                write(os.path.join(d, 'b.fa'), list(zip(names[cut:], p2)))
                out = os.path.join(d, 'o.fa')
                args = ['clustalw', '-PROFILE1=a.fa', '-PROFILE2=b.fa', '-PROFILE' if mode == 'profile' else '-SEQUENCES', f'-TYPE={typ.upper()}', '-OUTPUT=FASTA', '-OUTORDER=INPUT', f'-OUTFILE={out}']
                for k, v in o.items(): args.append(f'-{FLAGS[k]}' + ('' if v is None else f'={v}'))
                r = subprocess.run(args, capture_output=True, text=True, cwd=d)
                if not os.path.exists(out):
                    print('no output for set', i, r.stdout[-200:], file=sys.stderr); continue
                rows = read(out)
                try:
                    res = [rows[n] for n in names]
                except KeyError:
                    print('missing rows for set', i, file=sys.stderr); continue
            ostr = ','.join(k if v is None else f'{k}={v}' for k, v in o.items())
            f.write(f"{i}\t{mode}\t{typ}\t{';'.join(p1)}\t{';'.join(p2)}\t{';'.join(res)}\t{ostr}\n")
    print(nsets, 'sets')


if __name__ == '__main__':
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) > 2 else 1200, int(sys.argv[3]) if len(sys.argv) > 3 else 20261120)
