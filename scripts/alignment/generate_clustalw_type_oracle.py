#!/usr/bin/env python3
"""B09 audit A2-14: oracle for the ClustalW 2.1 residue-type guess (`-TYPE` omitted = auto): the installed `clustalw -INFILE=x.fa -ALIGN` prints "Sequence n: name  len bp|aa" with the unit
of the type it decided (FileReader::readFile: the flag is set from the SECOND sequence only, Sequence::checkDNAFlag: DNA when >= 85 % of its non-gap characters are A C G T U N).

    python3 scripts/alignment/generate_clustalw_type_oracle.py [n=500] > tests/.../TestData/Alignment/clustalw_type_oracle.tsv

Columns: sequences(';'-joined, as written to the FASTA file)  expected(dna|protein: bp / aa of the printed lines; they are equal for every sequence)
"""
import random, re, subprocess, sys, tempfile

rnd = random.Random(20261119)
PROT = list('ACDEFGHIKLMNPQRSTVWY')
n = int(sys.argv[1]) if len(sys.argv) > 1 else 500


def dna_like(frac_other):
    L = rnd.randint(8, 60)
    s = []
    for _ in range(L):
        r = rnd.random()
        s.append(rnd.choice('ACGT') if r > frac_other else rnd.choice('RYKMSWBDHVNXN'))
    if rnd.random() < 0.3:
        for _ in range(rnd.randint(1, 4)):
            s[rnd.randrange(L)] = '-'
    seq = ''.join(s)
    return seq.lower() if rnd.random() < 0.25 else seq


def prot_like():
    L = rnd.randint(8, 60)
    seq = ''.join(rnd.choice(PROT) for _ in range(L))
    if rnd.random() < 0.3:
        seq = ''.join(rnd.choice('ACGT') if rnd.random() < 0.5 else c for c in seq)
    return seq.lower() if rnd.random() < 0.2 else seq


print('sequences\texpected')
with tempfile.TemporaryDirectory() as d:
    for _ in range(n):
        k = rnd.choice([1, 2, 2, 3, 3, 4, 5, 6])
        seqs = []
        for i in range(k):
            kind = rnd.random()
            if kind < 0.45: seqs.append(dna_like(rnd.choice([0, 0, 0.05, 0.1, 0.15, 0.2])))
            elif kind < 0.55: seqs.append(dna_like(rnd.choice([0.12, 0.15, 0.16, 0.18])))      # around the 85 % threshold
            else: seqs.append(prot_like())
        if all(set(s.replace('-', '')) == set() for s in seqs): continue
        open(d + '/x.fa', 'w').write(''.join(f'>s{i}\n{s}\n' for i, s in enumerate(seqs)))
        p = subprocess.run(['clustalw', '-INFILE=' + d + '/x.fa', '-ALIGN', '-OUTFILE=' + d + '/o.aln', '-NEWTREE=' + d + '/o.dnd'], capture_output=True, text=True)
        units = set(re.findall(r'^Sequence \d+: \S+\s+\d+ (bp|aa)\s*$', p.stdout, re.M))
        if len(units) != 1:
            continue
        print(';'.join(seqs) + '\t' + ('dna' if units == {'bp'} else 'protein'))
