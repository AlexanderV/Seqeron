#!/usr/bin/env python3
"""ALIGN-MULTI-001 (audit A3-14g, B09 F143): MUSCLE 3.8.1551 `-maxiters 2` versus `-maxiters 8` on 12 random protein families (the families of
SequenceAligner_MuscleOptions_Tests.Iterations3AndAbove_ChangeSomeAlignmentsOfIteration2_AndNucleotidesAreSupported: C# `Random(5)`, 5 sequences of 30-49 residues, four mutated copies of the first).

Usage: python3 scripts/alignment/generate_muscle3_iterations_oracle.py INPUT.tsv > muscle3_iterations_oracle.tsv      (INPUT.tsv: trial <TAB> seqs joined by '|'; needs `muscle3` 3.8.1551 on PATH)
Output: trial <TAB> seqs <TAB> rows of -maxiters 2 (';', input order) <TAB> rows of -maxiters 8 (';', input order)
"""
import subprocess, sys, tempfile, os

def run(seqs, iters, d):
    with open(d + '/in.fa', 'w') as f:
        for i, s in enumerate(seqs): f.write(f'>s{i}\n{s}\n')
    r = subprocess.run(['muscle3', '-in', d + '/in.fa', '-out', d + '/out.fa', '-seqtype', 'protein', '-maxiters', str(iters), '-quiet'], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    rows = {}
    name = None
    for line in open(d + '/out.fa'):
        line = line.rstrip('\n')
        if line.startswith('>'): name = line[1:]; rows[name] = ''
        else: rows[name] += line
    return ';'.join(rows[f's{i}'] for i in range(len(seqs)))

print('trial\tseqs\trows_maxiters2\trows_maxiters8')
with tempfile.TemporaryDirectory() as d:
    for line in list(open(sys.argv[1]))[1:]:
        trial, seqs = line.rstrip('\n').split('\t')
        s = seqs.split('|')
        print('\t'.join([trial, seqs, run(s, 2, d), run(s, 8, d)]))
