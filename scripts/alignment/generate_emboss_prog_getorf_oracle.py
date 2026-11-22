#!/usr/bin/env python3
"""B09 WP61 (F280) oracle: EMBOSS 6.6.0 `getorf -auto -sequence in.dat [-table T] [-minsize N] [-maxsize N] [-find 0..6] [-methionine y/n] [-circular] [-noreverse] [-flanking N] [range / case qualifiers] -outseq stdout [-osformat F]`.
Inputs: one to four nucleotide sequences (planted ORFs, ambiguity codes, lower / mixed case, U, X, N runs, empty), some protein or damaged inputs, input formats (fasta, embl, genbank, gcg, nbrf, raw ...).
Usage: python3 scripts/alignment/generate_emboss_prog_getorf_oracle.py OUT.jsonl.gz [n=600] [seed]      (needs EMBOSS 6.6.0)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import emboss_prog_oracle_lib as L  # noqa: E402
import generate_seqret_oracle as S  # noqa: E402

TABLES = ['0', '1', '2', '3', '4', '5', '6', '9', '10', '11', '12', '13', '14', '15', '16', '21', '22', '23']
IN_FORMATS = ['fasta'] * 8 + ['embl', 'genbank', 'gcg', 'nbrf', 'raw', 'ig', 'ncbi', 'staden', 'codata', 'fastq']
OUT_FORMATS = ['fasta', 'embl', 'genbank', 'gcg', 'nbrf', 'raw', 'msf', 'clustal', 'phylip', 'ig', 'swiss', 'fastq', 'ncbi', 'text', 'staden', 'pir', 'codata', 'strider', 'mega', 'nexus', 'gff3', 'sam']


def gen(rnd, i, d):
    kind = 'n'
    r = rnd.random()
    opts = dict(orfs=True, iupac=rnd.choice([0, 0, 0, 0.02, 0.1]), case=rnd.choice(['upper', 'upper', 'upper', 'lower', 'mixed']), rna=rnd.random() < 0.1, x=rnd.random() < 0.1)
    seqs = L.rand_seqs(rnd, 'n', lens=(20, rnd.choice([100, 300, 600])), **opts)
    if r < 0.05:
        seqs = L.rand_seqs(rnd, 'p', nseq=1)
    text = L.fasta(rnd, seqs)
    fmt = rnd.choice(IN_FORMATS)
    text = L.convert(text, fmt, d)
    if rnd.random() < 0.06:
        text = S.damaged(rnd, text)
    args = ['-sequence', 'in.dat']
    if rnd.random() < 0.5: args += ['-table', rnd.choice(TABLES)]
    if rnd.random() < 0.8: args += ['-minsize', str(rnd.choice([0, 0, 1, 3, 3, 6, 6, 9, 12, 20, 30, 45, 60, 100]))]
    if rnd.random() < 0.3: args += ['-maxsize', str(rnd.choice([0, 3, 9, 30, 60, 150, 600, 1000000]))]
    if rnd.random() < 0.85: args += ['-find', str(rnd.randint(0, 6))]
    if rnd.random() < 0.3: args += ['-methionine', rnd.choice(['y', 'n', 'no', 'yes'])]
    if rnd.random() < 0.3: args += ['-circular']
    if rnd.random() < 0.25: args += ['-noreverse']
    if rnd.random() < 0.4: args += ['-flanking', str(rnd.choice([0, 1, 3, 10, 20, 50, 100, 150]))]
    if rnd.random() < 0.2: args += ['-sbegin1', str(rnd.randint(-20, 60))]
    if rnd.random() < 0.2: args += ['-send1', str(rnd.randint(-60, 200))]
    if rnd.random() < 0.08: args += ['-sreverse1']
    if rnd.random() < 0.05: args += [rnd.choice(['-slower1', '-supper1'])]
    if rnd.random() < 0.15 and fmt != 'fasta': args += ['-sformat1', fmt]
    if rnd.random() < 0.07: args += ['-snucleotide1']
    outfmt = 'fasta' if rnd.random() < 0.75 else rnd.choice(OUT_FORMATS)
    if outfmt != 'fasta': args += ['-osformat', outfmt]
    if rnd.random() < 0.04:
        args += rnd.choice([['-table', '7'], ['-find', '9'], ['-minsize', 'abc'], ['-table', 'x'], ['-flanking', '1e2']])
    args += ['-outseq', 'stdout']
    return args, {'in.dat': text}


if __name__ == '__main__':
    L.main('getorf', gen, 600)
