#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/tcoffee_score_oracle.tsv (B09 A2-27d): the T-Coffee 13.41 `score_ascii` evaluation of alignments, run with the installed program
(`t_coffee -infile aln.fa -mode evaluate -output score_ascii -outfile out -n_core 1 -quiet`).

Columns: id, cores, kind, fasta (the aligned rows, FASTA), output (the score_ascii file; "ERROR:" + the first lines of the message when the program refuses).
Escaping of the text columns: backslash -> \\\\, newline -> \\n, tab -> \\t, carriage return -> \\r.
Run: python3 scripts/alignment/generate_tcoffee_score_oracle.py
"""
import os
import random
import re
import shutil
import subprocess
import sys
import tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'tcoffee_score_oracle.tsv')
rng = random.Random(20261323)

NAMES = ['a', 'b', 'c', 'seq1', 'seq2', 'seq3', 'Human', 'mouse', 'rat_gene_2', 'another_long_name_xx', 'q', 'zeta', 'alpha', 'beta', 'sp_P12345_X', 'long_identifier_number_seven', 'S9', 'B2', 'AAA', 'aab']
AA = 'ACDEFGHIKLMNPQRSTVWY'


def esc(s):
    return s.replace('\\', '\\\\').replace('\n', '\\n').replace('\t', '\\t').replace('\r', '\\r')


def mutate(rng, base, alphabet, rate):
    out = []
    for c in base:
        r = rng.random()
        if r < rate / 2:
            continue
        if r < rate:
            out.append(rng.choice(alphabet))
            if rng.random() < 0.2:
                out.append(rng.choice(alphabet))
        else:
            out.append(c)
    return ''.join(out) or base[:3]


def run_tcoffee(args, cwd, env):
    return subprocess.run(['t_coffee'] + args, cwd=cwd, env=env, capture_output=True, text=True, timeout=120)


def read_fasta(text):
    rows = []
    for block in text.split('>')[1:]:
        lines = block.split('\n')
        rows.append((lines[0].split()[0], ''.join(lines[1:])))
    return rows


def random_alignment(rng, seqs):
    """Random gapped rows without all-gap columns."""
    n = len(seqs)
    longest = max(len(s) for s in seqs)
    cols = longest + rng.randrange(0, longest // 3 + 3)
    rows = []
    for s in seqs:
        pos = sorted(rng.sample(range(cols), len(s)))
        row = ['-'] * cols
        for p, c in zip(pos, s):
            row[p] = c
        rows.append(row)
    keep = [c for c in range(cols) if any(r[c] != '-' for r in rows)]
    return [''.join(r[c] for c in keep) for r in rows]


def main():
    tmp = tempfile.mkdtemp(prefix='tcscore')
    env = dict(os.environ, HOME=tmp)
    rows_out = []
    cases = int(sys.argv[1]) if len(sys.argv) > 1 else 280
    for i in range(cases):
        kind = rng.choice(['dna', 'dna', 'protein'])
        alphabet = 'ACGT' if kind == 'dna' else AA
        n = rng.choice([2, 3, 3, 4, 4, 5, 6, 8])
        length = rng.randrange(18, 90 if kind == 'dna' else 70)
        base = ''.join(rng.choice(alphabet) for _ in range(length))
        names = rng.sample(NAMES, n)
        seqs = [mutate(rng, base, alphabet, rng.choice([0.05, 0.15, 0.3])) for _ in range(n)]
        work = os.path.join(tmp, 'w%d' % i)
        os.makedirs(work)
        mode = rng.choice(['tcoffee', 'tcoffee', 'random'])
        if mode == 'tcoffee':
            with open(os.path.join(work, 'in.fa'), 'w') as fh:
                for nm, s in zip(names, seqs):
                    fh.write('>%s\n%s\n' % (nm, s))
            r = run_tcoffee(['-infile', 'in.fa', '-outfile', 'al.fa', '-output', 'fasta_aln', '-n_core', '1', '-quiet'], work, env)
            if not os.path.exists(os.path.join(work, 'al.fa')):
                continue
            aligned = dict(read_fasta(open(os.path.join(work, 'al.fa')).read()))
            rows = [aligned[nm] for nm in names]
        else:
            rows = random_alignment(rng, seqs)
        if rng.random() < 0.25:
            rows = [''.join(c.lower() if c != '-' and rng.random() < 0.3 else c for c in r) for r in rows]
        fasta = ''.join('>%s\n%s\n' % (nm, r) for nm, r in zip(names, rows))
        with open(os.path.join(work, 'x.fa'), 'w') as fh:
            fh.write(fasta)
        r = run_tcoffee(['-infile', 'x.fa', '-mode', 'evaluate', '-output', 'score_ascii', '-outfile', 'out.sa', '-n_core', '1', '-quiet'], work, env)
        path = os.path.join(work, 'out.sa')
        if os.path.exists(path):
            text = open(path).read()
        else:
            text = 'ERROR:' + '\n'.join((r.stdout + r.stderr).strip().split('\n')[-3:])
        rows_out.append((kind, fasta, text))
        shutil.rmtree(work)
    shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', newline='\n') as fh:
        fh.write('id\tcores\tkind\tfasta\toutput\n')
        for i, (kind, fasta, text) in enumerate(rows_out):
            fh.write('T%04d\t1\t%s\t%s\t%s\n' % (i + 1, kind, esc(fasta), esc(text)))
    print(len(rows_out), 'rows;', sum(1 for r in rows_out if r[2].startswith('ERROR')), 'refusals')


main()
