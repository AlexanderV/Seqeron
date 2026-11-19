#!/usr/bin/env python3
"""Builds tests/.../TestData/Alignment/kalign_oracle.tsv (B09 audit A2-20) with the installed Kalign 3.3.5 (`kalign`, Debian 1:3.3.5-1).

Every row = one run `kalign -i in.fa -o out.fa -q -n THREADS [--type T] [--gpo X] [--gpe Y] [--tgpe Z]` on a deterministic FASTA
(one line per sequence). The sets are generated here (python `random.Random`, fixed seeds) so the file can be rebuilt bit for bit:

  * DNA / RNA / protein families: a random ancestor, 2 - 12 descendants with substitutions and indels, a few IUPAC / ambiguity codes,
    lower-case letters, '-' / '.' characters (the program removes them), zero-length sequences (the program drops them);
  * every `--type` (none, dna, rna, internal, protein, divergent) including the type / alphabet mismatches the program refuses (rc != 0);
  * every gap-penalty combination (--gpo, --gpe, --tgpe alone and together; negative = the program's default);
  * -n 1 (serial dynamic programs) and -n 4 / -n 2 (the OpenMP build: `aln_runner`, the swap of the shorter sequence, ...);
  * long sets (> 500 columns: the divide-and-conquer branch of `aln_runner`) and sets of >= 100 / 1000 sequences (bisecting k-means guide tree).

Columns: id, kind, options, threads, names, seqs, rc, rows ('|' separated).  Run: python3 scripts/alignment/generate_kalign_oracle.py
"""
import os, random, subprocess, sys, tempfile

root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')
out_path = os.path.join(root, 'tests', 'Seqeron', 'Seqeron.Genomics.Tests', 'TestData', 'Alignment', 'kalign_oracle.tsv')
binary = 'kalign'

DNA = 'ACGT'
PROT = 'ACDEFGHIKLMNPQRSTVWY'


def mutate(rng, base, alpha, sub, indel):
    out = []
    for ch in base:
        r = rng.random()
        if r < indel / 2:
            continue
        if r < indel:
            out.append(ch)
            out.append(rng.choice(alpha))
            continue
        if rng.random() < sub:
            out.append(rng.choice(alpha))
        else:
            out.append(ch)
    return ''.join(out) or base[:1]


def family(rng, n, length, alpha, sub, indel):
    base = ''.join(rng.choice(alpha) for _ in range(length))
    return [mutate(rng, base, alpha, sub, indel) for _ in range(n)]


def decorate(rng, seqs, kind):
    """Lower case runs, ambiguity codes, gap characters - everything the program accepts and ignores / keeps."""
    res = []
    for s in seqs:
        s = list(s)
        if rng.random() < 0.25 and len(s) > 4:
            a = rng.randrange(len(s)); b = min(len(s), a + rng.randrange(1, 8))
            for i in range(a, b):
                s[i] = s[i].lower()
        if rng.random() < 0.25:
            amb = {'dna': 'NRYKM', 'rna': 'NRYKM', 'prot': 'XBZ'}[kind]
            for _ in range(rng.randrange(1, 3)):
                s[rng.randrange(len(s))] = rng.choice(amb)
        if rng.random() < 0.15:
            for _ in range(rng.randrange(1, 3)):
                s.insert(rng.randrange(len(s) + 1), rng.choice('-.'))
        res.append(''.join(s))
    return res


def to_rna(seqs):
    return [s.replace('T', 'U').replace('t', 'u') for s in seqs]


cases = []  # (id, kind, options list, threads, names, seqs)


def add(kind, options, threads, seqs):
    cid = 'K%04d' % (len(cases) + 1)
    cases.append((cid, kind, options, threads, ['s%d' % (i + 1) for i in range(len(seqs))], seqs))


rng = random.Random(20261119)
types_nuc = [None, 'dna', 'rna', 'internal']
types_prot = [None, 'protein', 'divergent']
gap_sets = [[], ['--gpo', '4'], ['--gpe', '2.5'], ['--tgpe', '3'], ['--gpo', '6', '--gpe', '1', '--tgpe', '0.5'],
            ['--gpo', '12', '--gpe', '4'], ['--gpe', '0'], ['--gpo', '0'], ['--gpo', '-1', '--gpe', '-1', '--tgpe', '-1'],
            ['--gpo', '30', '--gpe', '10', '--tgpe', '5'], ['--gpe', '1', '--tgpe', '-5']]


def opts(t, g):
    return (['--type', t] if t else []) + g


# 1. nucleotide families, every type and gap set, -n 4 and -n 1
for i in range(150):
    n = rng.choice([2, 3, 3, 4, 5, 6, 8, 10, 12])
    length = rng.choice([12, 20, 30, 45, 60, 90, 130, 200])
    seqs = family(rng, n, length, DNA, rng.choice([0.03, 0.1, 0.25]), rng.choice([0.0, 0.03, 0.08]))
    kind = 'dna'
    if i % 3 == 2:
        seqs = to_rna(seqs); kind = 'rna'
    seqs = decorate(rng, seqs, kind)
    t = types_nuc[i % 4]
    g = gap_sets[i % len(gap_sets)]
    add(kind, opts(t, g), rng.choice([4, 4, 1, 2]), seqs)
# 2. protein families
for i in range(120):
    n = rng.choice([2, 3, 4, 5, 6, 8, 10, 12])
    length = rng.choice([15, 25, 40, 60, 100, 150, 250])
    seqs = decorate(rng, family(rng, n, length, PROT, rng.choice([0.05, 0.2, 0.4]), rng.choice([0.0, 0.03, 0.1])), 'prot')
    t = types_prot[i % 3]
    g = gap_sets[(i * 7) % len(gap_sets)]
    add('prot', opts(t, g), rng.choice([4, 4, 1, 2]), seqs)
# 3. the same set under -n 1 / 2 / 4 (the serial and the OpenMP variants of do_align)
for i in range(30):
    n = rng.choice([3, 5, 7])
    seqs = family(rng, n, rng.choice([40, 80, 160]), DNA if i % 2 == 0 else PROT, 0.15, 0.06)
    kind = 'dna' if i % 2 == 0 else 'prot'
    for th in (1, 2, 4):
        add(kind, [], th, seqs)
# 4. type / alphabet mismatches (the program stops with rc != 0) and corner inputs
for t in ('protein', 'divergent'):
    add('dna', ['--type', t], 4, family(rng, 4, 40, DNA, 0.1, 0.02))
for t in ('dna', 'rna', 'internal'):
    add('prot', ['--type', t], 4, family(rng, 4, 40, PROT, 0.2, 0.02))
add('dna', [], 4, ['ACGTACGT', '----', 'ACGTCGT'])          # zero-length sequence dropped
add('dna', [], 4, ['ACGTACGT', '----', '....'])             # only one sequence left -> error
add('dna', [], 4, ['ACGTACGTAC', 'ACGTACGTAC'])             # identical
add('dna', [], 4, ['A', 'A', 'C'])                          # one-residue sequences
add('dna', [], 4, ['A', 'ACGTACGTACGT'])
add('prot', [], 4, ['M', 'MKV', 'MKVL'])
add('dna', [], 1, ['AC', 'GT', 'AT', 'CG'])
add('prot', [], 4, ['ACDEFGHIK', 'ACDEFGHIK', 'ACDEFGHIK', 'ACDEFGHIK'])
add('dna', [], 4, ['AAAAAAAAAAAAAA', 'AAAAAAAAAAAA', 'AAAAAAAAAAAAAAAAAA', 'AAAAAA'])
add('dna', [], 4, ['ACGTNNNNNNACGT', 'ACGTRYKMSWACGT', 'ACGTBDHVACGT', 'ACGTACGT'])
# 5. DNA / protein detection borderline: short mixed sequences made of the letters ACGTN only, and protein with letters that are also nucleotide codes
for i in range(8):
    add('prot', [], 4, [''.join(rng.choice('ACGT') for _ in range(rng.randrange(8, 20))) for _ in range(3)])
for i in range(6):
    add('prot', [], 4, [''.join(rng.choice('ACGTNRYMKACDEFGHIKLMNPQRSTVWY') for _ in range(rng.randrange(10, 30))) for _ in range(4)])
# 6. long sequences (the divide-and-conquer branch of aln_runner: more than 500 columns)
for i in range(14):
    n = rng.choice([3, 4, 6])
    if i % 2 == 0:
        seqs = family(rng, n, rng.choice([520, 640, 900, 1200]), DNA, 0.08, 0.02)
        kind = 'dna'
    else:
        seqs = family(rng, n, rng.choice([510, 600, 800]), PROT, 0.2, 0.02)
        kind = 'prot'
    for th in (1, 4):
        add(kind, [], th, seqs)
# 7. many sequences: bisecting k-means guide tree (>= 100 sequences) and the parallel variant (>= 1000)
for i, (n, length, alpha, kind) in enumerate([(100, 20, DNA, 'dna'), (101, 30, DNA, 'dna'), (120, 40, DNA, 'dna'), (150, 25, PROT, 'prot'),
                                              (130, 60, PROT, 'prot'), (200, 18, DNA, 'dna'), (110, 12, PROT, 'prot'), (99, 30, DNA, 'dna')]):
    seqs = family(rng, n, length, alpha, 0.2, 0.05)
    seqs = [s[:max(3, len(s) - rng.randrange(0, 6))] for s in seqs]
    for th in ((4, 1) if i < 4 else (4,)):
        add(kind, [], th, seqs)
seqs = family(rng, 1000, 14, DNA, 0.25, 0.05)
add('dna', [], 4, seqs)

# 8. the distance's 255-residue pattern limit (sequences longer than 255), larger UPGMA sets (30 - 90 sequences: profile-profile alignments), extreme gap penalties
for i in range(12):
    n = rng.choice([3, 4, 5])
    L = rng.choice([256, 270, 300, 400, 520])
    seqs = family(rng, n, L, DNA if i % 2 == 0 else PROT, 0.1, 0.03)
    kind = 'dna' if i % 2 == 0 else 'prot'
    # unequal lengths so that the shorter one (the pattern) is cut at 255 residues in the distance
    seqs = [x[:len(x) - rng.randrange(0, 40)] for x in seqs]
    add(kind, [], 4, seqs)
for i in range(24):
    n = rng.choice([30, 45, 60, 90])
    kind = 'dna' if i % 2 == 0 else 'prot'
    seqs = family(rng, n, rng.choice([25, 40, 70]), DNA if kind == 'dna' else PROT, rng.choice([0.05, 0.2]), rng.choice([0.02, 0.06]))
    add(kind, [], rng.choice([1, 4]), seqs)
extreme = [['--gpo', '0', '--gpe', '0', '--tgpe', '0'], ['--gpo', '100', '--gpe', '50', '--tgpe', '50'], ['--gpo', '1000', '--gpe', '0.001', '--tgpe', '2'],
           ['--gpo', '0.5', '--gpe', '0.5', '--tgpe', '0.5'], ['--gpo', '3', '--gpe', '3', '--tgpe', '-1'], ['--gpe', '7.77', '--tgpe', '9.99']]
for i in range(60):
    kind = ['dna', 'prot', 'rna'][i % 3]
    alpha = DNA if kind != 'prot' else PROT
    seqs = family(rng, rng.choice([3, 4, 6, 8]), rng.choice([30, 60, 120]), alpha, rng.choice([0.05, 0.2, 0.35]), rng.choice([0.02, 0.08]))
    if kind == 'rna':
        seqs = to_rna(seqs)
    t = {'dna': types_nuc, 'rna': types_nuc, 'prot': types_prot}[kind][i % len(types_nuc if kind != 'prot' else types_prot)]
    add(kind, opts(t, extreme[i % len(extreme)]), rng.choice([1, 4]), seqs)

rows = []
fail = 0
with tempfile.TemporaryDirectory() as d:
    for cid, kind, options, threads, names, seqs in cases:
        fin = os.path.join(d, 'in.fa'); fout = os.path.join(d, 'out.fa')
        if os.path.exists(fout):
            os.remove(fout)
        with open(fin, 'w', newline='\n') as f:
            for nme, s in zip(names, seqs):
                f.write('>' + nme + '\n' + s + '\n')
        cmd = [binary, '-i', fin, '-o', fout, '-q', '-n', str(threads)] + options
        p = subprocess.run(cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL, timeout=120)
        out_rows = []
        if p.returncode == 0 and os.path.exists(fout):
            cur = None
            for line in open(fout):
                line = line.rstrip('\n')
                if line.startswith('>'):
                    out_rows.append('')
                else:
                    out_rows[-1] += line
        else:
            fail += 1
        rows.append([cid, kind, ' '.join(options) or '-', str(threads), '|'.join(names), '|'.join(seqs), str(p.returncode), '|'.join(out_rows)])
ver = subprocess.run([binary, '--version'], capture_output=True, text=True).stdout.strip()
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='\n') as fh:
    fh.write('# ' + ver + '; uname -m ' + os.uname().machine + '\n')
    fh.write('id\tkind\toptions\tthreads\tnames\tseqs\trc\trows\n')
    for r in rows:
        fh.write('\t'.join(r) + '\n')
print(len(rows), 'rows;', fail, 'refused / failed;', ver)
