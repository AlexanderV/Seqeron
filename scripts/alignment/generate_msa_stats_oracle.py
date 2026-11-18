#!/usr/bin/env python3
"""ALIGN-STATS-001 (audit A1-9) oracle for MultipleAlignmentStatistics: EMBOSS 6.6.0 infoalign / distmat, ClustalW 2.1 conservation
line and Biopython 1.88 Bio.Phylo.TreeConstruction.DistanceCalculator('identity') on random gapped multiple alignments
(DNA with N/IUPAC, protein with X/B/Z, leading/trailing/internal gaps).

TSV (embedded in tests/.../TestData/Alignment/msa_stats_oracle.tsv): kind  rows(comma-joined)  param  expected
 infoalign : param = reference (0 = consensus, k = k-th row);  expected = per row "seqlen,alignlen,gaps,gaplen,ident,sim,diff,change" joined by ';'
 distmat   : param = "method:gapweight" (0 uncorrected, 1 Jukes-Cantor);  expected = upper triangle of the matrix (x100, 2 decimals) joined by ','
 clustal   : param = dna | protein (ClustalW -type; DNA columns only ever get '*');  expected = ClustalW conservation line ('|' stands for a blank column to keep the TSV field)
 bioid     : param = '';  expected = upper triangle of 100*(1-distance) of DistanceCalculator('identity') joined by ',' (rows ACGT only; Biopython
             scores matches over NON-gap columns but divides by the full alignment length, i.e. distmat gap weight 1)

Usage: python3 scripts/alignment/generate_msa_stats_oracle.py [n_alignments=120] > msa_stats_oracle.tsv
"""
import random, re, subprocess, sys, tempfile

rnd = random.Random(20261118)
DNA = ['A', 'C', 'G', 'T'] * 6 + ['N', 'R', 'Y']
PROT = list('ACDEFGHIKLMNPQRSTVWY') * 2 + ['X', 'B', 'Z']


def make_msa(kind, strict):
    L = rnd.randint(12, 70)
    n = rnd.randint(3, 8)
    alpha = (['A', 'C', 'G', 'T'] if strict else DNA) if kind == 'dna' else (list('ACDEFGHIKLMNPQRSTVWY') if strict else PROT)
    base = [rnd.choice(alpha) for _ in range(L)]
    rows = []
    for _ in range(n):
        rate = rnd.choice([0.0, 0.05, 0.15, 0.3, 0.5])
        row = [c if rnd.random() > rate else rnd.choice(alpha) for c in base]
        for _ in range(rnd.randint(0, 3)):
            s = rnd.randrange(L); e = min(L, s + rnd.randint(1, 5))
            for k in range(s, e): row[k] = '-'
        if rnd.random() < 0.35:
            for k in range(rnd.randint(1, 4)): row[k] = '-'
        if rnd.random() < 0.35:
            for k in range(1, rnd.randint(2, 5)): row[L - k] = '-'
        rows.append(row)
    for k in range(L):  # no all-gap column
        if all(r[k] == '-' for r in rows):
            rows[rnd.randrange(n)][k] = rnd.choice(alpha)
    rows = [''.join(r) for r in rows]
    if any(set(r) == {'-'} for r in rows): return make_msa(kind, strict)
    return rows


def write_fa(path, rows):
    open(path, 'w').write(''.join(f'>s{i + 1}\n{r}\n' for i, r in enumerate(rows)))


def infoalign(d, rows, ref):
    write_fa(d + '/a.fa', rows)
    subprocess.run(['infoalign', '-sequence', d + '/a.fa', '-outfile', d + '/o.txt', '-auto', '-refseq', str(ref), '-nousa', '-noheading', '-nodescription'], check=True, capture_output=True)
    out = []
    for line in open(d + '/o.txt'):
        if line.startswith('#'): continue
        f = line.split()
        # Name SeqLen AlignLen Gaps GapLen Ident Similar Differ %Change Weight
        out.append(','.join(f[1:9]))
    return ';'.join(out)


def distmat(d, rows, method, gapw):
    write_fa(d + '/a.fa', rows)
    subprocess.run(['distmat', '-sequence', d + '/a.fa', '-outfile', d + '/o.txt', '-auto', '-nucmethod', str(method), '-protmethod', str(method),
                    '-gapweight', str(gapw)], check=True, capture_output=True)
    n = len(rows); vals = []
    lines = open(d + '/o.txt').read().split('\n')
    start = [i for i, l in enumerate(lines) if re.match(r'^\s+1\s+2', l)][0] + 1
    for i in range(n - 1):
        f = lines[start + i].split('\t')
        nums = [x.strip() for x in f[1:] if x.strip() and not re.match(r'^s\d+ \d+$', x.strip())]
        # row i holds entries for columns i..n-1 (diagonal first); skip the name field
        row = [x for x in nums]
        vals += row[1:n - i]
    return ','.join(vals)


def clustal(d, rows, kind):
    write_fa(d + '/a.fa', rows)
    subprocess.run(['clustalw', '-convert', '-infile=' + d + '/a.fa', '-outfile=' + d + '/a.aln', '-output=clustal'] + (['-type=' + kind] if kind != 'auto' else []),
                   check=True, capture_output=True)
    lines = open(d + '/a.aln').read().split('\n')
    # first block: rows then conservation line; blocks of 60 columns
    res = []
    i = 3
    n = len(rows)
    while i < len(lines):
        if lines[i].strip() == '' and not (i + n < len(lines)):
            break
        blk = lines[i:i + n + 1]
        if len(blk) < n + 1 or not blk[0].strip(): i += 1; continue
        width = len(blk[0].split()[-1])
        cons = blk[n][16:16 + width] if len(blk[n]) >= 16 else ''
        res.append(cons.ljust(width))
        i += n + 2
    return ''.join(res)


def bioid(rows):
    from Bio.Align import MultipleSeqAlignment
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Phylo.TreeConstruction import DistanceCalculator
    msa = MultipleSeqAlignment([SeqRecord(Seq(r), id=f's{i + 1}') for i, r in enumerate(rows)])
    dm = DistanceCalculator('identity').get_distance(msa)
    n = len(rows)
    return ','.join(repr(round(100 * (1 - dm[i, j]), 6)) for i in range(n) for j in range(i + 1, n))


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 120
    print('kind\trows\tparam\texpected')
    with tempfile.TemporaryDirectory() as d:
        for t in range(n):
            kind = rnd.choice(['dna', 'protein'])
            rows = make_msa(kind, False)
            j = ','.join(rows)
            for ref in rnd.sample(range(0, len(rows) + 1), min(3, len(rows) + 1)) + [0]:
                print(f'infoalign\t{j}\t{ref}\t{infoalign(d, rows, ref)}')
            for method in (0, 1):
                for gw in (0, 0.5):
                    try:
                        print(f'distmat\t{j}\t{method}:{gw}\t{distmat(d, rows, method, gw)}')
                    except Exception as e:
                        print('# distmat failed', e, file=sys.stderr)
            cl = clustal(d, rows, kind)
            print(f'clustal\t{j}\t{kind}\t{cl.replace(" ", "|")}')
            srows = make_msa(kind, True)
            if kind == 'dna':
                print(f'bioid\t{",".join(srows)}\t\t{bioid(srows)}')
        # fixed protein conservation example
        rows = ['MKVSLAKAG-HST', 'MKVALAWAGQHSA', 'MRVGLADAGNHNV', 'MKVCLAEAG-HTA']
        print(f'clustal\t{",".join(rows)}\tprotein\t{clustal(d, rows, "protein").replace(" ", "|")}')


main()
