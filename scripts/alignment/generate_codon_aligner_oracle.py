#!/usr/bin/env python3
"""B09 WP67 (F307, R14-10b): oracle of Bio.Align.CodonAligner 1.88 (a nucleotide sequence aligned to its protein with frame shifts): score and every alignment as coordinates.
JSONL.GZ rows: id, protein, dna, table, match, mismatch, wildcard, frameshift (the four scores minus_two, minus_one, plus_one, plus_two), score, total (len(alignments)), alignments (the first 300 of the iterator), aligned_score, undefined, error.
Usage: python3 scripts/alignment/generate_codon_aligner_oracle.py OUT.jsonl.gz [n=500] [seed]      (needs Biopython 1.88)"""
import gzip, json, random, sys
from Bio import Align
from Bio.Data import CodonTable

def backtranslate(rnd, protein, table):
    t = CodonTable.generic_by_id[table]
    rev = {}
    for codon, aa in t.forward_table.items(): rev.setdefault(aa, []).append(codon)
    stops = list(t.stop_codons)
    out = []
    for aa in protein:
        if aa == '*': out.append(rnd.choice(stops))
        elif aa in rev: out.append(rnd.choice(rev[aa]))
        else: out.append(''.join(rnd.choice('ACGT') for _ in range(3)))
    return ''.join(out)

SPECIAL = [
    ('MSR', 'ATGTCTCGT', 1, [-3] * 4), ('MSR', 'AUGUCUCGU', 1, [-3] * 4), ('MSR', 'atgtctcgt', 1, [-3] * 4), ('', 'ATG', 1, [-3] * 4), ('M', '', 1, [-3] * 4), ('M', 'A', 1, [-3] * 4), ('M', 'AT', 1, [-3] * 4),
    ('MLGSNESRVCPKLSPS', 'ATGCTGGGCTCGAACGAGTCCGTGTATGCCCTAAGCTGAGCCCGTCG', 1, [-3] * 4), ('FKKKKF', 'TTTAAAAAAAAAAATTT', 1, [-1] * 4),
    ('MSR', 'AT?TCTCGT', 1, [-3] * 4), ('MSR', 'ATGTCTCGN', 1, [-3] * 4), ('MSR', 'ATGTCTCGTU', 1, [-3] * 4), ('MSR', 'ATGTCUCGT', 1, [-3] * 4), ('M*', 'ATGTAA', 1, [-3] * 4), ('MW', 'ATGTGA', 2, [-3] * 4),
    ('MX', 'ATGNNN', 1, [-3] * 4), ('KK', 'AAAAAAAAAAAA', 1, [0] * 4), ('FKKKKF', 'TTTAAAAAAAAAATTT', 1, [-2, -1, -1, -2]), ('MSR', 'GGATGTCTCGTGG', 1, [-3] * 4),
]


def run_row(row):
    protein, dna, table = row['protein'], row['dna'], row['table']
    al = Align.CodonAligner(codon_table=CodonTable.generic_by_id[table])
    al.match_score = row['match']; al.mismatch_score = row['mismatch']; al.wildcard = row['wildcard']
    (al.frameshift_minus_two_score, al.frameshift_minus_one_score, al.frameshift_plus_one_score, al.frameshift_plus_two_score) = row['frameshift']
    try:
        row['score'] = al.score(protein, dna)
        als = al.align(protein, dna)
        row['total'] = len(als)
        row['alignments'] = []
        for k, a in enumerate(als):
            if k >= 300: break
            row['alignments'].append([[int(x) for x in r] for r in a.coordinates.tolist()])
        row['aligned_score'] = als.score
        # the C iterator reads path cells its traceback never wrote when the optimal alignment skips leading residues (and then lists a vertex of the start twice, or a one-vertex 'alignment', for ever): such a row has no defined result
        lst = row['alignments']
        row['undefined'] = any(len(a[0]) < 2 or a[0][0] != 0 or a[0][-1] != len(protein) for a in lst[:max(row['total'], 1)])
    except Exception as e:
        row['error'] = '%s: %s' % (type(e).__name__, e)
    return row

def main():
    out = sys.argv[1]; n = int(sys.argv[2]) if len(sys.argv) > 2 else 500; seed = int(sys.argv[3]) if len(sys.argv) > 3 else 20262501
    rnd = random.Random(seed)
    rows = []
    for k, (protein, dna, table, frameshift) in enumerate(SPECIAL):
        rows.append(run_row({'id': k, 'protein': protein, 'dna': dna, 'table': table, 'match': 1, 'mismatch': 0, 'wildcard': 'X', 'frameshift': frameshift}))
    while len(rows) < n:
        table = rnd.choice([1, 1, 1, 2, 11])
        aas = 'ACDEFGHIKLMNPQRSTVWY'
        protein = ''.join(rnd.choice(aas) for _ in range(rnd.randint(1, 9)))
        if rnd.random() < 0.1: protein += rnd.choice(['*', 'X'])
        dna = backtranslate(rnd, protein, table)
        for _ in range(rnd.choice([0, 0, 1, 1, 2, 3])):
            p = rnd.randrange(len(dna) + 1)
            op = rnd.random()
            if op < 0.4: dna = dna[:p] + ''.join(rnd.choice('ACGT') for _ in range(rnd.choice([1, 2]))) + dna[p:]
            elif op < 0.7: dna = dna[:p] + dna[p + rnd.choice([1, 2]):]
            else: dna = dna[:p] + rnd.choice('ACGT') + dna[p + 1:]
        if rnd.random() < 0.2: dna = ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(0, 4))) + dna + ''.join(rnd.choice('ACGT') for _ in range(rnd.randint(0, 4)))
        if rnd.random() < 0.1: dna = dna.replace('T', 'U')
        if rnd.random() < 0.1 and dna: p = rnd.randrange(len(dna)); dna = dna[:p] + 'N' + dna[p + 1:]
        row = {'id': len(rows), 'protein': protein, 'dna': dna, 'table': table, 'match': rnd.choice([1, 2, 3]), 'mismatch': rnd.choice([0, 0, -1, -2, 0.5]), 'wildcard': rnd.choice(['X', 'X', 'X', 'A']),
               'frameshift': [rnd.choice([-3, -2, -1.5, -1, 0, -5]) for _ in range(4)] if rnd.random() < 0.6 else [-3, -3, -3, -3]}
        if rnd.random() < 0.15: row['frameshift'] = [row['frameshift'][0]] * 4
        rows.append(run_row(row))
    with gzip.open(out, 'wt', encoding='utf-8', newline='\n', compresslevel=9) as fh:
        for r in rows: fh.write(json.dumps(r) + '\n')
    print(len(rows), 'rows; errors', sum(1 for r in rows if 'error' in r), 'multi', sum(1 for r in rows if r.get('total', 0) > 1), 'truncated', sum(1 for r in rows if r.get('total', 0) > 300))

if __name__ == '__main__':
    main()
