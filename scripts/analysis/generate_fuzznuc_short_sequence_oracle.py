#!/usr/bin/env python3
"""B09 audit A2-18: what EMBOSS fuzznuc 6.6.0 does when the sequence is SHORTER than the pattern and -pmismatch k > 0 is given (the "overhanging hit" of Known_Motif_Search.md).

    python3 scripts/analysis/generate_fuzznuc_short_sequence_oracle.py [n=300] > tests/.../TestData/Analysis/fuzznuc_short_sequence_oracle.tsv

Columns: sequence  pattern  k  complement(Y|N)  hits   where hits = ';'-joined  start,end,strand,mismatches  of the GFF lines (1-based; '-' = no hit).
Pattern: A/C/G/T with, in 20 % of the sequences, IUPAC codes R Y N. Sequence length 1 .. m - 1 (always shorter than the pattern).
"""
import random, re, subprocess, sys, tempfile

rnd = random.Random(20261119)
n = int(sys.argv[1]) if len(sys.argv) > 1 else 300
print('sequence\tpattern\tk\tcomplement\thits')
with tempfile.TemporaryDirectory() as d:
    for _ in range(n):
        m = rnd.randint(3, 8)
        ln = rnd.randint(1, m - 1)
        k = rnd.randint(1, m - 1)
        seq = ''.join(rnd.choice('ACGT') for _ in range(ln))
        iupac = rnd.random() < 0.2
        pat = ''.join(rnd.choice('ACGTRYN') if (iupac and rnd.random() < 0.4) else rnd.choice('ACGT') for _ in range(m))
        compl = rnd.random() < 0.5
        open(d + '/s.fa', 'w').write('>s\n' + seq + '\n')
        subprocess.run(['fuzznuc', '-sequence', d + '/s.fa', '-pattern', pat, '-pmismatch', str(k), '-complement', 'Y' if compl else 'N', '-rformat', 'gff', '-outfile', d + '/o.gff', '-auto'], capture_output=True)
        hits = []
        for line in open(d + '/o.gff'):
            if line.startswith('#') or not line.strip():
                continue
            f = line.split('\t')
            mm = re.search(r'\*mismatch (\d+)', line)
            hits.append(f'{f[3]},{f[4]},{f[6]},{mm.group(1) if mm else 0}')
        print('\t'.join([seq, pat, str(k), 'Y' if compl else 'N', ';'.join(hits) if hits else '-']))
