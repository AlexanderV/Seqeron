#!/usr/bin/env python3
# Independent check of a BWT file of E9 (shares no code with the program): BWT[j] = dense(T[SA[j] - 1]) for every row j, the dummy byte 0 in the row of the suffix 0, and the meta file (n, primary row, symbol counts).
# usage: python3 -I check_bwt_e9.py <text file (one byte per symbol)> <sa file> <bwt file>      (the meta is <bwt file>.meta)
# Output: one line "bwt_check <yes|NO> n=... primary=... rows=... mismatches=..." ; exit code 0 only when everything is equal.
import array, sys

text_path, sa_path, bwt_path = sys.argv[1:4]
text = open(text_path, 'rb').read()
n = len(text)
letters = sorted(set(text))
dense = [0] * 256
for i, b in enumerate(letters):
    dense[b] = i
meta = dict(l.split('=', 1) for l in open(bwt_path + '.meta').read().split('\n') if '=' in l)
counts = [0] * len(letters)
primary = -1
rows = 0
mismatches = 0
CH = 1 << 20
with open(sa_path, 'rb') as fs, open(bwt_path, 'rb') as fb:
    while True:
        a = array.array('I')
        try:
            a.fromfile(fs, CH)
        except EOFError:
            pass
        if not a:
            break
        got = fb.read(len(a))
        exp = bytearray(dense[text[p - 1]] for p in a)   # p = 0 reads text[-1]: patched below
        for k, p in enumerate(a):
            if p == 0:
                if primary >= 0:
                    mismatches += 1
                primary = rows + k
                exp[k] = 0
            else:
                counts[exp[k]] += 1
        if bytes(exp) != got:
            mismatches += sum(1 for x, y in zip(exp, got) if x != y)
        rows += len(a)
ok = rows == n + 1 and mismatches == 0 and int(meta['n']) == n and int(meta['primary']) == primary and meta['counts'] == ','.join(map(str, counts)) and int(meta['sigma']) == len(letters)
import os
ok = ok and os.path.getsize(bwt_path) == n + 1
print('bwt_check %s n=%d primary=%d rows=%d mismatches=%d sigma=%d' % ('yes' if ok else 'NO', n, primary, rows, mismatches, len(letters)))
sys.exit(0 if ok else 1)
