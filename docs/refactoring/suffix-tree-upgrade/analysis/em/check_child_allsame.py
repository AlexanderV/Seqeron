#!/usr/bin/env python3
"""The closed form of the child table of a^n (used where the in-RAM reference does not fit, a^n 200 M): the tagged codec of layout A holds, for the text a^n, the field Down | (k + 1) at every rank k in 0 .. n - 1
(the interval pushed at rank k + 1 has left boundary k and first l-index k + 1; the root has k = 0) and 0 at rank n. Layout F: PSV[0] = PSV[1] = 0xFFFFFFFF, PSV[k] = k - 1 for k >= 2; NSV[k] = n + 1 for k <= 1, and NSV[k] = n + 1
for the strictly increasing LCP of a^n (there is no smaller value to the right); the check takes the layout from the file extension (.cld tagged, .psv, .nsv).
usage: check_child_allsame.py <file> <n>"""
import sys
from array import array

path, n = sys.argv[1], int(sys.argv[2])
ext = path.rsplit('.', 1)[-1]
DOWN = 3 << 30
bad = 0
chunk = 1 << 20
with open(path, 'rb') as f:
    total = (n + 1) * 4
    f.seek(0, 2)
    if f.tell() != total:
        print('length %d, expected %d' % (f.tell(), total))
        sys.exit(1)
    f.seek(0)
    k = 0
    while k <= n:
        count = min(chunk, n + 1 - k)
        got = array('I')
        got.fromfile(f, count)
        if ext == 'cld':
            exp = array('I', (DOWN | (j + 1) for j in range(k, min(k + count, n))))
            if k + count > n:
                exp.append(0)
        elif ext == 'psv':
            exp = array('I', (0xFFFFFFFF if j <= 1 else j - 1 for j in range(k, k + count)))
        elif ext == 'nsv':
            exp = array('I', [n + 1]) * count
        else:
            print('unknown extension')
            sys.exit(2)
        if got != exp:
            for j in range(count):
                if got[j] != exp[j]:
                    print('first difference at rank %d: %d, expected %d' % (k + j, got[j], exp[j]))
                    break
            bad = 1
            break
        k += count
print('equal' if not bad else 'DIFFERENT')
sys.exit(bad)
