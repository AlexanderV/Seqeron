#!/usr/bin/env python3
"""Exact check of an E4 output for the text a^n (no oracle needed): the suffix array is n, n-1, ..., 0 and the LCP is 0, 0, 1, 2, ..., n-1 (LCP[r] = r - 1 for r >= 1, LCP[0] = 0).
usage: check_allsame.py <n> <sa file> <lcp base file> [<overflow file>]   exit code 0 when every entry is as expected (4-byte base entries = kind U32, no overflow file content);
for a width of 1 or 2 bytes the escape slots are resolved through the overflow pairs (ADR-003 kind 5)."""
import array, sys

n = int(sys.argv[1])
B = 1 << 20


def read_sa(path):
    with open(path, 'rb') as f:
        r = 0
        while r <= n:
            m = min(B, n + 1 - r)
            a = array.array('I')
            a.fromfile(f, m)
            exp = array.array('I', range(n - r, n - r - m, -1))
            if a != exp:
                for i in range(m):
                    if a[i] != exp[i]:
                        print(f'sa differs at rank {r + i}: {a[i]} instead of {exp[i]}')
                        return False
            r += m
    return True


def read_lcp(path, overflow):
    import os
    size = os.path.getsize(path)
    width = size // (n + 1)
    if size != (n + 1) * width or width not in (1, 2, 4):
        print('unexpected base file size', size)
        return False
    code = {1: 'B', 2: 'H', 4: 'I'}[width]
    escape = {1: 255, 2: 65535, 4: None}[width]
    pairs = array.array('I')
    if overflow and os.path.getsize(overflow) > 0:
        with open(overflow, 'rb') as f:
            pairs.fromfile(f, os.path.getsize(overflow) // 4)
    k = 0
    with open(path, 'rb') as f:
        r = 0
        while r <= n:
            m = min(B, n + 1 - r)
            a = array.array(code)
            a.fromfile(f, m)
            exp = array.array('I', [0, 0] + list(range(1, m - 1))) if r == 0 else array.array('I', range(r - 1, r - 1 + m))
            if escape is not None:
                got = array.array('I', a)
                for i in range(m):
                    if got[i] == escape:
                        if 2 * k + 1 >= len(pairs) or pairs[2 * k] != r + i:
                            print(f'missing or misplaced overflow pair for rank {r + i}')
                            return False
                        got[i] = pairs[2 * k + 1]
                        k += 1
                a = got
            if array.array('I', a) != exp:
                for i in range(m):
                    if a[i] != exp[i]:
                        print(f'lcp differs at rank {r + i}: {a[i]} instead of {exp[i]}')
                        return False
            r += m
    return k == len(pairs) // 2


ok = read_sa(sys.argv[2]) and read_lcp(sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
print('yes' if ok else 'NO')
sys.exit(0 if ok else 1)
