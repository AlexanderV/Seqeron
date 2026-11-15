#!/usr/bin/env python3
"""E7 helper: the ASCII text (ACGT, one byte per symbol) of a 2-bit packed DNA file of `packrand` (four symbols per byte, first symbol in the top two bits, codes 0..3 = A C G T).
usage: python3 -I e7_unpack_packed.py <packed file> <n symbols> <out text file>"""
import sys
src, n, dst = sys.argv[1], int(sys.argv[2]), sys.argv[3]
data = open(src, 'rb').read((n + 3) // 4)
out = bytearray(4 * len(data))
for k in range(4):
    table = bytes(ord('ACGT'[(b >> (6 - 2 * k)) & 3]) for b in range(256))
    out[k::4] = data.translate(table)
open(dst, 'wb').write(bytes(out[:n]))
