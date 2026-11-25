"""Line-by-line port of EMBOSS 6.x emboss/merger.c (merger_Merge, merger_bestquality,
merger_quality) for an ungapped suffix/prefix alignment: m = a[start1:], n = b[:olen]."""
import sys

GOOD = "aAcCgGtTuU"
AMB = "mrwsykvhdbMRWSYKVHDB"

def quality(seq, pos, window):
    value = 0
    i = pos
    while i < pos + window and i < len(seq):
        if seq[i] in GOOD: value += 2
        elif seq[i] in AMB: value += 1
        i += 1
    jlast = pos - window
    j = pos - 1
    while j > jlast and j >= 0:
        if seq[j] in GOOD: value += 2
        elif seq[j] in AMB: value += 1
        j -= 1
    import numpy as np
    return np.float32(value / (window * 2.0 + 1.0))

def bestquality(a, b, apos, bpos):
    qa, qb = quality(a, apos, 5), quality(b, bpos, 5)
    if qa == qb: qa, qb = quality(a, apos, 10), quality(b, bpos, 10)
    if qa == qb: qa, qb = quality(a, apos, 20), quality(b, bpos, 20)
    return qa >= qb

def merge(c1, c2, l):
    a, b = c1.lower(), c2.lower()           # ajSeqFmtLower
    start1, start2 = len(a) - l, 0
    p, q = list(a[start1:]), list(b[:l])   # ungapped walk alignment
    ms = []
    if start1 > start2: ms.extend(a[:start1])
    apos, bpos = start1, start2
    for i in range(l):
        if p[i] in 'nN': ms.append(q[i].upper())
        elif q[i] in 'nN': ms.append(p[i].upper())
        elif p[i] != q[i]:
            ms.append(p[i].upper() if bestquality(a, b, apos, bpos) else q[i].upper())
        else: ms.append(p[i])
        apos += 1; bpos += 1
    alen, blen = len(a) - apos, len(b) - bpos
    if alen > blen: ms.append(a[apos:])
    if blen > alen: ms.append(b[bpos:])
    elif alen and blen: ms.append(a[apos:])
    return "".join(ms)

if __name__ == "__main__":
    for line in sys.stdin:
        if not line.strip(): continue
        c1, c2, l = line.split()
        r = merge(c1, c2, int(l))
        print(c1, c2, l, "->", r, r.upper())
