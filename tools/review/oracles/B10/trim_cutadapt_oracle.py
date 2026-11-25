"""R1-WP3 A8-1 oracle: cutadapt 5.2 quality_trim_index, 3'-only default (-q C = (0,C)) vs both ends (-q C,C).
Cross-checks the independent running-sum formulation used by AssemblyTrimFuzzTests.ReferenceTrim on 20000 random reads.
Run: python3 -I trim_cutadapt_oracle.py   (needs pip cutadapt==5.2)"""
import random
from cutadapt.qualtrim import quality_trim_index as q

def ref(qual, c, f):
    n = len(qual); end = n; s = 0; m = 0
    if c >= 1:
        for i in range(n - 1, -1, -1):
            s += ord(qual[i]) - 33 - c
            if s > 0: break
            if s < m: m = s; end = i
    start = 0; s = 0; m = 0
    if f >= 1:
        for i in range(n):
            s += ord(qual[i]) - 33 - f
            if s > 0: break
            if s < m: m = s; start = i + 1
    return (0, 0) if start >= end else (start, end)

if __name__ == "__main__":
    for qs, c in [('$#%,()<;IK', 10), ('IIII##', 10), ('KI;<)(,%#$', 10), ('##I##', 10), ('5&&', 20), ('###IIII##', 20)]:
        print(qs, c, "3'-only (0,C):", q(qs, 0, c, 33), " both (C,C):", q(qs, c, c, 33))
    r = random.Random(1); bad = 0
    for _ in range(20000):
        n = r.randint(0, 60); lo = r.randint(0, 44); hi = min(93, lo + r.randint(0, 44))
        qs = ''.join(chr(33 + r.randint(lo, hi)) for _ in range(n)); c = r.randint(-3, 44)
        for f in (0, c):
            if ref(qs, c, f) != q(qs, f, c, 33): bad += 1
    print("mismatches vs cutadapt:", bad)
