#!/usr/bin/env python3
"""ADR-005 section 5: bytes per input character of the container sections for each stored-width option, at n = 10^8, 10^9, 3x10^9 and 6.2x10^9.
Arithmetic from the section layouts (ADR-003 section 4, ADR-004 section 5a) with the escape shares measured in results/child_stats.txt; the measured anchors are printed below the table:
a container of 10^8 random DNA symbols is 1,325,019,393 bytes = 13.250 B/char (text 0.25, SA 4, LCP 1, CHILD 4, PSI 4), the genome-like text 13.486 and the repeat text 14.250
(analysis/s1b-handoff/reference_hashes.tsv). Core = text + SA + LCP + CHILD; Full = Core + PSI (the navigation pack on disk; the LCP pyramid, 0.13 B/char, is built in memory today and not counted).
usage: bchar.py [--tsv out.tsv]   (markdown on stdout)"""
import sys

NS = [10**8, 10**9, 3 * 10**9, 6.2 * 10**9]
TEXT = 0.25            # 2 bits per symbol (b = 2); 0.375 with b = 3
LCP = 1.0              # U8 + escapes (random DNA: ~0; the genome-like container is 0.24 B/char larger and the repeat container 1.0, measured)
D16_ESC = 0.001        # 2-byte delta cells: 7.2-7.4 thousand escapes per 10^8 cells (0.007 %), 8 bytes each in the hash table (grows with n; rounded up)
D8_ESC = 0.05          # 1-byte delta cells: 0.46-0.81 % escapes x 8 bytes (in a hash table of load 0.5: x2), measured

def frac_long(n):
    """Share of the position values that do not fit a 32-bit slot with one escape value (a value at or above 2^32 - 1)."""
    return max(0.0, 1.0 - (2**32 - 1) / n)

def w(n):   return 4.0 if n + 1 <= 2**32 - 1 else 5.0        # the fixed width the format would choose from n (never 8)
def hyb(n): return 4.0 + 16.0 / 64 + 5.0 * frac_long(n)       # slot + 16-byte record per 64 slots + 40-bit long values
OPTIONS = [
    # name, SA, CHILD, PSI (bytes per entry as a function of n), can hold n
    ('w4 fixed 32-bit (today)',          lambda n: 4.0,           lambda n: 4.0,           lambda n: 4.0,           lambda n: n + 1 <= 2**32 - 1),
    ('w5 fixed 40-bit',                  lambda n: 5.0,           lambda n: 5.0,           lambda n: 5.0,           lambda n: n < 2**40),
    ('w8 fixed 64-bit',                  lambda n: 8.0,           lambda n: 8.0,           lambda n: 8.0,           lambda n: True),
    ('hybrid, all three arrays',         hyb,                     hyb,                     hyb,                     lambda n: True),
    ('hybrid SA + PSI, fixed CHILD',     hyb,                     w,                       hyb,                     lambda n: True),
    ('p34 (34 bits, 15 per line)',       lambda n: 64.0 / 15,     lambda n: 64.0 / 15,     lambda n: 64.0 / 15,     lambda n: n < 2**34),
    ('w by n (4, then 5) + CHILD d16',   w,                       lambda n: 2.0 + D16_ESC, w,                       lambda n: n < 2**40),
    ('w by n (4, then 5) + CHILD d8',    w,                       lambda n: 1.0 + D8_ESC,  w,                       lambda n: n < 2**40),
    ('w by n (4, then 5), CHILD as is',  w,                       w,                       w,                       lambda n: n < 2**40),
]
out = []
print('| option | ' + ' | '.join(f'n = {n:.3g}: Core / Full' for n in NS) + ' |')
print('|---|' + '---|' * len(NS))
for name, sa, ch, psi, ok in OPTIONS:
    cells = []
    for n in NS:
        core = TEXT + sa(n) + LCP + ch(n)
        full = core + psi(n)
        fits = ok(n)
        cells.append(f'{core:.2f} / {full:.2f}' + ('' if fits else ' (cannot hold n)'))
        out.append((name, n, round(core, 3), round(full, 3), round(sa(n), 3), round(ch(n), 3), round(psi(n), 3), fits))
    print(f'| {name} | ' + ' | '.join(cells) + ' |')
print()
print('Per section at the recommended option (w by n + CHILD d16): text 0.25, SA 4 (n <= 4.29x10^9) or 5, LCP 1.0, CHILD 2.0, PSI 4 or 5; + 0.125 text with b = 3.')
if '--tsv' in sys.argv:
    with open(sys.argv[sys.argv.index('--tsv') + 1], 'w') as f:
        f.write('option\tn\tcore_B_per_char\tfull_B_per_char\tsa_B\tchild_B\tpsi_B\tcan_hold_n\n')
        for r in out:
            f.write('\t'.join(str(x) for x in r) + '\n')
