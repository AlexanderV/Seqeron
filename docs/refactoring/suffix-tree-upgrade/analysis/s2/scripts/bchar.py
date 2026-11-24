#!/usr/bin/env python3
"""ADR-005 section 5: bytes per input character of the container sections for each stored-width option, at n = 10^8, 10^9, 3x10^9 and 6.2x10^9.
Pure arithmetic from the section layouts (ADR-003 section 4, ADR-004 section 5a); the measured anchors are printed beside it: a container of 10^8 random DNA symbols is 1,325,019,393 bytes = 13.250 B/char
(text 0.25, SA 4, LCP 1, CHILD 4, PSI 4), the genome-like text 13.486 and the repeat text 14.250 (analysis/s1b-handoff/reference_hashes.tsv).
usage: bchar.py [--tsv out.tsv]"""
import sys, math

NS = [10**8, 10**9, 3 * 10**9, 6.2 * 10**9]
TEXT = 0.25            # 2 bits per symbol (b = 2); 0.375 with b = 3
LCP = 1.0              # U8 + escapes (random DNA: ~0; genome-like +0.24 and repeat text +1.0 in the measured containers)
RMQ = 0.13             # the LCP pyramid of the navigation pack (fanout 32), built in memory today
BUCKET = 0.0

def frac_long(n):
    """Share of the position values that do not fit a 32-bit slot with one escape value (a value at or above 2^32 - 1)."""
    return max(0.0, 1.0 - (2**32 - 1) / n)

def per_array(opt, n):
    """bytes per entry of one position array (SA, CHILD or PSI) and whether it can hold n"""
    if opt == 'w4':   return 4.0, n + 1 <= 2**32 - 1
    if opt == 'w5':   return 5.0, n < 2**40
    if opt == 'w8':   return 8.0, True
    if opt == 'hyb':  # 32-bit slots + a 16-byte record (bitmap, cumulative count) per 64 slots + 5-byte long values
        return 4.0 + 16.0 / 64 + 5.0 * frac_long(n), True
    if opt == 'p34':  return 64.0 / 15, n < 2**34
    raise ValueError(opt)

def row(opt, n):
    e, ok = per_array(opt, n)
    core = TEXT + e + LCP + e          # text + SA + LCP + CHILD
    full = core + e + RMQ              # + PSI + pyramid (the navigation pack on disk)
    return core, full, ok, e

out = []
print(f"{'option':8}" + ''.join(f"  n={n:.3g}".rjust(26) for n in NS))
print(' ' * 8 + ''.join('       core /  full / SA'.rjust(26) for _ in NS))
for opt in ['w4', 'w5', 'w8', 'hyb', 'p34']:
    line = f"{opt:8}"
    for n in NS:
        core, full, ok, e = row(opt, n)
        cell = f"{core:6.2f} / {full:5.2f} / {e:4.2f}" + ('' if ok else ' (n/a)')
        line += cell.rjust(26)
        out.append((opt, n, round(core, 3), round(full, 3), round(e, 3), ok))
    print(line)
if '--tsv' in sys.argv:
    with open(sys.argv[sys.argv.index('--tsv') + 1], 'w') as f:
        f.write('option\tn\tcore_B_per_char\tfull_B_per_char\tbytes_per_position_entry\tcan_hold_n\n')
        for r in out:
            f.write('\t'.join(str(x) for x in r) + '\n')
# the extra rows the option "delta-coded CHILD" needs measured data for (stats command): not computed here
