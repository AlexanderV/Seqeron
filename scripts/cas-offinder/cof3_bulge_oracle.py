#!/usr/bin/env python3
"""Source-faithful port of the Cas-OFFinder 3.0.0b3 bulge enumeration and matcher (oracle for B08 H2/H3).

Ported line by line from `cas-offinder.cpp` of Cas-OFFinder 3.0.0b3 (https://github.com/snugel/cas-offinder):
`parseInput` (the compare-string construction loops, DNA gap offsets j in [first, last+1], RNA removed-run starts
j in [first, last-size]; a DNA compare equal to the ungapped query is reported as type X, size 0) and the output
formatting block (offset / substr arithmetic for the PAM-after and the reversed-PAM branches). The OpenCL kernels
are replaced by their documented semantics: the finder keeps every window of length dnab + len(pattern) whose
pattern positions match (N = wildcard, IUPAC) on either strand; the comparer counts mismatches at the non-N
positions of a compare string and keeps the window when the count is <= the threshold.

Validation: with the toy genomes of tests/.../CrisprDesigner_OffTargetBulge_Tests.cs it reproduces all 8 rows the
real 3.0.0b3 binary printed (docs/Evidence/CRISPR-OFF-001-Evidence.md section D19): 4 exact rows (DNA/RNA, both
strands), the 1-mismatch DNA row and the 3 two-base rows. The differential in B08.md F92 compares this port with
CrisprDesigner.FindOffTargetsWithBulges on 8000 random cases (four systems, 0-4 mismatches, 0-2 bulges).

Usage: cof3_bulge_oracle.py GUIDE GENOME MAX_MISMATCHES MAX_DNA_BULGE MAX_RNA_BULGE
Output record: (is_forward, forward PAM start, kind, size, gap offset in guide coordinates, mismatches,
aligned guide, aligned target); kind in {X, DNA, RNA}; gap = -1 for X; the PAM columns are dropped.
The window padding is the tool's own: a PAM whose padded window crosses a sequence end is not matched.
"""
import sys

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'N': 'ACGT', 'R': 'AG', 'Y': 'CT', 'S': 'CG', 'W': 'AT',
         'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG'}


def rc(s):
    return ''.join(COMP[c] for c in reversed(s))


def compares(s0, dnab, rnab, reversed_pam):
    """Port of parseInput's loops (cas-offinder.cpp). s0 = guide + 'N'*pamlen (PAM after) or 'N'*pamlen + guide
    (PAM before = the tool's is_reversed_pam). Returns [(compare, (bulge_str, index))] in the tool's order."""
    out = []
    first = len(s0) - len(s0.lstrip('N'))
    last = len(s0.rstrip('N')) - 1  # find_last_not_of('N')
    tmp = (s0 + 'N' * dnab) if reversed_pam else ('N' * dnab + s0)
    if dnab == 0:
        out.append((s0, ("0", 0)))
    else:
        for i in range(1, dnab + 1):
            pre = dnab - i
            for j in range(first, last + 2):  # j < find_last_not_of('N') + 2
                if reversed_pam:
                    compare = s0[:j] + 'N' * i + s0[j:] + 'N' * pre
                else:
                    compare = 'N' * pre + s0[:j] + 'N' * i + s0[j:]
                bi = ("0", -j if reversed_pam else j) if compare == tmp else (str(i), -j if reversed_pam else j)
                out.append((compare, bi))
    for i in range(1, rnab + 1):
        pre = dnab + i
        for j in range(first, last + 1 - i):  # j < find_last_not_of('N') + 1 - i
            bi = (s0[j:j + i], -j if reversed_pam else j)
            if reversed_pam:
                compare = s0[:j] + s0[j + i:] + 'N' * pre
            else:
                compare = 'N' * pre + s0[:j] + s0[j + i:]
            out.append((compare, bi))
    return out


def match(window, compare):
    mm = 0
    for w, c in zip(window, compare):
        if c == 'N':
            continue
        if w not in IUPAC[c]:
            mm += 1
    return mm


def run(guide, genome, threshold, dnab, rnab, pam='NGG', pam_after=True):
    pamlen = len(pam)
    reversed_pam = not pam_after
    s0 = (guide + 'N' * pamlen) if pam_after else ('N' * pamlen + guide)
    pattern = ('N' * dnab + 'N' * len(guide) + pam) if pam_after else (pam + 'N' * len(guide) + 'N' * dnab)
    plen = len(pattern)
    comps = compares(s0, dnab, rnab, reversed_pam)
    L = len(genome)
    recs = []
    for strand in (0, 1):
        seq = genome if strand == 0 else rc(genome)
        for w in range(0, len(seq) - plen + 1):
            win = seq[w:w + plen]
            if any(win[k] not in IUPAC[pattern[k]] for k in range(plen)):
                continue
            for compare, (bstr, jj) in comps:
                mm = match(win, compare)
                if mm > threshold:
                    continue
                j = abs(jj)
                if bstr.isdigit():
                    size = int(bstr)
                    if reversed_pam:
                        seq_rna = compare[:j] + '-' * size + compare[j + size:]
                        seq_dna = win
                        if size == 0:
                            seq_rna = seq_rna[:len(seq_rna) - dnab]
                            seq_dna = seq_dna[:len(seq_dna) - dnab]
                    else:
                        offset = dnab - size
                        seq_rna = compare[offset:offset + j] + '-' * size + compare[offset + j + size:]
                        seq_dna = win[offset:]
                    kind = 'X' if size == 0 else 'DNA'
                else:
                    size = len(bstr)
                    if reversed_pam:
                        seq_rna = compare[:j] + bstr + compare[j:]
                        seq_dna = win[:j] + '-' * size + win[j:]
                        seq_rna = seq_rna[:len(seq_rna) - dnab - size]
                        seq_dna = seq_dna[:len(seq_dna) - dnab - size]
                    else:
                        offset = dnab + size
                        seq_rna = compare[offset:offset + j] + bstr + compare[offset + j:]
                        seq_dna = win[offset:offset + j] + '-' * size + win[offset + j:]
                    kind = 'RNA'
                # keep the guide columns only: drop the PAM columns; for PAM-before also drop the DNA-bulge padding
                if pam_after:
                    ag, at = seq_rna[:-pamlen], seq_dna[:-pamlen]
                else:
                    ag, at = seq_rna[pamlen:], seq_dna[pamlen:]
                    if kind == 'DNA':  # drop the trailing padding columns the tool leaves on reversed-PAM DNA rows
                        ag, at = ag[:len(guide) + size], at[:len(guide) + size]
                pam_start = (w + plen - pamlen) if pam_after else w
                fwd_pam = pam_start if strand == 0 else L - pam_start - pamlen
                recs.append((strand == 0, fwd_pam, kind, size, -1 if kind == 'X' else (j if pam_after else j - pamlen), mm, ag, at))
    return recs


if __name__ == '__main__':
    guide, genome = sys.argv[1], sys.argv[2]
    t, d, r = map(int, sys.argv[3:6])
    for rec in sorted(run(guide, genome, t, d, r)):
        print(rec)
