#!/usr/bin/env python3
"""B09 WP67 (F308, R14-13): an INDEPENDENT computation of the information-content values of Biopython 1.79's test suite (Tests/test_AlignInfo.py), written from the definition (not from SummaryInfo): per column the frequency f of each
letter (gaps and the ignored characters not counted; with a pseudo count k: (count + e k) / (total + k)) and sum f log(f / e) / log(base). It reproduces the published 7.32029999423075, 133.061475107 and 7.546 and the per-column vector of the
PSSM theory example, and the hand-derived values of AlignInfo_Tests.cs. Usage: python3 scripts/alignment/alignment_summary_info_reference.py"""
import math


def ic(rows, e, base=2.0, ignore=('-',), pseudo=0.0):
    total_ic = 0.0
    vector = []
    for c in range(len(rows[0])):
        col = [r[c] for r in rows if r[c] not in ignore]
        letters = sorted(set(''.join(rows)) - set(ignore))
        n = float(len(col))
        col_ic = 0.0
        if n == 0:                      # a column of ignored characters adds nothing
            vector.append(col_ic)
            continue
        for a in letters:
            f = col.count(a) / n if not pseudo else (col.count(a) + e[a] * pseudo) / (n + pseudo)
            if f > 0:
                col_ic += f * math.log(f / e[a]) / math.log(base)
        vector.append(col_ic)
        total_ic += col_ic
    return total_ic, vector


if __name__ == '__main__':
    print(ic(['ACGTCGCG', 'GGGGCCCC', 'AAACACAC'], dict.fromkeys('ACGT', 0.25))[0], '(published 7.32029999423075)')
    prot = ['MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-', 'MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*', 'MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*']
    print(ic(prot, dict.fromkeys('ACDEFGHIKLMNPQRSTVWY', 0.05), ignore=('-', '*'))[0], '(published 133.061475107)')
    dna = ['AACCACGTTTAA', 'CACCACGTGGGT', 'CACCACGTTCGC', 'GCGCACGTGGGG', 'TCGCACGTTGTG', 'TGGCACGTGTTT', 'TGACACGTGGGA', 'TTACACGTGCGC']
    total, vec = ic(dna, {'A': 0.325, 'G': 0.175, 'T': 0.325, 'C': 0.175}, base=math.e, pseudo=1)
    print(total, '(published 7.546)', [round(v, 2) for v in vec])
    print(ic(['A', 'A', 'C'], {'A': 0.5, 'C': 0.5})[0], '(hand-derived 0.08170416594551044)')
