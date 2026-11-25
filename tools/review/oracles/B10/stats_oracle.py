"""Oracle for ASSEMBLY-STATS-001 CalculateNxCurve / CalculateAuN and ASSEMBLY-GAPS-001 ExtractContigs numbering.

Reference code (fetched 2026-11-25, ablab/quast master): quast_libs/N50.py (NG50_and_LG50, au_metric,
transcribed below) and quast_libs/qutils.py split_by_ns (numbering: every piece increments total_contigs,
--min-contig filters afterwards in correct_fasta).
Usage: python3 -I stats_oracle.py   (prints the oracle values locked in the unit tests)
"""
from fractions import Fraction


def ng50_and_lg50(numlist, reference_length, percentage=50.0):
    s = reference_length
    limit = reference_length * (100.0 - percentage) / 100.0
    lg = 0
    for l in numlist:
        s -= l
        lg += 1
        if s <= limit:
            return l, lg
    return None, None


def au_metric(numlist):  # QUAST: float(sum n**2) / float(sum)
    return float(sum(n ** 2 for n in numlist)) / float(sum(numlist))


def au_exact(numlist):  # correctly rounded exact rational
    return float(Fraction(sum(n * n for n in numlist), sum(numlist)))


def curve(lengths, thresholds=(10, 20, 30, 40, 50, 60, 70, 80, 90)):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    out = []
    for t in sorted(thresholds):
        nx, lx = ng50_and_lg50(lengths, total, float(t))
        out.append((t, nx or 0, lx or 0))
    return out


if __name__ == "__main__":
    asm_b = [80, 70, 50, 40, 30, 20, 10, 5]
    print("curve AssemblyB", curve(asm_b))
    print("curve 500..100", curve([500, 400, 300, 200, 100]))
    big = [2_000_000_000, 1_999_999_999, 1_500_000_001, 3, 1]       # total > 2^32, sum of squares > 2^53 (and > 2^64)
    print("big total", sum(big), "sumsq", sum(n * n for n in big))
    print("auN big quast", repr(au_metric(big)), "exact", repr(au_exact(big)))
    big2 = [2_147_483_647] * 50 + [1]
    print("auN big2 quast", repr(au_metric(big2)), "exact", repr(au_exact(big2)))
    # float-sum precision loss example: a 2^53+1 style sum of squares
    small = [94_906_267, 94_906_267, 1]                          # 94906267^2 = 9007199326062755 > 2^53
    print("auN small quast", repr(au_metric(small)), "exact", repr(au_exact(small)), sum(n * n for n in small))
