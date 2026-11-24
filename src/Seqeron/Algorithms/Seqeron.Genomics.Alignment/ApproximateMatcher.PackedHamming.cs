using System.Numerics;

namespace Seqeron.Genomics.Alignment
{
    public static partial class ApproximateMatcher
    {
        // Bit-parallel Hamming scan (B09 WP78, R14-18): the window of m <= 32 nucleotides is kept in one 64-bit word, two bits per base, and compared
        // with the pattern by XOR + popcount, the "2-bit XOR / popcount" engine of read mappers (e.g. Navarro 2001 §6, bit-parallel counters; Hyyrö 2003).
        // Exactly the result of the character loop of FindWithMismatchesCore (same windows, same counts, same positions, same order); only the cost of one window
        // changes from O(m) character tests to O(1) word operations.
        //   * literal pattern of A/C/G/T (and an IUPAC pattern of A/C/G/T only, which is the same thing): window code (2 bits per base) XOR pattern code; a window
        //     position holding a symbol other than A/C/G/T is a mismatch in both modes (the pattern symbol is A/C/G/T, the sequence side is literal), tracked
        //     by a second rolling word of "bad" flags.
        //   * IUPAC pattern with m <= 16: one nibble per position, the window holds the one-hot code of its base (A = 1, C = 2, G = 4, T = 8, any other symbol = 0), the
        //     pattern holds the base set of its code; a position matches when the two nibbles share a bit (IupacHelper.MatchesIupac).
        // Everything else (a literal pattern with a symbol other than A/C/G/T, patterns longer than the word) keeps the character loop.
        private const int PackedMaxLiteral = 32;
        private const int PackedMaxIupac = 16;
        private const ulong LowBitOfPair = 0x5555555555555555UL;
        private const ulong LowBitOfNibble = 0x1111111111111111UL;

        // code of a base by table lookup (A = 0, C = 1, G = 2, T = 3, anything else -1); the scan reads one character per window, so the lookup has no branch
        private static readonly sbyte[] BaseCodes = BuildBaseCodes();

        private static sbyte[] BuildBaseCodes()
        {
            var table = new sbyte[128];
            Array.Fill(table, (sbyte)-1);
            table['A'] = 0;
            table['C'] = 1;
            table['G'] = 2;
            table['T'] = 3;
            return table;
        }

        private static int BaseCode(char c) => c < 128 ? BaseCodes[c] : -1;

        private static int OneHotBase(char c)
        {
            int code = BaseCode(c);
            return code < 0 ? 0 : 1 << code;
        }

        private static int IupacSet(char code) => (IupacHelper.MatchesIupac('A', code) ? 1 : 0) | (IupacHelper.MatchesIupac('C', code) ? 2 : 0) | (IupacHelper.MatchesIupac('G', code) ? 4 : 0) | (IupacHelper.MatchesIupac('T', code) ? 8 : 0);

        /// <summary>True when the pattern can be scanned by the packed engines (see above).</summary>
        private static bool CanScanPacked(string pat, bool iupacPattern, out bool onlyBases)
        {
            onlyBases = true;
            foreach (char c in pat)
            {
                if (BaseCode(c) < 0)
                    onlyBases = false;
            }

            if (onlyBases)
                return pat.Length <= PackedMaxLiteral;
            return iupacPattern && pat.Length <= PackedMaxIupac;
        }

        // The scan state: the hot loop of NextHit keeps the window words in locals and returns at every hit (hits are rare), so the iterator that wraps it (lazy, as the character loop was) does not pay
        // for iterator fields in the per-position work.
        private struct PackedScan
        {
            public int M, MaxMismatches, T, NextCheck;
            public ulong Pattern, WindowMask, LowMask, Window, Bad;
            public bool Nibbles;

            /// <summary>The start of the next window within the mismatch bound, or -1 at the end of the sequence.</summary>
            public int NextHit(string seq, CancellationToken cancellationToken)
            {
                const int checkInterval = 1000;
                ulong w = Window, bad = Bad;
                int t = T, nextCheck = NextCheck, m = M, maxMm = MaxMismatches;
                ulong p = Pattern, windowMask = WindowMask, lowMask = LowMask;
                var table = BaseCodes;
                int found = -1;
                while (t < seq.Length)
                {
                    char c = seq[t];
                    int i = t - m + 1;
                    t++;
                    int count;
                    if (Nibbles)
                    {
                        w = ((w << 4) | (uint)OneHotBase(c)) & windowMask;
                        if (i < 0)
                            continue;
                        ulong both = w & p;
                        count = m - BitOperations.PopCount((both | (both >> 1) | (both >> 2) | (both >> 3)) & lowMask);
                    }
                    else
                    {
                        int code = c < 128 ? table[c] : -1;
                        w = ((w << 2) | (uint)(code & 3)) & windowMask;                 // -1 & 3 = 3: the window code of a bad symbol is irrelevant, the bad word marks it
                        bad = ((bad << 2) | (uint)((code >> 31) & 1)) & windowMask;     // 1 for code -1
                        if (i < 0)
                            continue;
                        ulong d = w ^ p;
                        count = BitOperations.PopCount(((d | (d >> 1)) | bad) & lowMask);
                    }

                    if (i == nextCheck)
                    {
                        cancellationToken.ThrowIfCancellationRequested();
                        nextCheck += checkInterval;
                    }

                    if (count <= maxMm)
                    {
                        found = i;
                        break;
                    }
                }

                Window = w;
                Bad = bad;
                T = t;
                NextCheck = nextCheck;
                return found;
            }
        }

        private static IEnumerable<ApproximateMatchResult> FindWithMismatchesPacked(string seq, string pat, int maxMismatches, bool iupacPattern, bool onlyBases, CancellationToken cancellationToken)
        {
            int m = pat.Length;
            var scan = new PackedScan { M = m, MaxMismatches = maxMismatches, Nibbles = !onlyBases };
            if (onlyBases)
            {
                scan.WindowMask = m == 32 ? ulong.MaxValue : (1UL << (2 * m)) - 1;
                scan.LowMask = LowBitOfPair & scan.WindowMask;
                foreach (char c in pat)
                    scan.Pattern = (scan.Pattern << 2) | (uint)BaseCode(c);
            }
            else
            {
                scan.WindowMask = m == 16 ? ulong.MaxValue : (1UL << (4 * m)) - 1;
                scan.LowMask = LowBitOfNibble & scan.WindowMask;
                foreach (char c in pat)
                    scan.Pattern = (scan.Pattern << 4) | (uint)IupacSet(c);
            }

            int hit;
            while ((hit = scan.NextHit(seq, cancellationToken)) >= 0)
                yield return PackedHit(seq, pat, hit, maxMismatches, iupacPattern);
        }

        private static ApproximateMatchResult PackedHit(string seq, string pat, int i, int maxMismatches, bool iupacPattern)
        {
            var positions = CollectMismatchPositions(seq, i, pat, maxMismatches, iupacPattern);
            return new ApproximateMatchResult(i, seq.Substring(i, pat.Length), positions.Count, positions.AsReadOnly(), MismatchType.Substitution);
        }
    }
}
