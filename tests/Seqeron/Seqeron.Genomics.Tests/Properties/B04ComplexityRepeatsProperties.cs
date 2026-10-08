using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the behaviours introduced by review batch B04
/// (docs/Validation/review-2026-09/B04.md, fixes F1–F20) in <see cref="SequenceComplexity"/> and
/// <see cref="RepeatFinder"/>: LZ76 exhaustive-history complexity (F1), SDUST masking (F3),
/// low-complexity regions as the union of flagged windows (F4), maximal primitive ACGT-only
/// microsatellites (F5/F6), ACGT-only pairing in inverted / direct repeats and palindromes
/// (F10–F12), the palindrome length cap (F13), the TRF copy-number rule (F14/F15), the linguistic-
/// complexity alphabet (F20) and the entropy kernel (F19).
/// Each property compares against a brute-force oracle written from the sourced definition, or asserts
/// a mathematically exact invariant. Sizes are small so the fixture stays fast.
///
/// Test Units: SEQ-COMPLEX-COMPRESS-001, SEQ-COMPLEX-DUST-001, SEQ-COMPLEX-KMER-001, SEQ-COMPLEX-WINDOW-001,
/// REP-STR-001, REP-INV-001, REP-DIRECT-001, REP-PALIN-001, REP-APPROX-001, SEQ-COMPLEX-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Complexity")]
public class B04ComplexityRepeatsProperties
{
    private static Arbitrary<string> StringOver(string alphabet, int minLen, int maxLen) =>
        (from n in Gen.Choose(minLen, maxLen)
         from chars in Gen.Elements(alphabet.ToCharArray()).ArrayOf(n)
         select new string(chars)).ToArbitrary();

    // Low-entropy (repeat-rich) DNA: each string draws from a random 1–4 letter sub-alphabet.
    private static Gen<string> RepeatRichDna(int minLen, int maxLen) =>
        from alpha in Gen.Elements("A", "AC", "AT", "CG", "ACG", "ACGT")
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
        select new string(chars);

    private static bool IsAcgt(char c) => c is 'A' or 'C' or 'G' or 'T';

    private static char Complement(char c) => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => '?' };

    private static bool Pairs(char x, char y) => IsAcgt(x) && IsAcgt(y) && Complement(x) == y;

    private static string ReverseComplementAcgt(string s) => new(s.Reverse().Select(Complement).ToArray());

    #region F1 — LZ76 exhaustive history

    // Lempel & Ziv (1976) exhaustive history, straight from the definition: the component starting at p is
    // extended while S[p..p+len) occurs in S[0..p+len−1) (overlap allowed). Cubic; test sizes only.
    private static int BruteForceLz76(string s)
    {
        int n = s.Length, p = 0, c = 0;
        while (p < n)
        {
            int len = 1;
            while (p + len <= n && s[..(p + len - 1)].Contains(s.Substring(p, len), StringComparison.Ordinal))
                len++;
            c++;
            p += len;
        }
        return c;
    }

    /// <summary>F1: c(S) equals the brute-force exhaustive-history count on binary / 2- / 4-letter strings, and 0 ≤ c ≤ n.</summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property LempelZiv_EqualsExhaustiveHistoryDefinition_AndAtMostN()
    {
        var gen = from alpha in Gen.Elements("01", "AC", "ACGT")
                  from n in Gen.Choose(0, 60)
                  from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
                  select new string(chars);
        return Prop.ForAll(gen.ToArbitrary(), s =>
        {
            int c = SequenceComplexity.CalculateLempelZivComplexity(s);
            return (c == BruteForceLz76(s) && c >= 0 && c <= s.Length && (s.Length == 0 || c >= 1))
                .Label($"LZ76('{s}') = {c}, brute force {BruteForceLz76(s)}");
        });
    }

    /// <summary>
    /// F1: prefix extension — appending one symbol either lengthens the last component or opens one new
    /// component, so c(S[..k+1]) ∈ {c(S[..k]), c(S[..k]) + 1} for every k (an LZ78 phrase parse also has this,
    /// but together with the brute-force equality it pins the LZ76 history).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property LempelZiv_PrefixExtension_MonotoneByAtMostOne()
    {
        return Prop.ForAll(StringOver("ACGT", 1, 50), s =>
        {
            bool ok = true;
            int prev = 0;
            for (int k = 1; k <= s.Length; k++)
            {
                int c = SequenceComplexity.CalculateLempelZivComplexity(s[..k]);
                ok &= c == prev || c == prev + 1;
                prev = c;
            }
            return ok.Label($"LZ76 prefix sequence not unit-step monotone for '{s}'");
        });
    }

    #endregion

    #region F3 — SDUST masking soundness (Morgulis et al. 2006; lh3/sdust)

    /// <summary>
    /// Coverage of all "high-scoring" intervals: bases lying in some interval x of 4..W bases (ℓ ≥ 2 triplets)
    /// whose DUST score r(x)/(ℓ−1), r = Σ c(c−1)/2, is strictly greater than T.
    /// </summary>
    private static bool[] HighScoringCoverage(string s, int w, double threshold)
    {
        int n = s.Length;
        var cover = new bool[n];
        for (int st = 0; st < n; st++)
        {
            var counts = new Dictionary<string, int>();
            int r = 0;
            for (int len = 3; len <= Math.Min(w, n - st); len++)
            {
                string t = s.Substring(st + len - 3, 3);
                counts.TryGetValue(t, out int c);
                r += c;
                counts[t] = c + 1;
                int l = len - 3; // ℓ − 1
                if (l >= 1 && r > threshold * l)
                    for (int i = st; i < st + len; i++) cover[i] = true;
            }
        }
        return cover;
    }

    /// <summary>
    /// F3 (soundness of the SDUST port): every base masked by MaskLowComplexity lies inside an interval of at most
    /// W bases whose DUST score (divisor ℓ − 1, F2) strictly exceeds T — SDUST masks only unions of perfect
    /// intervals, each of which is such an interval — and every unmasked base is returned unchanged (same length).
    /// A sequence without any such interval is therefore returned verbatim.
    /// (Completeness against a naive "all perfect intervals" enumeration is NOT asserted: sdust's suffix pruning
    /// <c>cv[t] &gt; 2T</c> is part of the reference algorithm; exact agreement with compiled sdust is locked by the
    /// unit rows in SequenceComplexity_CalculateDustScore_Tests.)
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property Sdust_MaskedBases_LieInIntervalsScoringAboveThreshold()
    {
        var gen = from s in RepeatRichDna(0, 60)
                  from w in Gen.Choose(3, 24)
                  from t in Gen.Elements(0.5, 1.0, 1.5, 2.0, 3.0)
                  select (s, w, t);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            string masked = SequenceComplexity.MaskLowComplexity(new DnaSequence(x.s), x.w, x.t, 'x');
            bool[] allowed = HighScoringCoverage(x.s, x.w, x.t);
            bool ok = masked.Length == x.s.Length;
            for (int i = 0; ok && i < x.s.Length; i++)
                ok = masked[i] == 'x' ? allowed[i] : masked[i] == x.s[i];
            return ok.Label($"SDUST('{x.s}', W={x.w}, T={x.t}) = '{masked}' masks a base outside every interval scoring > T");
        });
    }

    #endregion

    #region F4 — low-complexity regions = union of flagged windows (BBDuk maskLowEntropy)

    /// <summary>
    /// F4: FindLowComplexityRegions returns exactly the maximal runs of positions covered by windows whose
    /// entropy is strictly below the threshold: regions ascending, disjoint and non-adjacent, every flagged
    /// window inside a region, every region base inside some flagged window, End inclusive and consistent
    /// with Length/Sequence, MinEntropy = lowest flagged-window entropy of the region.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property LowComplexityRegions_EqualUnionOfFlaggedWindows()
    {
        var gen = from s in RepeatRichDna(1, 70)
                  from w in Gen.Choose(1, 12)
                  from t in Gen.Elements(0.5, 0.9, 1.0, 1.3, 1.6, 2.1)
                  select (s, w, t);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            var (s, w, t) = x;
            var regions = SequenceComplexity.FindLowComplexityRegions(new DnaSequence(s), w, t).ToList();

            // Reference: per-window entropy from the public scalar metric, then a boolean coverage mask.
            var cover = new bool[s.Length];
            var windowEntropy = new double[Math.Max(0, s.Length - w + 1)];
            for (int i = 0; i + w <= s.Length; i++)
            {
                windowEntropy[i] = SequenceComplexity.CalculateShannonEntropy(s.Substring(i, w));
                if (windowEntropy[i] < t)
                    for (int k = i; k < i + w; k++) cover[k] = true;
            }
            var expected = new List<(int Start, int End, double Min)>();
            for (int i = 0; i < s.Length;)
            {
                if (!cover[i]) { i++; continue; }
                int j = i;
                while (j + 1 < s.Length && cover[j + 1]) j++;
                double min = double.MaxValue;
                for (int k = i; k + w - 1 <= j; k++)
                    if (windowEntropy[k] < t) min = Math.Min(min, windowEntropy[k]);
                expected.Add((i, j, min));
                i = j + 1;
            }

            bool ok = regions.Count == expected.Count
                && regions.Zip(expected).All(p =>
                    p.First.Start == p.Second.Start && p.First.End == p.Second.End
                    && p.First.Length == p.Second.End - p.Second.Start + 1
                    && p.First.Sequence == s.Substring(p.Second.Start, p.First.Length)
                    && p.First.MinEntropy == p.Second.Min)
                && regions.Zip(regions.Skip(1)).All(p => p.Second.Start > p.First.End + 1);
            return ok.Label($"LCR('{s}', w={w}, t={t}): got [{string.Join(";", regions.Select(r => $"{r.Start}-{r.End}"))}] " +
                            $"expected [{string.Join(";", expected.Select(r => $"{r.Start}-{r.End}"))}]");
        });
    }

    #endregion

    #region F5/F6 — microsatellites: maximal primitive ACGT-only runs

    private static bool IsPrimitive(string u) =>
        Enumerable.Range(1, u.Length - 1).All(d => u.Length % d != 0
            || string.Concat(Enumerable.Repeat(u[..d], u.Length / d)) != u);

    /// <summary>
    /// F5/F6: on mixed-case ACGT+N input the string overload returns exactly the brute-force set of maximal
    /// primitive runs with ≥ minRepeats complete copies of an A/C/G/T unit — each locus once (no rotations),
    /// never an N run — ordered by unit length, then position.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property Microsatellites_EqualBruteForceMaximalPrimitiveAcgtRuns()
    {
        var gen = from alpha in Gen.Elements("AC", "ACN", "ACGTN", "ATn", "acgt", "N")
                  from n in Gen.Choose(0, 50)
                  from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
                  from minRep in Gen.Choose(2, 4)
                  select (s: new string(chars), minRep);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            string u = x.s.ToUpperInvariant();
            var expected = new List<(int Pos, string Unit, int Count)>();
            for (int p = 1; p <= 6; p++)
                for (int i = 0; i + p <= u.Length; i++)
                {
                    if (i > 0 && i - 1 + p < u.Length && u[i - 1] == u[i - 1 + p]) continue; // not left-maximal
                    int e = i + p;
                    while (e < u.Length && u[e] == u[e - p]) e++;
                    int copies = (e - i) / p;
                    string unit = u.Substring(i, p);
                    if (copies >= x.minRep && unit.All(IsAcgt) && IsPrimitive(unit))
                        expected.Add((i, unit, copies));
                }

            var actual = RepeatFinder.FindMicrosatellites(x.s, 1, 6, x.minRep)
                .Select(r => (r.Position, r.RepeatUnit, r.RepeatCount)).ToList();
            return actual.SequenceEqual(expected)
                .Label($"'{x.s}' minRep={x.minRep}: got [{string.Join(";", actual)}] expected [{string.Join(";", expected)}]");
        });
    }

    #endregion

    #region F10/F11/F12 — ACGT-only pairing / matching

    /// <summary>
    /// F10: on input containing N / IUPAC / U / gaps, every reported inverted repeat has ACGT-only arms, the right
    /// arm is the exact Watson–Crick reverse complement of the left arm, and coordinates match the reported strings.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property InvertedRepeats_NonAcgtNeverPairs()
    {
        return Prop.ForAll(StringOver("ACGTNRYSU-acgtn", 0, 60), s =>
        {
            string u = s.ToUpperInvariant();
            var reps = RepeatFinder.FindInvertedRepeats(s, minArmLength: 2, maxLoopLength: 12, minLoopLength: 0).ToList();
            bool ok = reps.All(r =>
                r.LeftArm.All(IsAcgt) && r.RightArm == ReverseComplementAcgt(r.LeftArm)
                && u.Substring(r.LeftArmStart, r.ArmLength) == r.LeftArm
                && u.Substring(r.RightArmStart, r.ArmLength) == r.RightArm
                && r.RightArmStart - r.LeftArmStart - r.ArmLength == r.LoopLength);
            return ok.Label($"non-ACGT pairing or inconsistent coordinates in '{s}'");
        });
    }

    /// <summary>
    /// F11: FindDirectRepeats returns exactly the brute-force maximal repeated pairs (i &lt; j, left- and right-
    /// maximal, non-ACGT symbols never match, length in [min, max], spacing j − i − L ≥ minSpacing), sorted by
    /// (FirstPosition, SecondPosition) — MUMmer repeat-match -f semantics.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property DirectRepeats_EqualBruteForceMaximalPairs_AcgtOnly()
    {
        var gen = from alpha in Gen.Elements("AC", "ACN", "ACGT", "ACGTNR", "acgT")
                  from n in Gen.Choose(0, 50)
                  from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
                  from minLen in Gen.Choose(2, 5)
                  from maxLen in Gen.Elements(6, 10, int.MaxValue)
                  from minSpacing in Gen.Elements(int.MinValue, -3, 0, 1, 4)
                  select (s: new string(chars), minLen, maxLen, minSpacing);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            string u = x.s.ToUpperInvariant();
            bool Eq(int a, int b) => IsAcgt(u[a]) && u[a] == u[b];
            var expected = new List<(int, int, int)>();
            for (int i = 0; i < u.Length; i++)
                for (int j = i + 1; j < u.Length; j++)
                {
                    if (i > 0 && Eq(i - 1, j - 1)) continue;
                    int len = 0;
                    while (j + len < u.Length && Eq(i + len, j + len)) len++;
                    if (len >= x.minLen && len <= x.maxLen && (long)j - i - len >= x.minSpacing)
                        expected.Add((i, j, len));
                }
            var actual = RepeatFinder.FindDirectRepeats(x.s, x.minLen, x.maxLen, x.minSpacing)
                .Select(r => (r.FirstPosition, r.SecondPosition, r.Length)).ToList();
            return actual.SequenceEqual(expected)
                .Label($"'{x.s}' ({x.minLen},{x.maxLen},{x.minSpacing}): got {actual.Count}, expected {expected.Count}");
        });
    }

    /// <summary>
    /// F12: FindPalindromes (string overload, mixed case, N/IUPAC/U/gap) returns exactly every even window of
    /// length in [min, max] that is an A/C/G/T reverse-complement palindrome (Rosalind REVP), ordered by
    /// (position, length) — complete and sound against brute force.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property Palindromes_EqualBruteForceAcgtReversePalindromes()
    {
        var gen = from alpha in Gen.Elements("AT", "ACGT", "ACGTN", "ACGTNSWRY-U", "gcAT")
                  from n in Gen.Choose(0, 60)
                  from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
                  from minLen in Gen.Elements(4, 6)
                  from extra in Gen.Choose(0, 10)
                  select (s: new string(chars), minLen, maxLen: minLen + extra);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            string u = x.s.ToUpperInvariant();
            var expected = new List<(int, int)>();
            for (int p = 0; p < u.Length; p++)
                for (int len = x.minLen; len <= x.maxLen && p + len <= u.Length; len += 2)
                {
                    bool pal = true;
                    for (int k = 0; pal && k < len / 2; k++) pal = Pairs(u[p + k], u[p + len - 1 - k]);
                    if (pal) expected.Add((p, len));
                }
            var result = RepeatFinder.FindPalindromes(x.s, x.minLen, x.maxLen).ToList();
            bool ok = result.Select(r => (r.Position, r.Length)).SequenceEqual(expected)
                      && result.All(r => r.Sequence == u.Substring(r.Position, r.Length));
            return ok.Label($"palindromes of '{x.s}' [{x.minLen},{x.maxLen}] differ from brute force");
        });
    }

    /// <summary>F13: maxLength beyond the sequence length (up to int.MaxValue) is equivalent to maxLength = n, on both overloads.</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property Palindromes_HugeMaxLength_EquivalentToSequenceLengthCap()
    {
        return Prop.ForAll(StringOver("AT", 4, 40), s =>
        {
            var capped = RepeatFinder.FindPalindromes(s, 4, s.Length).ToList();
            bool ok = RepeatFinder.FindPalindromes(s, 4, int.MaxValue).SequenceEqual(capped)
                      && RepeatFinder.FindPalindromes(new DnaSequence(s), 4, int.MaxValue).SequenceEqual(capped)
                      && RepeatFinder.FindPalindromes(s, 4, int.MaxValue - 1).SequenceEqual(capped);
            return ok.Label($"int.MaxValue maxLength differs from maxLength = n for '{s}'");
        });
    }

    #endregion

    #region F14/F15 — TRF copy-number rule and statistics ranges

    // Planted (possibly mutated) tandem repeat inside random flanks.
    private static Arbitrary<string> PlantedApproximateRepeat() =>
        (from preLen in Gen.Choose(0, 12)
         from pre in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(preLen)
         from unit in Gen.Elements("CAG", "AC", "GATA", "TTAGGG", "A", "ACGTC")
         from k in Gen.Choose(3, 12)
         from mutPos in Gen.Choose(0, 200)
         from mutBase in Gen.Elements('A', 'C', 'G', 'T', 'N')
         from post in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(8)
         let core = string.Concat(Enumerable.Repeat(unit, k))
         let mutated = core.Length > 0 ? core.Remove(mutPos % core.Length, 1).Insert(mutPos % core.Length, mutBase.ToString()) : core
         select new string(pre) + mutated + new string(post)).ToArbitrary();

    /// <summary>
    /// F14/F15/F18: every reported approximate repeat satisfies the TRF minimum copy number (≥ 1.9 aligned copies
    /// for consensus ≤ 50 bp; never below 1.8), lies inside the sequence, has Start ascending, percentages in
    /// [0, 100] with matches + indels ≤ 100, A/C/G/T composition ≤ 100 %, entropy in [0, 2] and score ≥ minScore.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property ApproximateRepeats_TrfCopyRule_AndStatisticsRanges()
    {
        const double eps = 1e-9;
        return Prop.ForAll(PlantedApproximateRepeat(), s =>
        {
            var results = RepeatFinder.FindApproximateTandemRepeats(s, 1, 12, 20).ToList();
            bool ok = results.All(r =>
                r.CopyNumber >= (r.ConsensusSize <= 50 ? 1.9 : 1.8) - eps
                && r.Start >= 0 && r.SpanLength > 0 && r.Start + r.SpanLength <= s.Length
                && r.Period is >= 1 and <= 12
                && r.PercentMatches is >= 0 and <= 100 + eps && r.PercentIndels is >= 0 and <= 100 + eps
                && r.PercentMatches + r.PercentIndels <= 100 + eps
                && r.PercentA + r.PercentC + r.PercentG + r.PercentT <= 100 + eps
                && r.Entropy is >= 0 and <= 2 + eps
                && r.AlignmentScore >= 20
                && r.Consensus.Length == r.ConsensusSize);
            ok &= results.Zip(results.Skip(1)).All(p => p.First.Start <= p.Second.Start);
            return ok.Label($"TRF contract violated for '{s}': {string.Join("; ", results)}");
        });
    }

    #endregion

    #region F19/F20 — entropy kernel and linguistic-complexity alphabet

    /// <summary>
    /// F20: LC ∈ (0, 1] for ANY symbol alphabet (a = |symbols ∪ {A,C,G,T/U}|), on the hash path (m ≤ 12) and the
    /// all-length suffix-tree path; `ACGTN`-style inputs never exceed 1.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property LinguisticComplexity_AnyAlphabet_InUnitInterval()
    {
        var gen = from s in StringOver("ACGTNRYU-xz#", 1, 40).Generator
                  from m in Gen.Elements(1, 2, 3, 6, 12, 13, 40)
                  select (s, m);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            double lc = SequenceComplexity.CalculateLinguisticComplexity(x.s, x.m);
            return (lc > 0 && lc <= 1.0).Label($"LC('{x.s}', m={x.m}) = {lc:R} outside (0,1]");
        });
    }

    /// <summary>
    /// F19: k-mer entropy is bounded by log₂ of the number of distinct k-mers (uniform distribution maximises
    /// Shannon entropy) and equals it exactly when every k-mer is distinct.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property KmerEntropy_AtMostLog2DistinctKmers()
    {
        var gen = from s in RepeatRichDna(0, 60)
                  from k in Gen.Choose(1, 6)
                  select (s, k);
        return Prop.ForAll(gen.ToArbitrary(), x =>
        {
            double h = SequenceComplexity.CalculateKmerEntropy(x.s, x.k);
            int total = Math.Max(0, x.s.Length - x.k + 1);
            int distinct = Enumerable.Range(0, total).Select(i => x.s.Substring(i, x.k)).Distinct().Count();
            double bound = distinct == 0 ? 0 : Math.Log2(distinct);
            bool ok = h >= 0 && h <= bound + 1e-12 && (distinct != total || Math.Abs(h - bound) <= 1e-12);
            return ok.Label($"H_{x.k}('{x.s}') = {h:R}, log2(distinct={distinct}) = {bound:R}");
        });
    }

    #endregion

}
