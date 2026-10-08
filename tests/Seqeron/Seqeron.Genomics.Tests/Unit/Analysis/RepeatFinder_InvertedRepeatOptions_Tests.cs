namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// REP-INV-001 options and the scored variant:
/// <list type="bullet">
/// <item><c>maxMismatches</c> — EMBOSS <c>palindrome -nummismatches k</c> (palindrome.c: inward walk from a pairing
/// outer pair until the (k+1)-th mismatch, trailing mismatches trimmed, <c>palindrome_AInB</c> nesting filter);</item>
/// <item><c>maxArmLength</c> — EMBOSS <c>palindrome -maxpallen</c> (start pairs limited to span
/// 2·maxpallen + gaplimit + 1; longer stems suppress their sub-stems but are not printed);</item>
/// <item><c>allowWobble</c> — G·U (G·T) wobble pairs (Crick 1966; Varani &amp; McClain 2000) via the canonical
/// <see cref="RnaSecondaryStructure.CanPair"/>;</item>
/// <item><see cref="RepeatFinder.FindInvertedRepeatsScored(string,int,int,int,int,int)"/> — EMBOSS <c>einverted</c>
/// (Durbin &amp; Thierry-Mieg 1993).</item>
/// </list>
/// Values marked "EMBOSS" were produced by the EMBOSS 6.6.0 binaries (<c>palindrome -overlap Y</c>, <c>einverted</c>);
/// coordinates converted to 0-based. Evidence: docs/Evidence/REP-INV-001-Evidence.md.
/// </summary>
[TestFixture]
public class RepeatFinder_InvertedRepeatOptions_Tests
{
    private static List<(int Left, int Right, int Arm, int Mismatches)> Stems(IEnumerable<InvertedRepeatResult> results) =>
        results.Select(r => (r.LeftArmStart, r.RightArmStart, r.ArmLength, r.Mismatches)).ToList();

    #region maxMismatches (EMBOSS palindrome -nummismatches)

    /// <summary>
    /// EMBOSS palindrome -minpallen 4 -gaplimit 10 -maxpallen 11 on GAATTCAGG·AAAA·CCTCAATTC (right arm has one
    /// mismatch): k = 0 → nothing; k = 1 → gaattcagg/cttaactcc (1..9 / 22..14) with one interior mismatch;
    /// k = 2 → five stems.
    /// </summary>
    [Test]
    public void MaxMismatches_EmbossPalindromeWorkedExample()
    {
        const string seq = "GAATTCAGGAAAACCTCAATTC";

        Assert.Multiple(() =>
        {
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxMismatches: 0)), Is.Empty);
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxMismatches: 1)),
                Is.EqualTo(new[] { (0, 13, 9, 1) }));
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxMismatches: 2)),
                Is.EqualTo(new[] { (0, 8, 6, 2), (0, 10, 5, 2), (0, 12, 5, 2), (0, 13, 9, 1), (12, 17, 4, 2) }));
        });
    }

    /// <summary>
    /// EMBOSS palindrome -minpallen 4 -maxpallen 9 -gaplimit 10 -nummismatches 1 on GGGGAGAAAATTCTCCCC:
    /// ggggagaa/cccctctt (exact; the inner a·a mismatch is trimmed) and gaaaa/ctctt with one interior mismatch.
    /// </summary>
    [Test]
    public void MaxMismatches_TrailingMismatchTrimmed_InteriorCounted()
    {
        var stems = RepeatFinder.FindInvertedRepeats("GGGGAGAAAATTCTCCCC", 4, 10, 0, maxMismatches: 1).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(Stems(stems), Is.EqualTo(new[] { (0, 10, 8, 0), (5, 10, 5, 1) }));
            Assert.That(stems[1].LeftArm, Is.EqualTo("GAAAA"));
            Assert.That(stems[1].RightArm, Is.EqualTo("TTCTC"));
            Assert.That(stems[1].Loop, Is.Empty);
        });
    }

    /// <summary>Every reported mismatch stem pairs at both ends and has at most k interior mismatches.</summary>
    [Test]
    public void MaxMismatches_EndsPair_MismatchCountBounded()
    {
        var rnd = new Random(20260930);
        for (int t = 0; t < 200; t++)
        {
            string s = RandomSequence(rnd, rnd.Next(10, 120), "ACGT");
            int k = rnd.Next(1, 4);
            foreach (var r in RepeatFinder.FindInvertedRepeats(s, 3, 12, 0, maxMismatches: k))
            {
                int observed = Enumerable.Range(0, r.ArmLength)
                    .Count(i => !WatsonCrick(s[r.LeftArmStart + i], s[r.RightArmStart + r.ArmLength - 1 - i]));
                Assert.That(WatsonCrick(s[r.LeftArmStart], s[r.RightArmStart + r.ArmLength - 1]), s);
                Assert.That(WatsonCrick(s[r.LeftArmStart + r.ArmLength - 1], s[r.RightArmStart]), s);
                Assert.That(r.Mismatches, Is.EqualTo(observed).And.LessThanOrEqualTo(k), s);
            }
        }
    }

    /// <summary>maxMismatches = 0 (explicit) is the default exact search.</summary>
    [Test]
    public void MaxMismatches_Zero_IdenticalToDefault()
    {
        var rnd = new Random(7);
        for (int t = 0; t < 100; t++)
        {
            string s = RandomSequence(rnd, rnd.Next(0, 200), "ACGTN");
            Assert.That(RepeatFinder.FindInvertedRepeats(s, 3, 20, 2, maxMismatches: 0, maxArmLength: int.MaxValue, allowWobble: false),
                Is.EqualTo(RepeatFinder.FindInvertedRepeats(s, 3, 20, 2)), s);
        }
    }

    #endregion

    #region maxArmLength (EMBOSS palindrome -maxpallen)

    /// <summary>
    /// EMBOSS palindrome -minpallen 4 -gaplimit 10: a perfect arm-10 stem G×10·AAA·C×10 is reported with
    /// -maxpallen 11 (= len/2, unbounded) and not at all with -maxpallen 4, 5 or 6 — longer stems are not split.
    /// </summary>
    [Test]
    public void MaxArmLength_LongerStemNotSplit()
    {
        const string seq = "GGGGGGGGGGAAACCCCCCCCCC";

        Assert.Multiple(() =>
        {
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0)), Is.EqualTo(new[] { (0, 13, 10, 0) }));
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxArmLength: 11)), Is.EqualTo(new[] { (0, 13, 10, 0) }));
            foreach (int maxArm in new[] { 4, 5, 6 })
                Assert.That(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxArmLength: maxArm), Is.Empty, $"maxArm {maxArm}");
        });
    }

    /// <summary>
    /// EMBOSS palindrome -minpallen 4 -gaplimit 10 on TTGCATGCAAAAAATTTTTTTGCATGCAA: unbounded (-maxpallen 14) →
    /// (0,5,5) (0,15,14) (8,14,6) (19,24,5); -maxpallen 6 → the arm-14 stem disappears, the others remain;
    /// -maxpallen 5 → only the two arm-5 stems; -maxpallen 4 → nothing.
    /// </summary>
    [Test]
    public void MaxArmLength_EmbossWorkedExample()
    {
        const string seq = "TTGCATGCAAAAAATTTTTTTGCATGCAA";

        Assert.Multiple(() =>
        {
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxArmLength: 14)),
                Is.EqualTo(new[] { (0, 5, 5, 0), (0, 15, 14, 0), (8, 14, 6, 0), (19, 24, 5, 0) }));
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxArmLength: 6)),
                Is.EqualTo(new[] { (0, 5, 5, 0), (8, 14, 6, 0), (19, 24, 5, 0) }));
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxArmLength: 5)),
                Is.EqualTo(new[] { (0, 5, 5, 0), (19, 24, 5, 0) }));
            Assert.That(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 0, maxArmLength: 4), Is.Empty);
        });
    }

    #endregion

    #region allowWobble (G·U)

    /// <summary>
    /// GG(T|U)GC·AAAA·GCA(T|U)C: the stem G-C, G·U, U-A, G-C, C-G needs the G·U wobble; without it only a 3-bp
    /// stem remains (below minArm 4). Wobble pairing follows <see cref="RnaSecondaryStructure.CanPair"/> (T read as U).
    /// </summary>
    [TestCase("GGTGCAAAAGCATC")]
    [TestCase("GGUGCAAAAGCAUC")]
    [TestCase("ggugcaaaagcauc")]
    public void AllowWobble_GuPairCompletesStem(string seq)
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 3), Is.Empty);
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(seq, 4, 10, 3, allowWobble: true)),
                Is.EqualTo(new[] { (0, 9, 5, 0) }));
        });
    }

    /// <summary>G/T arms pair only through wobble: G×6·NNNN·T×6 is one maximal arm-6 stem (N never pairs).</summary>
    [Test]
    public void AllowWobble_GtOnlyStem()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindInvertedRepeats("GGGGGGNNNNTTTTTT", 4, 10, 3), Is.Empty);
            Assert.That(Stems(RepeatFinder.FindInvertedRepeats("GGGGGGNNNNTTTTTT", 4, 10, 3, allowWobble: true)),
                Is.EqualTo(new[] { (0, 10, 6, 0) }));
        });
    }

    #endregion

    #region Differential: literal palindrome.c oracle (all options) and all-stems brute force (wobble)

    /// <summary>
    /// Literal transcription of the EMBOSS palindrome.c search (current/rev loops, inward walk with mismatch
    /// counter and trailing-mismatch trim, <c>palindrome_AInB</c> list filter, <c>-maxpallen</c> start-span bound and
    /// print filter), generalised by an inner stop at loop &lt; minLoop and a pairing predicate; compared on random
    /// sequences over every option combination.
    /// </summary>
    [Test]
    public void AllOptions_MatchLiteralPalindromeOracle()
    {
        var rnd = new Random(42);
        string[] alphabets = { "ACGT", "ACGT", "GT", "ACGTN", "ACGU", "AGTU", "ACGTRY-" };
        for (int t = 0; t < 400; t++)
        {
            string s = RandomSequence(rnd, rnd.Next(0, 45), alphabets[rnd.Next(alphabets.Length)]);
            int minArm = rnd.Next(2, 6), minLoop = rnd.Next(0, 7), maxLoop = minLoop + rnd.Next(0, 16);
            int k = rnd.Next(0, 4);
            int maxArm = rnd.Next(3) == 0 ? rnd.Next(minArm, minArm + 7) : int.MaxValue;
            bool wobble = rnd.Next(2) == 0;

            var expected = LiteralPalindrome(s.ToUpperInvariant(), minArm, maxLoop, minLoop, k, maxArm, wobble);
            var actual = Stems(RepeatFinder.FindInvertedRepeats(s, minArm, maxLoop, minLoop, k, maxArm, wobble));
            Assert.That(actual, Is.EqualTo(expected), $"{s} {minArm} {maxLoop} {minLoop} k={k} maxArm={maxArm} wobble={wobble}");
        }
    }

    /// <summary>Exact wobble stems = every exact stem under the RNA pair set, minus those nested in both arms of another.</summary>
    [Test]
    public void AllowWobble_MatchesAllStemsBruteForce()
    {
        var rnd = new Random(99);
        for (int t = 0; t < 300; t++)
        {
            string s = RandomSequence(rnd, rnd.Next(0, 40), rnd.Next(2) == 0 ? "ACGT" : "GTAU");
            int minArm = rnd.Next(2, 5), minLoop = rnd.Next(0, 5), maxLoop = minLoop + rnd.Next(0, 12);

            var all = new List<(int I, int J, int A)>();
            for (int i = 0; i < s.Length; i++)
                for (int a = minArm; i + 2 * a + minLoop <= s.Length; a++)
                    for (int loop = minLoop; loop <= maxLoop && i + 2 * a + loop <= s.Length; loop++)
                    {
                        int j = i + a + loop;
                        if (Enumerable.Range(0, a).All(x => RnaPair(s[i + x], s[j + a - 1 - x])))
                            all.Add((i, j, a));
                    }
            var expected = all
                .Where(x => !all.Any(y => y != x && y.I <= x.I && y.I + y.A >= x.I + x.A && y.J <= x.J && y.J + y.A >= x.J + x.A))
                .OrderBy(x => x.I).ThenBy(x => x.J)
                .Select(x => (x.I, x.J, x.A, 0)).ToList();

            Assert.That(Stems(RepeatFinder.FindInvertedRepeats(s, minArm, maxLoop, minLoop, allowWobble: true)),
                Is.EqualTo(expected), s);
        }
    }

    internal static List<(int, int, int, int)> LiteralPalindrome(
        string s, int minLen, int maxGap, int minLoop, int maxMismatches, int maxLen, bool wobble)
    {
        int n = s.Length;
        var list = new List<(int Fs, int Fe, int Rs, int Re)>();
        for (int current = 0; current < n - 1; current++)
        {
            long iend = Math.Min((long)current + 2L * maxLen + maxGap, n - 1);
            for (int rev = (int)iend; rev > current + minLen; rev--)
            {
                if (!Pair(s[current], s[rev], wobble))
                    continue;
                int count = 0, mismatches = 0, mismatchAtEnd = 0, ic = current, ir = rev;
                while (mismatches <= maxMismatches && ir - ic - 1 >= minLoop)
                {
                    if (Pair(s[ic], s[ir], wobble)) mismatchAtEnd = 0;
                    else { mismatches++; mismatchAtEnd++; }
                    count++; ic++; ir--;
                }
                count -= mismatchAtEnd;
                int gap = rev - current - 2 * count + 1;
                if (count < minLen || gap > maxGap)
                    continue;
                var a = (Fs: current, Fe: current + count, Rs: rev, Re: rev - count);
                if (!list.Any(b => a.Fs >= b.Fs && a.Fe <= b.Fe && a.Rs <= b.Rs && a.Re >= b.Re))
                    list.Add(a);
            }
        }
        return list
            .Where(p => p.Fe - p.Fs <= maxLen)
            .Select(p => (p.Fs, p.Re + 1, p.Fe - p.Fs,
                Enumerable.Range(0, p.Fe - p.Fs).Count(x => !Pair(s[p.Fs + x], s[p.Rs - x], wobble))))
            .OrderBy(p => p.Item1).ThenBy(p => p.Item2).ToList();
    }

    private static bool Pair(char a, char b, bool wobble) => wobble ? RnaPair(a, b) : WatsonCrick(a, b);

    private static bool WatsonCrick(char a, char b) =>
        (char.ToUpperInvariant(a), char.ToUpperInvariant(b)) is ('A', 'T') or ('T', 'A') or ('C', 'G') or ('G', 'C');

    private static bool RnaPair(char a, char b)
    {
        char x = char.ToUpperInvariant(a) == 'U' ? 'T' : char.ToUpperInvariant(a);
        char y = char.ToUpperInvariant(b) == 'U' ? 'T' : char.ToUpperInvariant(b);
        return (x, y) is ('A', 'T') or ('T', 'A') or ('C', 'G') or ('G', 'C') or ('G', 'T') or ('T', 'G');
    }

    private static string RandomSequence(Random rnd, int n, string alphabet) =>
        new(Enumerable.Range(0, n).Select(_ => alphabet[rnd.Next(alphabet.Length)]).ToArray());

    #endregion

    #region Validation

    [Test]
    public void Options_InvalidParameters_ThrowEagerly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeats("ACGT", maxMismatches: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeats("ACGT", minArmLength: 5, maxArmLength: 4));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeats((string)null!, maxMismatches: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeats(new DnaSequence("ACGT"), maxArmLength: 1));
            Assert.That(RepeatFinder.FindInvertedRepeats("", maxMismatches: 3, maxArmLength: 4, allowWobble: true), Is.Empty);
        });
    }

    #endregion

    #region FindInvertedRepeatsScored (EMBOSS einverted)

    private const string PlantedImperfect =
        "CCCAACCCATGCGTACGTTAGCCTAGGATCCATTTTTTTTTGGATACTAGGCAACGTACGCATGGGAAGGG";

    /// <summary>
    /// EMBOSS einverted (defaults gap 12, threshold 50, match 3, mismatch −4) on a 24-bp arm planted with one mismatch
    /// and one deleted base in the right arm: "Score 60: 28/31 (90%) matches, 1 gaps", left 1..32, right 71..41
    /// (the flanks CCC·AACCC / GGGAAGGG extend the alignment); 60 = 28·3 − 3·4 − 12.
    /// </summary>
    [Test]
    public void Scored_EinvertedWorkedExample_Defaults()
    {
        var r = RepeatFinder.FindInvertedRepeatsScored(PlantedImperfect).Single();

        Assert.Multiple(() =>
        {
            Assert.That((r.LeftArmStart, r.LeftArmEnd, r.RightArmStart, r.RightArmEnd), Is.EqualTo((0, 31, 40, 70)));
            Assert.That((r.Score, r.Matches, r.Mismatches, r.Gaps), Is.EqualTo((60, 28, 3, 1)));
            Assert.That(r.LeftArmAlignment, Is.EqualTo("CCCAACCCATGCGTACGTTAGCCTAGGATCCA"));
            Assert.That(r.MatchLine, Is.EqualTo("|||  |||||||||||||| |||||| |||||"));
            Assert.That(r.RightArmAlignment, Is.EqualTo("GGGAAGGGTACGCATGCAA-CGGATCATAGGT"));
            Assert.That((r.LeftArmLength, r.RightArmLength, r.LoopLength), Is.EqualTo((32, 31, 8)));
            Assert.That((int)r.PercentMatches, Is.EqualTo(90));
        });
    }

    /// <summary>EMBOSS einverted -gap 8 -threshold 30 on the same sequence: Score 64 (= 84 − 12 − 8), same alignment.</summary>
    [Test]
    public void Scored_EinvertedWorkedExample_CustomGap()
    {
        var r = RepeatFinder.FindInvertedRepeatsScored(PlantedImperfect, gapPenalty: 8, threshold: 30).Single();

        Assert.That((r.Score, r.Matches, r.Mismatches, r.Gaps, r.LeftArmStart, r.RightArmEnd),
            Is.EqualTo((64, 28, 3, 1, 0, 70)));
    }

    /// <summary>
    /// EMBOSS einverted -threshold 20 on T×8(GCAT)×4GCAA·A×12·TT(GCAT)×4GCA×8 … : "Score 84: 28/28 (100%) matches,
    /// 0 gaps", 1..28 / 68..41.
    /// </summary>
    [Test]
    public void Scored_EinvertedPerfectRepeat()
    {
        const string seq = "TTTTTTTTGCATGCATGCATGCATGCAAAAAAAAAAAAAATTGCATGCATGCATGCATGCAAAAAAAA";

        var r = RepeatFinder.FindInvertedRepeatsScored(seq, threshold: 20).Single();

        Assert.Multiple(() =>
        {
            Assert.That((r.LeftArmStart, r.LeftArmEnd, r.RightArmStart, r.RightArmEnd), Is.EqualTo((0, 27, 40, 67)));
            Assert.That((r.Score, r.Matches, r.Mismatches, r.Gaps), Is.EqualTo((84, 28, 0, 0)));
            Assert.That(r.MatchLine, Is.EqualTo(new string('|', 28)));
            Assert.That(r.RightArmAlignment, Is.EqualTo("AAAAAAAACGTACGTACGTACGTACGTT"));
        });
    }

    /// <summary>EMBOSS einverted defaults (einverted.acd, EMBOSS 6.6.0).</summary>
    [Test]
    public void Scored_DefaultsAreEinverted()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.EinvertedDefaultGapPenalty, Is.EqualTo(12));
            Assert.That(RepeatFinder.EinvertedDefaultThreshold, Is.EqualTo(50));
            Assert.That(RepeatFinder.EinvertedDefaultMatchScore, Is.EqualTo(3));
            Assert.That(RepeatFinder.EinvertedDefaultMismatchScore, Is.EqualTo(-4));
            Assert.That(RepeatFinder.EinvertedDefaultMaxRepeatLength, Is.EqualTo(2000));
        });
    }

    /// <summary>Only a/c/g/t pair in einverted (N scores as a mismatch); input case is irrelevant.</summary>
    [Test]
    public void Scored_NeverPairsN_CaseInsensitive()
    {
        Assert.Multiple(() =>
        {
            Assert.That(RepeatFinder.FindInvertedRepeatsScored(new string('N', 300), threshold: 0), Is.Empty);
            Assert.That(RepeatFinder.FindInvertedRepeatsScored(PlantedImperfect.ToLowerInvariant()),
                Is.EqualTo(RepeatFinder.FindInvertedRepeatsScored(PlantedImperfect)));
            Assert.That(RepeatFinder.FindInvertedRepeatsScored(new DnaSequence(PlantedImperfect)),
                Is.EqualTo(RepeatFinder.FindInvertedRepeatsScored(PlantedImperfect)));
        });
    }

    /// <summary>Structural invariants of every scored result on random sequences with planted imperfect repeats.</summary>
    [Test]
    public void Scored_ResultInvariants()
    {
        var rnd = new Random(5);
        for (int t = 0; t < 60; t++)
        {
            string s = PlantRepeat(rnd, RandomSequence(rnd, rnd.Next(50, 400), "ACGT"));
            foreach (var r in RepeatFinder.FindInvertedRepeatsScored(s, threshold: 30, maxRepeatLength: rnd.Next(20, 300)))
            {
                Assert.That(r.LeftArmStart, Is.LessThanOrEqualTo(r.LeftArmEnd), s);
                Assert.That(r.LeftArmEnd, Is.LessThan(r.RightArmStart), s);
                Assert.That(r.RightArmStart, Is.LessThanOrEqualTo(r.RightArmEnd), s);
                Assert.That(r.Score, Is.GreaterThanOrEqualTo(30), s);
                Assert.That(r.MatchLine.Length, Is.EqualTo(r.LeftArmAlignment.Length).And.EqualTo(r.RightArmAlignment.Length));
                Assert.That(r.LeftArmAlignment.Count(c => c != '-'), Is.EqualTo(r.LeftArmLength), s);
                Assert.That(r.RightArmAlignment.Count(c => c != '-'), Is.EqualTo(r.RightArmLength), s);
                // einverted quirks (reproduced): a final trace-back gap step is counted but not drawn (its column then
                // shows both bases, possibly '|'), and an innermost loop-0 cell ("blunt join") is scored but not counted.
                Assert.That(r.MatchLine.Count(c => c == '|') - r.Matches, Is.InRange(0, 1), s);
                Assert.That(r.Gaps - r.LeftArmAlignment.Count(c => c == '-') - r.RightArmAlignment.Count(c => c == '-'), Is.InRange(0, 1), s);
                Assert.That(r.Score, Is.GreaterThanOrEqualTo(
                    RepeatFinder.EinvertedDefaultMatchScore * r.Matches + RepeatFinder.EinvertedDefaultMismatchScore * r.Mismatches
                    - RepeatFinder.EinvertedDefaultGapPenalty * r.Gaps), s);
            }
        }
    }

    [Test]
    public void Scored_InvalidParameters_ThrowEagerly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", gapPenalty: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", threshold: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", matchScore: -1));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", mismatchScore: 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored("ACGT", maxRepeatLength: 1));
            Assert.Throws<ArgumentOutOfRangeException>(() => RepeatFinder.FindInvertedRepeatsScored((string)null!, threshold: 1_000_000));
            Assert.Throws<ArgumentNullException>(() => RepeatFinder.FindInvertedRepeatsScored((DnaSequence)null!));
            Assert.That(RepeatFinder.FindInvertedRepeatsScored((string)null!), Is.Empty);
            Assert.That(RepeatFinder.FindInvertedRepeatsScored(""), Is.Empty);
        });
    }

    internal static string PlantRepeat(Random rnd, string background)
    {
        int arm = rnd.Next(12, 40), loop = rnd.Next(0, 20);
        if (2 * arm + loop + 4 > background.Length)
            return background;
        var left = RandomSequence(rnd, arm, "ACGT");
        var right = new System.Text.StringBuilder(left.Reverse().Select(c => c switch
        {
            'A' => 'T', 'T' => 'A', 'C' => 'G', _ => 'C'
        }).ToArray().AsSpan().ToString());
        for (int m = rnd.Next(0, 3); m > 0; m--)
            right[rnd.Next(right.Length)] = "ACGT"[rnd.Next(4)];
        if (rnd.Next(2) == 0)
            right.Remove(rnd.Next(right.Length), 1);
        string ir = left + RandomSequence(rnd, loop, "ACGT") + right;
        int p = rnd.Next(0, background.Length - ir.Length + 1);
        return background[..p] + ir + background[(p + ir.Length)..];
    }

    #endregion
}
