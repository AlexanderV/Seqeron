using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Tests.Unit.Analysis;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier properties for the REP-INV-001 options added by review batch B04 (F21–F24):
/// mismatch-tolerant stems and the maximum arm length (EMBOSS palindrome -nummismatches / -maxpallen),
/// G·U wobble, and the einverted-equivalent scored search.
/// Test Unit: REP-INV-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Repeats")]
public class RepInvOptionsProperties
{
    private static Gen<string> Dna(string alphabets, int minLen, int maxLen) =>
        from alpha in Gen.Elements(alphabets.Split('|'))
        from n in Gen.Choose(minLen, maxLen)
        from chars in Gen.Elements(alpha.ToCharArray()).ArrayOf(n)
        select new string(chars);

    /// <summary>
    /// For every option combination the result equals the literal palindrome.c transcription (inward walk,
    /// trailing-mismatch trim, AInB list, -maxpallen span bound and print filter; generalised to minLoop/wobble).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property AllOptions_EqualLiteralPalindromeOracle()
    {
        var gen = from s in Dna("ACGT|AT|GC|GT|ACGTN|ACGU", 0, 50)
                  from minArm in Gen.Choose(2, 5)
                  from minLoop in Gen.Choose(0, 5)
                  from extra in Gen.Choose(0, 15)
                  from k in Gen.Choose(0, 4)
                  from bounded in Gen.Elements(true, false)
                  from armExtra in Gen.Choose(0, 6)
                  from wobble in Gen.Elements(true, false)
                  select (s, minArm, minLoop, maxLoop: minLoop + extra, k, maxArm: bounded ? minArm + armExtra : int.MaxValue, wobble);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var expected = RepeatFinder_InvertedRepeatOptions_Tests.LiteralPalindrome(
                c.s, c.minArm, c.maxLoop, c.minLoop, c.k, c.maxArm, c.wobble);
            var actual = RepeatFinder.FindInvertedRepeats(c.s, c.minArm, c.maxLoop, c.minLoop, c.k, c.maxArm, c.wobble)
                .Select(r => (r.LeftArmStart, r.RightArmStart, r.ArmLength, r.Mismatches)).ToList();
            return actual.SequenceEqual(expected).Label($"{c}");
        });
    }

    /// <summary>No reported stem (any option) lies inside another reported stem in both arms; arm ≤ maxArmLength.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Options_NonNested_ArmBounded()
    {
        var gen = from s in Dna("ACGT|AT|GT", 0, 120)
                  from k in Gen.Choose(0, 3)
                  from maxArm in Gen.Choose(3, 12)
                  from wobble in Gen.Elements(true, false)
                  select (s, k, maxArm, wobble);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var r = RepeatFinder.FindInvertedRepeats(c.s, 3, 15, 0, c.k, c.maxArm, c.wobble).ToList();
            bool nested = r.Any(a => r.Any(b => !a.Equals(b)
                && b.LeftArmStart <= a.LeftArmStart && b.LeftArmStart + b.ArmLength >= a.LeftArmStart + a.ArmLength
                && b.RightArmStart <= a.RightArmStart && b.RightArmStart + b.ArmLength >= a.RightArmStart + a.ArmLength));
            return (!nested && r.All(x => x.ArmLength <= c.maxArm && x.Mismatches <= c.k)).Label($"{c}");
        });
    }

    /// <summary>
    /// Scored results: ordered arm coordinates, score ≥ threshold, alignment rows consistent with coordinates,
    /// and the einverted count quirks bounded (≤ 1 undrawn final gap, ≤ 1 uncounted '|').
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property Scored_StructuralInvariants()
    {
        var gen = from seed in Gen.Choose(0, int.MaxValue - 1)
                  from n in Gen.Choose(20, 300)
                  from threshold in Gen.Elements(20, 30, 50)
                  from gap in Gen.Elements(4, 8, 12)
                  from maxRep in Gen.Elements(10, 60, 2000)
                  select (seed, n, threshold, gap, maxRep);
        return Prop.ForAll(gen.ToArbitrary(), c =>
        {
            var rnd = new Random(c.seed);
            string bg = new(Enumerable.Range(0, c.n).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            string s = RepeatFinder_InvertedRepeatOptions_Tests.PlantRepeat(rnd, bg);
            var res = RepeatFinder.FindInvertedRepeatsScored(s, c.gap, c.threshold, maxRepeatLength: c.maxRep).ToList();
            bool ok = res.All(r =>
                r.LeftArmStart <= r.LeftArmEnd && r.LeftArmEnd < r.RightArmStart && r.RightArmStart <= r.RightArmEnd
                && r.Score >= c.threshold
                && r.LeftArmAlignment.Length == r.MatchLine.Length && r.RightArmAlignment.Length == r.MatchLine.Length
                && r.LeftArmAlignment.Count(ch => ch != '-') == r.LeftArmLength
                && r.RightArmAlignment.Count(ch => ch != '-') == r.RightArmLength
                && r.MatchLine.Count(ch => ch == '|') - r.Matches is 0 or 1
                && r.Gaps - r.LeftArmAlignment.Count(ch => ch == '-') - r.RightArmAlignment.Count(ch => ch == '-') is 0 or 1);
            return ok.Label($"{c} {s}");
        });
    }
}
