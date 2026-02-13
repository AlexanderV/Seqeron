using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for Edit Distance (Levenshtein distance).
/// Verifies metric space axioms that must hold for ALL string pairs.
///
/// Test Unit: PAT-APPROX-002 (Property Extension)
/// Evidence: Wikipedia — Levenshtein distance properties (metric axioms).
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Alignment")]
public class EditDistanceProperties
{
    private static Arbitrary<string> SafeStringArbitrary() =>
        ArbMap.Default.ArbFor<NonNull<string>>().Generator
            .Select(s => s.Get ?? "")
            .ToArbitrary();

    private static Arbitrary<string> DnaStringArbitrary() =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Select(a => new string(a))
            .ToArbitrary();

    /// <summary>
    /// Metric axiom: d(a, a) == 0 — distance from a string to itself is always zero.
    /// Evidence: Wikipedia — "the distance between two strings is zero if and only if the strings are equal."
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Identity_DistanceToSelfIsZero()
    {
        return Prop.ForAll(SafeStringArbitrary(), s =>
        {
            return (ApproximateMatcher.EditDistance(s, s) == 0)
                .Label($"d(\"{Truncate(s)}\", \"{Truncate(s)}\") should be 0");
        });
    }

    /// <summary>
    /// Metric axiom: d(a, b) == d(b, a) — edit distance is symmetric.
    /// Evidence: Wikipedia — "the Levenshtein distance is a metric."
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Symmetry_DistanceIsCommutative()
    {
        return Prop.ForAll(SafeStringArbitrary(), SafeStringArbitrary(), (a, b) =>
        {
            int dAB = ApproximateMatcher.EditDistance(a, b);
            int dBA = ApproximateMatcher.EditDistance(b, a);
            return (dAB == dBA)
                .Label($"d(a,b)={dAB}, d(b,a)={dBA}");
        });
    }

    /// <summary>
    /// Metric axiom: d(a, c) &lt;= d(a, b) + d(b, c) — triangle inequality.
    /// Evidence: Wikipedia — metric property of Levenshtein distance.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TriangleInequality()
    {
        return Prop.ForAll(SafeStringArbitrary(), SafeStringArbitrary(), SafeStringArbitrary(),
            (a, b, c) =>
            {
                int ab = ApproximateMatcher.EditDistance(a, b);
                int bc = ApproximateMatcher.EditDistance(b, c);
                int ac = ApproximateMatcher.EditDistance(a, c);
                return (ac <= ab + bc)
                    .Label($"d(a,c)={ac} <= d(a,b)={ab} + d(b,c)={bc}");
            });
    }

    /// <summary>
    /// Non-negativity: d(a, b) >= 0 — distance is never negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property NonNegativity()
    {
        return Prop.ForAll(SafeStringArbitrary(), SafeStringArbitrary(), (a, b) =>
        {
            return (ApproximateMatcher.EditDistance(a, b) >= 0)
                .Label("Distance must be non-negative");
        });
    }

    /// <summary>
    /// Upper bound: d(a, b) &lt;= max(|a|, |b|) — can always transform by deleting all + inserting all.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property UpperBound_MaxLength()
    {
        return Prop.ForAll(SafeStringArbitrary(), SafeStringArbitrary(), (a, b) =>
        {
            int d = ApproximateMatcher.EditDistance(a, b);
            int maxLen = Math.Max(a.Length, b.Length);
            return (d <= maxLen)
                .Label($"d={d}, max(|a|,|b|)={maxLen}");
        });
    }

    /// <summary>
    /// Empty string property: d("", s) == |s| — transforming empty to s requires |s| insertions.
    /// Evidence: Wikipedia definition — base case of the recursive formula.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property EmptyString_DistanceEqualsLength()
    {
        return Prop.ForAll(SafeStringArbitrary(), s =>
        {
            int d = ApproximateMatcher.EditDistance("", s);
            return (d == s.Length)
                .Label($"d(\"\", \"{Truncate(s)}\")={d}, expected={s.Length}");
        });
    }

    /// <summary>
    /// DNA-specific: edit distance between a DNA sequence and its reverse is bounded by length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DnaSequence_DistanceToReverseIsBounded()
    {
        return Prop.ForAll(DnaStringArbitrary(), seq =>
        {
            if (seq.Length == 0) return true.ToProperty();

            string reversed = new string(seq.Reverse().ToArray());
            int d = ApproximateMatcher.EditDistance(seq, reversed);
            return (d <= seq.Length)
                .Label($"d(seq, rev)={d}, |seq|={seq.Length}");
        });
    }

    private static string Truncate(string s, int maxLen = 20) =>
        s.Length <= maxLen ? s : s[..maxLen] + "...";
}
