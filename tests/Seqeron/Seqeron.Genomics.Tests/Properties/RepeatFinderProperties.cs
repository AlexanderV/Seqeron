using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for repeat finder algorithms.
/// Verifies invariants for palindromes, inverted repeats, and direct repeats.
///
/// Test Units: REP-PALIN-001, REP-INV-001, REP-DIRECT-001, REP-STR-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class RepeatFinderProperties
{
    /// <summary>
    /// Palindromes found must be within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Palindrome_Positions_WithinBounds()
    {
        string seq = "ACGTGAATTCACGTACGTGATATCACGT";
        var palindromes = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();

        foreach (var p in palindromes)
        {
            Assert.That(p.Position, Is.GreaterThanOrEqualTo(0),
                $"Palindrome position {p.Position} must be ≥ 0");
            Assert.That(p.Position + p.Length, Is.LessThanOrEqualTo(seq.Length),
                $"Palindrome at {p.Position}+{p.Length} exceeds sequence length {seq.Length}");
        }
    }

    /// <summary>
    /// Each palindrome sequence equals its own reverse complement.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Palindrome_EqualsReverseComplement()
    {
        string seq = "ACGTGAATTCACGTACGTGATATCACGT";
        var palindromes = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();

        foreach (var p in palindromes)
        {
            string palSeq = seq.Substring(p.Position, p.Length);
            string revComp = ReverseComplement(palSeq);
            Assert.That(palSeq, Is.EqualTo(revComp),
                $"Palindrome '{palSeq}' at {p.Position} must equal its reverse complement '{revComp}'");
        }
    }

    /// <summary>
    /// Palindrome lengths are even (DNA palindromes pair symmetrically).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Palindrome_LengthsAreEven()
    {
        string seq = "ACGTGAATTCACGTACGTGATATCACGT";
        var palindromes = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();

        foreach (var p in palindromes)
            Assert.That(p.Length % 2, Is.EqualTo(0),
                $"Palindrome length {p.Length} must be even");
    }

    /// <summary>
    /// Inverted repeat positions are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void InvertedRepeat_Positions_WithinBounds()
    {
        string seq = "ACGTACGTAAAAAAATGCATGCA";
        var repeats = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4).ToList();

        foreach (var r in repeats)
        {
            Assert.That(r.LeftArmStart, Is.GreaterThanOrEqualTo(0));
            Assert.That(r.RightArmStart + r.ArmLength, Is.LessThanOrEqualTo(seq.Length));
        }
    }

    /// <summary>
    /// Direct repeat positions are within sequence bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DirectRepeat_Positions_WithinBounds()
    {
        string seq = "ACGTACGTTTTTTTTTACGTACGT";
        var repeats = RepeatFinder.FindDirectRepeats(seq, minLength: 4).ToList();

        foreach (var r in repeats)
        {
            Assert.That(r.FirstPosition, Is.GreaterThanOrEqualTo(0));
            Assert.That(r.SecondPosition, Is.GreaterThanOrEqualTo(0));
            Assert.That(r.FirstPosition, Is.LessThan(seq.Length));
            Assert.That(r.SecondPosition, Is.LessThan(seq.Length));
        }
    }

    /// <summary>
    /// Microsatellite repeat unit length is within requested range.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Microsatellite_UnitLength_InRange()
    {
        string seq = "ACACACACACACGTGTGTGTGTAAAAAAAAAA";
        int minUnit = 1, maxUnit = 4;
        var results = RepeatFinder.FindMicrosatellites(seq, minUnit, maxUnit, minRepeats: 3).ToList();

        foreach (var r in results)
            Assert.That(r.RepeatUnit.Length, Is.InRange(minUnit, maxUnit),
                $"Unit length {r.RepeatUnit.Length} must be in [{minUnit}, {maxUnit}]");
    }

    /// <summary>
    /// Microsatellite repeat count is at least the requested minimum.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Microsatellite_RepeatCount_MeetsMinimum()
    {
        string seq = "ACACACACACACGTGTGTGTGTAAAAAAAAAA";
        int minRepeats = 3;
        var results = RepeatFinder.FindMicrosatellites(seq, 1, 4, minRepeats).ToList();

        foreach (var r in results)
            Assert.That(r.RepeatCount, Is.GreaterThanOrEqualTo(minRepeats),
                $"Repeat count {r.RepeatCount} must be ≥ {minRepeats}");
    }

    private static string ReverseComplement(string dna)
    {
        var comp = dna.Select(c => c switch
        {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            _ => c
        }).Reverse().ToArray();
        return new string(comp);
    }
}
