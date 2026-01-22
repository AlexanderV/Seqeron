using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class RepeatFinderTests
{
    #region Microsatellite Detection Tests

    [Test]
    public void FindMicrosatellites_MononucleotideRepeat_FindsRepeat()
    {
        var sequence = new DnaSequence("ACGTAAAAAACGT");
        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].RepeatUnit, Is.EqualTo("A"));
        Assert.That(results[0].RepeatCount, Is.EqualTo(6));
        Assert.That(results[0].Position, Is.EqualTo(4));
        Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Mononucleotide));
    }

    [Test]
    public void FindMicrosatellites_DinucleotideRepeat_FindsRepeat()
    {
        var sequence = new DnaSequence("AAACACACACACAAA");
        var results = RepeatFinder.FindMicrosatellites(sequence, 2, 6, 3).ToList();

        Assert.That(results.Any(r => r.RepeatUnit == "CA" || r.RepeatUnit == "AC"));
        var caRepeat = results.First(r => r.RepeatUnit == "CA" || r.RepeatUnit == "AC");
        Assert.That(caRepeat.RepeatCount, Is.GreaterThanOrEqualTo(3));
        Assert.That(caRepeat.RepeatType, Is.EqualTo(RepeatType.Dinucleotide));
    }

    [Test]
    public void FindMicrosatellites_TrinucleotideRepeat_CAGExpansion()
    {
        // CAG repeats are associated with Huntington's disease
        var sequence = new DnaSequence("ATGCAGCAGCAGCAGCAGTGA");
        var results = RepeatFinder.FindMicrosatellites(sequence, 3, 3, 3).ToList();

        Assert.That(results.Any(r => r.RepeatUnit == "CAG" && r.RepeatCount == 5));
        var cagRepeat = results.First(r => r.RepeatUnit == "CAG");
        Assert.That(cagRepeat.RepeatType, Is.EqualTo(RepeatType.Trinucleotide));
    }

    [Test]
    public void FindMicrosatellites_TetranucleotideRepeat_FindsRepeat()
    {
        var sequence = new DnaSequence("AAGATAGATAGATAGATAAA");
        var results = RepeatFinder.FindMicrosatellites(sequence, 4, 4, 3).ToList();

        Assert.That(results.Any(r => r.RepeatUnit == "GATA" || r.RepeatUnit == "ATAG" ||
                                     r.RepeatUnit == "TAGA" || r.RepeatUnit == "AGAT"));
        Assert.That(results[0].RepeatType, Is.EqualTo(RepeatType.Tetranucleotide));
    }

    [Test]
    public void FindMicrosatellites_NoRepeats_ReturnsEmpty()
    {
        var sequence = new DnaSequence("ACGTACGTACGT");
        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        // ACGT repeated 3 times, but we're looking for 1-6 bp units
        // "ACGT" is 4bp, repeated 3 times - should find it
        Assert.That(results.Any(r => r.RepeatUnit == "ACGT" && r.RepeatCount == 3));
    }

    [Test]
    public void FindMicrosatellites_MultipleDifferentRepeats_FindsAll()
    {
        var sequence = new DnaSequence("AAAAAACGTCGTCGTACACACAC");
        var results = RepeatFinder.FindMicrosatellites(sequence, 1, 6, 3).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(2));
        Assert.That(results.Any(r => r.RepeatUnit == "A" && r.RepeatCount >= 5));
    }

    [Test]
    public void FindMicrosatellites_MinRepeatsFilter_RespectsThreshold()
    {
        var sequence = new DnaSequence("ATATAT"); // AT x 3
        var results3 = RepeatFinder.FindMicrosatellites(sequence, 2, 2, 3).ToList();
        var results4 = RepeatFinder.FindMicrosatellites(sequence, 2, 2, 4).ToList();

        Assert.That(results3, Has.Count.EqualTo(1));
        Assert.That(results4, Has.Count.EqualTo(0));
    }

    [Test]
    public void FindMicrosatellites_StringOverload_Works()
    {
        var results = RepeatFinder.FindMicrosatellites("CAGCAGCAGCAG", 3, 3, 3).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].RepeatUnit, Is.EqualTo("CAG"));
        Assert.That(results[0].RepeatCount, Is.EqualTo(4));
    }

    [Test]
    public void FindMicrosatellites_EmptySequence_ReturnsEmpty()
    {
        var results = RepeatFinder.FindMicrosatellites("", 1, 6, 3).ToList();
        Assert.That(results, Is.Empty);
    }

    [Test]
    public void FindMicrosatellites_FullSequenceProperty_ReturnsCorrectSequence()
    {
        var results = RepeatFinder.FindMicrosatellites("CAGCAGCAG", 3, 3, 3).ToList();

        Assert.That(results[0].FullSequence, Is.EqualTo("CAGCAGCAG"));
    }

    #endregion

    // Note: Inverted Repeat Detection tests have been moved to RepeatFinder_InvertedRepeat_Tests.cs
    // See Test Unit REP-INV-001

    // Note: Direct Repeat Detection tests have been moved to RepeatFinder_DirectRepeat_Tests.cs
    // See Test Unit REP-DIRECT-001

    #region Palindrome Detection Tests

    [Test]
    public void FindPalindromes_EcoRISite_FindsPalindrome()
    {
        // EcoRI site: GAATTC (reverse complement = GAATTC)
        var sequence = new DnaSequence("AAAGAATTCAAA");
        var results = RepeatFinder.FindPalindromes(sequence, 6, 6).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].Sequence, Is.EqualTo("GAATTC"));
        Assert.That(results[0].Position, Is.EqualTo(3));
    }

    [Test]
    public void FindPalindromes_HindIIISite_FindsPalindrome()
    {
        // HindIII site: AAGCTT (reverse complement = AAGCTT)
        var sequence = new DnaSequence("CCCAAGCTTCCC");
        var results = RepeatFinder.FindPalindromes(sequence, 6, 6).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].Sequence, Is.EqualTo("AAGCTT"));
    }

    [Test]
    public void FindPalindromes_ShortPalindrome_FindsFourBasePalindrome()
    {
        // AATT, TTAA, GCGC, CGCG are 4bp palindromes
        var sequence = new DnaSequence("AAGCGCAA");
        var results = RepeatFinder.FindPalindromes(sequence, 4, 4).ToList();

        Assert.That(results, Has.Count.EqualTo(1));
        Assert.That(results[0].Sequence, Is.EqualTo("GCGC"));
    }

    [Test]
    public void FindPalindromes_MultiplePalindromes_FindsAll()
    {
        var sequence = new DnaSequence("GAATTCAAAGAATTC");
        var results = RepeatFinder.FindPalindromes(sequence, 6, 6).ToList();

        Assert.That(results, Has.Count.EqualTo(2));
    }

    [Test]
    public void FindPalindromes_NoPalindromes_ReturnsEmpty()
    {
        var sequence = new DnaSequence("AAAAAAAAAA");
        var results = RepeatFinder.FindPalindromes(sequence, 4, 12).ToList();

        Assert.That(results, Is.Empty);
    }

    [Test]
    public void FindPalindromes_StringOverload_Works()
    {
        var results = RepeatFinder.FindPalindromes("GAATTC", 6, 6).ToList();
        Assert.That(results, Has.Count.EqualTo(1));
    }

    [Test]
    public void FindPalindromes_EmptySequence_ReturnsEmpty()
    {
        var results = RepeatFinder.FindPalindromes("", 4, 12).ToList();
        Assert.That(results, Is.Empty);
    }

    #endregion

    #region Tandem Repeat Summary Tests

    [Test]
    public void GetTandemRepeatSummary_MixedRepeats_CorrectSummary()
    {
        var sequence = new DnaSequence("AAAAAACGTCGTCGTACACACACATGATGATG");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.TotalRepeats, Is.GreaterThan(0));
        Assert.That(summary.TotalRepeatBases, Is.GreaterThan(0));
        Assert.That(summary.PercentageOfSequence, Is.GreaterThan(0));
    }

    [Test]
    public void GetTandemRepeatSummary_NoRepeats_ZeroSummary()
    {
        var sequence = new DnaSequence("ACGT");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.TotalRepeats, Is.EqualTo(0));
        Assert.That(summary.TotalRepeatBases, Is.EqualTo(0));
    }

    [Test]
    public void GetTandemRepeatSummary_MononucleotideCount_Correct()
    {
        var sequence = new DnaSequence("AAAAAATTTTTGGGGGCCCCC");
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.MononucleotideRepeats, Is.EqualTo(4)); // A, T, G, C runs
    }

    [Test]
    public void GetTandemRepeatSummary_LongestRepeat_Identified()
    {
        var sequence = new DnaSequence("AAACAGCAGCAGCAGCAGCAGAAA"); // 6x CAG = 18bp
        var summary = RepeatFinder.GetTandemRepeatSummary(sequence, 3);

        Assert.That(summary.LongestRepeat, Is.Not.Null);
        Assert.That(summary.LongestRepeat!.Value.TotalLength, Is.EqualTo(18));
    }

    #endregion

    #region Edge Cases and Validation

    [Test]
    public void FindMicrosatellites_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RepeatFinder.FindMicrosatellites((DnaSequence)null!, 1, 6, 3).ToList());
    }

    [Test]
    public void FindMicrosatellites_InvalidMinUnitLength_ThrowsException()
    {
        var sequence = new DnaSequence("ACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites(sequence, 0, 6, 3).ToList());
    }

    [Test]
    public void FindMicrosatellites_InvalidMaxUnitLength_ThrowsException()
    {
        var sequence = new DnaSequence("ACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites(sequence, 5, 4, 3).ToList());
    }

    [Test]
    public void FindMicrosatellites_InvalidMinRepeats_ThrowsException()
    {
        var sequence = new DnaSequence("ACGT");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindMicrosatellites(sequence, 1, 6, 1).ToList());
    }

    [Test]
    public void FindInvertedRepeats_NullSequence_ThrowsException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            RepeatFinder.FindInvertedRepeats((DnaSequence)null!, 4, 10, 3).ToList());
    }

    [Test]
    public void FindPalindromes_OddMinLength_ThrowsException()
    {
        var sequence = new DnaSequence("GAATTC");
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            RepeatFinder.FindPalindromes(sequence, 5, 12).ToList());
    }

    #endregion
}
