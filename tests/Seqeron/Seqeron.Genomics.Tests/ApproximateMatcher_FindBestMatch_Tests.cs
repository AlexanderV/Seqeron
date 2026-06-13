// PAT-APPROX-003 — Best Match and Frequency Analysis
// Evidence: docs/Evidence/PAT-APPROX-003-Evidence.md
// TestSpec: tests/TestSpecs/PAT-APPROX-003.md
// Source: Compeau P, Pevzner P (2015). Bioinformatics Algorithms, ch.1. ROSALIND BA1H/BA1I/BA1N.
using System;
using System.Linq;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Test suite for PAT-APPROX-003: Best Match and Frequency Analysis.
///
/// Canonical methods:
/// - ApproximateMatcher.FindBestMatch (minimum-Hamming-distance window, leftmost tie-break)
/// - ApproximateMatcher.CountApproximateOccurrences (Count_d, BA1H/BA1I)
/// - ApproximateMatcher.FindFrequentKmersWithMismatches (Frequent Words with Mismatches, BA1I)
///
/// Evidence sources (all retrieved 2026-06-13):
/// - ROSALIND BA1H: https://rosalind.info/problems/ba1h/
/// - ROSALIND BA1I: https://rosalind.info/problems/ba1i/
/// - ROSALIND BA1N: https://rosalind.info/problems/ba1n/
/// </summary>
[TestFixture]
[Category("PAT-APPROX-003")]
[Category("Pattern Matching")]
[Category("Approximate Matching")]
public class ApproximateMatcher_FindBestMatch_Tests
{
    #region FindFrequentKmersWithMismatches (BA1I)

    // M1 — ROSALIND BA1I sample: ACGTTGCATGTCGCATGATGCATGAGAGCT, k=4, d=1
    // Expected most-frequent 4-mers (set): GATG, ATGC, ATGT; max Count_d = 5.
    [Test]
    [Description("BA1I sample: returns the published set of most-frequent 4-mers with d=1, each count 5")]
    public void FrequentKmers_Ba1iSample_ReturnsThreeMostFrequentKmers()
    {
        var result = ApproximateMatcher
            .FindFrequentKmersWithMismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
            .ToList();

        var kmers = result.Select(r => r.Kmer).OrderBy(s => s).ToArray();
        Assert.Multiple(() =>
        {
            Assert.That(kmers, Is.EqualTo(new[] { "ATGC", "ATGT", "GATG" }),
                "BA1I sample most-frequent 4-mers (d=1) must equal the published set {GATG, ATGC, ATGT}");
            Assert.That(result.Select(r => r.Count), Is.All.EqualTo(5),
                "Each most-frequent 4-mer must have Count_1 = 5 per the BA1I sample");
        });
    }

    // M2 — INV-03: all tied maxima are returned (the BA1I sample has exactly three).
    [Test]
    [Description("BA1I: every k-mer achieving the maximum Count_d is returned (3 ties)")]
    public void FrequentKmers_Ba1iSample_ReturnsAllTiedKmers()
    {
        var result = ApproximateMatcher
            .FindFrequentKmersWithMismatches("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
            .ToList();

        Assert.That(result, Has.Count.EqualTo(3),
            "All k-mers tying for the maximum Count_d must be returned; BA1I sample has exactly 3");
    }

    // S1 — d=0 reduces FrequentWords to exact frequent k-mer (Neighbors(P,0)={P}).
    [Test]
    [Description("d=0 degenerates to the exact most-frequent k-mer")]
    public void FrequentKmers_ZeroMismatch_ReducesToExactFrequentKmer()
    {
        // AAAAAA has three overlapping exact 4-mers AAAA at positions 0,1,2.
        var result = ApproximateMatcher.FindFrequentKmersWithMismatches("AAAAAA", 4, 0).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(result, Has.Count.EqualTo(1),
                "With d=0 only the single exact 4-mer is most frequent");
            Assert.That(result[0].Kmer, Is.EqualTo("AAAA"),
                "The exact most-frequent 4-mer of AAAAAA is AAAA");
            Assert.That(result[0].Count, Is.EqualTo(3),
                "AAAA occurs exactly 3 times (positions 0,1,2) so Count_0 = 3");
        });
    }

    // C2 — invalid k or d.
    [Test]
    [Description("k <= 0 throws ArgumentOutOfRangeException")]
    public void FrequentKmers_InvalidKOrD_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ApproximateMatcher.FindFrequentKmersWithMismatches("ACGT", 0, 1).ToList(),
                "k = 0 is invalid (k must be positive)");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => ApproximateMatcher.FindFrequentKmersWithMismatches("ACGT", 2, -1).ToList(),
                "d < 0 is invalid (mismatch budget cannot be negative)");
        });
    }

    // C3 — empty sequence yields empty result.
    [Test]
    [Description("Empty sequence yields no frequent k-mers")]
    public void FrequentKmers_EmptySequence_ReturnsEmpty()
    {
        var result = ApproximateMatcher.FindFrequentKmersWithMismatches("", 4, 1).ToList();
        Assert.That(result, Is.Empty, "An empty sequence contains no k-mers");
    }

    #endregion

    #region CountApproximateOccurrences (Count_d, BA1H/BA1I)

    // M3 — ROSALIND BA1H sample Count_d.
    [Test]
    [Description("BA1H sample: Count_3(Text, ATTCTGGA) = 5")]
    public void CountApproximateOccurrences_Ba1hSample_ReturnsFive()
    {
        const string text =
            "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC";
        int count = ApproximateMatcher.CountApproximateOccurrences(text, "ATTCTGGA", 3);
        Assert.That(count, Is.EqualTo(5),
            "BA1H sample: ATTCTGGA occurs at 5 positions with at most 3 mismatches");
    }

    // M4 — ROSALIND BA1H sample positions (0-based).
    [Test]
    [Description("BA1H sample: approximate-occurrence start positions are 6 7 26 27 78")]
    public void FindWithMismatches_Ba1hSample_ReturnsExpectedPositions()
    {
        const string text =
            "CGCCCGAATCCAGAACGCATTCCCATATTTCGGGACCACTGGCCTCCACGGTACGGACGTCAATCAAATGCCTAGCGGCTTGTGGTTTCTCCTACGCTCC";
        var positions = ApproximateMatcher.FindWithMismatches(text, "ATTCTGGA", 3)
            .Select(r => r.Position).OrderBy(p => p).ToArray();
        Assert.That(positions, Is.EqualTo(new[] { 6, 7, 26, 27, 78 }),
            "BA1H published sample output: starting positions 6 7 26 27 78");
    }

    // M5 — ROSALIND BA1I Count_1 worked example.
    [Test]
    [Description("BA1I worked example: Count_1(AACAAGCTGATAAACATTTAAAGAG, AAAAA) = 4")]
    public void CountApproximateOccurrences_Count1WorkedExample_ReturnsFour()
    {
        int count = ApproximateMatcher.CountApproximateOccurrences(
            "AACAAGCTGATAAACATTTAAAGAG", "AAAAA", 1);
        Assert.That(count, Is.EqualTo(4),
            "Worked example: AAAAA matches at AACAA, ATAAA, AAACA, AAAGA (4 windows within 1 mismatch)");
    }

    // M6 — Count_0 equals exact occurrence count (INV-01).
    [Test]
    [Description("Count_0 equals exact occurrence count")]
    public void CountApproximateOccurrences_ExactZeroMismatch_EqualsExactCount()
    {
        int count = ApproximateMatcher.CountApproximateOccurrences("ACGTACGT", "ACGT", 0);
        Assert.That(count, Is.EqualTo(2),
            "ACGT occurs exactly twice in ACGTACGT; with d=0 Count_d = exact count");
    }

    // S3 — case-insensitivity (input is upper-cased internally).
    [Test]
    [Description("Lowercase input is matched case-insensitively")]
    public void CountApproximateOccurrences_LowercaseInput_MatchesUppercase()
    {
        int count = ApproximateMatcher.CountApproximateOccurrences("acgtacgt", "ACGT", 0);
        Assert.That(count, Is.EqualTo(2),
            "Lowercase sequence must produce the same Count_d as its uppercase form");
    }

    // C4 — empty / oversized pattern.
    [Test]
    [Description("Empty inputs or pattern longer than sequence give Count 0")]
    public void CountApproximateOccurrences_EmptyOrTooLongPattern_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.CountApproximateOccurrences("", "ACGT", 1), Is.EqualTo(0),
                "Empty sequence has no windows");
            Assert.That(ApproximateMatcher.CountApproximateOccurrences("ACGT", "", 1), Is.EqualTo(0),
                "Empty pattern yields no occurrences");
            Assert.That(ApproximateMatcher.CountApproximateOccurrences("AC", "ACGT", 1), Is.EqualTo(0),
                "Pattern longer than sequence has no equal-length window");
        });
    }

    #endregion

    #region FindBestMatch (minimum-Hamming-distance window)

    // M7 — exact match: distance 0 at position 0 (INV-04).
    [Test]
    [Description("Exact match returns distance 0 at the leftmost position and IsExact")]
    public void FindBestMatch_ExactMatch_ReturnsZeroDistanceAtPositionZero()
    {
        var result = ApproximateMatcher.FindBestMatch("ACGTACGT", "ACGT");

        Assert.That(result, Is.Not.Null, "An exact match exists, so a result must be returned");
        Assert.Multiple(() =>
        {
            Assert.That(result!.Value.Distance, Is.EqualTo(0),
                "An exact occurrence has Hamming distance 0");
            Assert.That(result!.Value.IsExact, Is.True,
                "IsExact must be true when the minimum distance is 0");
            Assert.That(result!.Value.Position, Is.EqualTo(0),
                "Leftmost exact occurrence of ACGT in ACGTACGT is at position 0");
        });
    }

    // M8 — no exact match: leftmost minimum-distance window (INV-04, INV-05).
    [Test]
    [Description("No exact match: returns the leftmost minimum-Hamming-distance window")]
    public void FindBestMatch_NoExactMatch_ReturnsLeftmostMinimumDistanceWindow()
    {
        // TTTTTTTT vs ACGT: every length-4 window is TTTT; HammingDistance(ACGT,TTTT)=3
        // (only the last position T matches). Leftmost minimal window is position 0.
        var result = ApproximateMatcher.FindBestMatch("TTTTTTTT", "ACGT");

        Assert.That(result, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(result!.Value.Distance, Is.EqualTo(3),
                "HammingDistance(ACGT, TTTT) = 3 (only the final T matches)");
            Assert.That(result!.Value.Position, Is.EqualTo(0),
                "Leftmost window achieving the minimum distance is at position 0");
            Assert.That(result!.Value.MatchedSequence, Is.EqualTo("TTTT"),
                "The reported best window is the 4-mer TTTT");
            Assert.That(result!.Value.IsExact, Is.False,
                "A non-zero minimum distance is not an exact match");
        });
    }

    // S2 — tie-break convention: leftmost window among equal minima (INV-05).
    [Test]
    [Description("Tied minimum-distance windows: the leftmost is returned")]
    public void FindBestMatch_TiedMinima_ReturnsLeftmostWindow()
    {
        // ACGTACGA vs TTTT: window ACGT (pos 0) has Hamming distance 3 (only T at pos 3 matches),
        // which is the minimum; the leftmost such window is position 0.
        var result = ApproximateMatcher.FindBestMatch("ACGTACGA", "TTTT");

        Assert.That(result, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(result!.Value.Distance, Is.EqualTo(3),
                "Minimum Hamming distance over all 4-mer windows of ACGTACGA vs TTTT is 3");
            Assert.That(result!.Value.Position, Is.EqualTo(0),
                "Leftmost window achieving the minimum distance is at position 0 (ACGT)");
            Assert.That(result!.Value.MatchedSequence, Is.EqualTo("ACGT"),
                "Leftmost minimal window is the 4-mer ACGT");
        });
    }

    // C1 — empty / too-short inputs return null.
    [Test]
    [Description("Empty inputs or pattern longer than sequence return null")]
    public void FindBestMatch_EmptyOrTooShort_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ApproximateMatcher.FindBestMatch("", "ACGT"), Is.Null,
                "Empty sequence has no window to match");
            Assert.That(ApproximateMatcher.FindBestMatch("ACGT", ""), Is.Null,
                "Empty pattern yields no match");
            Assert.That(ApproximateMatcher.FindBestMatch("AC", "ACGT"), Is.Null,
                "Pattern longer than sequence has no equal-length window");
        });
    }

    #endregion
}
