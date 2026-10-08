using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// frequent_kmers_with_mismatches_and_revcomp — delegation to
/// <c>ApproximateMatcher.FindFrequentKmersWithMismatchesAndReverseComplements</c> (ROSALIND BA1J).
/// Expected values: BA1J sample + brute-force cases locked in ApproximateMatcher_FindBestMatch_Tests.
/// </summary>
[TestFixture]
public class FrequentKmersWithMismatchesAndRevcompTests
{
    [Test]
    public void FrequentKmersWithMismatchesAndRevcomp_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AlignmentTools.FrequentKmersWithMismatchesAndRevcomp("ACGT", 2, 0));
        Assert.Throws<ArgumentException>(() => AlignmentTools.FrequentKmersWithMismatchesAndRevcomp("", 2, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FrequentKmersWithMismatchesAndRevcomp("ACGT", 0, 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.FrequentKmersWithMismatchesAndRevcomp("ACGT", 2, -1));
    }

    [TestCase("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1, new[] { "ACAT", "ATGT" }, 9)]
    [TestCase("ACGT", 4, 0, new[] { "ACGT" }, 2)]
    [TestCase("AAAAA", 2, 0, new[] { "AA", "TT" }, 4)]
    [TestCase("ACGTTNCATGTCGCATGATGCANGAGAGCT", 4, 1, new[] { "CATG" }, 8)]
    public void FrequentKmersWithMismatchesAndRevcomp_MatchesBa1jReference(string text, int k, int d, string[] expected, int score)
    {
        var r = AlignmentTools.FrequentKmersWithMismatchesAndRevcomp(text, k, d);
        Assert.That(r.Items.Select(i => i.Kmer), Is.EqualTo(expected), "sorted ordinal");
        Assert.That(r.Items.Select(i => i.Count), Is.All.EqualTo(score));

        var lib = global::Seqeron.Genomics.Alignment.ApproximateMatcher
            .FindFrequentKmersWithMismatchesAndReverseComplements(text, k, d)
            .OrderBy(t => t.Kmer, StringComparer.Ordinal).ToArray();
        Assert.That(r.Items.Select(i => (i.Kmer, i.Count)), Is.EqualTo(lib));
    }
}
