using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

/// <summary>
/// R14-17 (B09 F298): <c>find_longest_repeat</c> / <c>find_longest_common_region</c> / <c>calculate_similarity</c> take N / IUPAC codes / lower case (string overloads of GenomicAnalyzer, no strict
/// ACGT DnaSequence), the both-strand option, the strands and maxMismatches options and a cancellation token. Expected values: independent brute force (Python) over every window.
/// </summary>
[TestFixture]
[Category("McpCore")]
public class FindLongestRepeatCommonRegionOptionsTests
{
    [Test]
    public void LongestCommonRegion_WithNRuns_EndsTheRegionAtTheN()
    {
        var r = SuffixTreeGenomicsTools.FindLongestCommonRegion("acgtNacgtac", "TTNACGTAA");
        Assert.That((r.Region, r.Position1, r.Position2, r.Length), Is.EqualTo(("ACGTA", 5, 3, 5)));
        Assert.That(r.IsReverseComplement, Is.False);
        Assert.That(r.Mismatches, Is.Null);
    }

    [Test]
    public void LongestCommonRegion_ReverseComplementAndBothStrands()
    {
        // GGGATCCCAT is the reverse complement of ATGGGATCCC: the whole 10-mer on the reverse strand; forward only 8
        var fwd = SuffixTreeGenomicsTools.FindLongestCommonRegion("GGGATCCCAT", "ATGGGATCCC");
        Assert.That((fwd.Region, fwd.Position1, fwd.Position2), Is.EqualTo(("GGGATCCC", 0, 2)));
        var rev = SuffixTreeGenomicsTools.FindLongestCommonRegion("GGGATCCCAT", "ATGGGATCCC", "ReverseComplement");
        Assert.That((rev.Region, rev.Position1, rev.Position2, rev.Length, rev.IsReverseComplement), Is.EqualTo(("GGGATCCCAT", 0, 0, 10, true)));
        var both = SuffixTreeGenomicsTools.FindLongestCommonRegion("GGGATCCCAT", "ATGGGATCCC", "both");
        Assert.That((both.Length, both.IsReverseComplement), Is.EqualTo((10, true)));
    }

    [Test]
    public void LongestCommonRegion_WithMismatches_IsTheKMismatchLcs()
    {
        // brute force: 1 mismatch gives the 8-mer pair ACGTACGT / ACGAACGT, 0 mismatches only ACGT at (0, 4)
        var k1 = SuffixTreeGenomicsTools.FindLongestCommonRegion("ACGTACGTTTGA", "ACGAACGTCTGA", maxMismatches: 1);
        Assert.That((k1.Region, k1.Region2, k1.Position1, k1.Position2, k1.Length, k1.Mismatches), Is.EqualTo(("ACGTACGT", "ACGAACGT", 0, 0, 8, 1)));
        var k0 = SuffixTreeGenomicsTools.FindLongestCommonRegion("ACGTACGTTTGA", "ACGAACGTCTGA");
        Assert.That((k0.Region, k0.Position1, k0.Position2), Is.EqualTo(("ACGT", 0, 4)));
    }

    [Test]
    public void LongestRepeat_BothStrands_ReportsTheReversePositions()
    {
        // brute force over S # revcomp(S): the longest class is TGGATCCA (a reverse-complement palindrome), once on the forward strand at 9
        var r = SuffixTreeGenomicsTools.FindLongestRepeat("ACGGATCCTTGGATCCAGT", bothStrands: true);
        Assert.That(r.Repeat, Is.EqualTo("TGGATCCA"));
        Assert.That(r.Length, Is.EqualTo(8));
        Assert.That(r.Positions, Is.EqualTo(new[] { 9 }));
        Assert.That(r.ReversePositions, Is.EqualTo(new[] { 9 }));
        var plain = SuffixTreeGenomicsTools.FindLongestRepeat("ACGGATCCTTGGATCCAGT");
        Assert.That(plain.ReversePositions, Is.Null);
    }

    [Test]
    public void LongestRepeat_NRunsNeverRepeat_AndLowerCaseIsFolded()
    {
        var r = SuffixTreeGenomicsTools.FindLongestRepeat("acgtNacgt");
        Assert.That((r.Repeat, r.Positions), Is.EqualTo(("ACGT", new[] { 0, 5 })));
        var broken = SuffixTreeGenomicsTools.FindLongestRepeat("ACNGTACNGT");
        Assert.That(broken.Repeat.Contains('N'), Is.False);
    }

    [Test]
    public void CancelledToken_StopsTheTools()
    {
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => SuffixTreeGenomicsTools.FindLongestRepeat("ACGTACGT", cancellationToken: cts.Token));
        Assert.Throws<OperationCanceledException>(() => SuffixTreeGenomicsTools.FindLongestCommonRegion("ACGT", "ACGT", cancellationToken: cts.Token));
        Assert.Throws<OperationCanceledException>(() => SuffixTreeGenomicsTools.CalculateSimilarity("ACGTACGT", "ACGTACGA", cancellationToken: cts.Token));
    }

    [Test]
    public void CalculateSimilarity_IupacSequences_SkipWindowsWithAnN()
    {
        // 4-mers of ACGTNACGT: ACGT, CGTN (skipped), GTNA (skipped), TNAC (skipped), NACG (skipped), ACGT; so the distinct set equals that of ACGT
        var r = SuffixTreeGenomicsTools.CalculateSimilarity("acgtNacgt", "ACGT", 4);
        Assert.That(r.Similarity, Is.EqualTo(100.0));
    }
}
