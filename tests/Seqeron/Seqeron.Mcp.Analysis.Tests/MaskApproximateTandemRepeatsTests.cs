using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>mask_approximate_tandem_repeats</c> MCP tool. Expected = the masked file of compiled TRF 4.10.0
/// (`trf craft.fa 2 7 7 80 10 50 500 -m`), locked in RepeatFinder_TrfParameters_Tests. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class MaskApproximateTandemRepeatsTests
{
    private const string U3 =
        "CACACACACACACACACGCACACACACACCGACACACACACACACACACACACACACACAAGATCGTCAACGACCATCGTTCTTCTGCTTTAAGTGTGAGTTCTCTCTTACAAGCTAATACACCCCAGCGTTCTCCGTAC";

    [Test]
    public void MaskApproximateTandemRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.MaskApproximateTandemRepeats(U3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MaskApproximateTandemRepeats(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MaskApproximateTandemRepeats(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MaskApproximateTandemRepeats(U3, matchProbability: 85));
    }

    [Test]
    public void MaskApproximateTandemRepeats_Binding_EqualsTrfMaskedFile()
    {
        var masked = AnalysisTools.MaskApproximateTandemRepeats(U3).Masked;
        Assert.That(masked, Is.EqualTo(new string('N', 60) + U3[60..]));
    }

    [Test]
    public void MaskApproximateTandemRepeats_SoftMask_LowercasesTheRepeat()
    {
        var masked = AnalysisTools.MaskApproximateTandemRepeats(U3, softMask: true).Masked;
        Assert.That(masked, Is.EqualTo(U3[..60].ToLowerInvariant() + U3[60..]));
    }

    // RepeatFinder_TrfDetection_Tests D2: a table of zeros admits the period-28 candidate (148..209) TRF rejects.
    private const string Spread28 =
        "CAGTCACGGGCTCTGGATCCAGCAGCAGTGCAGCATGTTGGTACCCTATCCCCATACGACACTGTTTGGCGCTGTTGGTTTATGCACGAGTCGTTACTAT"
        + "ATAAAGACCTCGAAGTGCCAGAATTCATCTTTGACCTCAGCGCGTTCGTACTCCGATCGGAACCGCCCGTTCACTGTACTCCGATCGGAACCGCCCCGAT"
        + "ATGTACTCCATTAATCGTCCCTTTGAATTCGGAGATACGCGTGACGGACGTATCGCGTCTCCATTCTTAGCCGACTCCACGACCTCCTTAATGGTTAATC"
        + "AACATAAGAATATTCCCAGGAG";

    [Test]
    public void MaskApproximateTandemRepeats_ApparentSizeTable_IsApplied()
    {
        string zeros = string.Join(",", Enumerable.Repeat(0, 2001));
        Assert.That(AnalysisTools.MaskApproximateTandemRepeats(Spread28).Masked, Is.EqualTo(Spread28));
        Assert.That(AnalysisTools.MaskApproximateTandemRepeats(Spread28, apparentSizeTable: zeros).Masked,
            Is.EqualTo(Spread28[..147] + new string('N', 62) + Spread28[209..]));
    }
}
