using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>kmer_distance</c> MCP tool.
/// Expected values derived from the Euclidean word-composition distance definition
/// (identity => 0; orthogonal monomer vectors => sqrt(2)), NOT the wrapper output.
/// </summary>
[TestFixture]
public class KmerDistanceTests
{
    [Test]
    public void KmerDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.KmerDistance("AAAA", "TTTT", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("", "TTTT", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("AAAA", "", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance(null!, "TTTT", 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("AAAA", null!, 1));
        Assert.Throws<ArgumentException>(() => AnalysisTools.KmerDistance("AAAA", "TTTT", 0));
    }

    [Test]
    public void KmerDistance_Binding_InvokesSuccessfully()
    {
        // Identical sequences => distance 0.
        var same = AnalysisTools.KmerDistance("ACGTACGT", "ACGTACGT", 2).Distance;
        Assert.That(same, Is.EqualTo(0.0).Within(1e-12));

        // {A:1} vs {T:1} => sqrt(1^2 + 1^2) = sqrt(2).
        var orth = AnalysisTools.KmerDistance("AAAA", "TTTT", 1).Distance;
        Assert.That(orth, Is.EqualTo(Math.Sqrt(2)).Within(1e-12));
    }
}
