using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_clumps</c> MCP tool.
/// Expected results derived from the clump definition in KmerAnalyzer.FindClumps: a k-mer is a
/// clump if it occurs >= minOccurrences times within some window of size windowSize (Compeau &amp;
/// Pevzner Ch.1). Order is unspecified, so results are compared as sets. NOT the wrapper's output.
/// </summary>
[TestFixture]
public class FindClumpsTests
{
    [Test]
    public void FindClumps_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindClumps("AAAA", 2, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumps("", 2, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumps(null!, 2, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumps("AAAA", 0, 4, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumps("AAAA", 3, 2, 3)); // window < k
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindClumps("AAAA", 2, 4, 0)); // minOcc <= 0
    }

    [Test]
    public void FindClumps_Binding_InvokesSuccessfully()
    {
        // "AAAA", k=2, window=4, minOcc=3: window "AAAA" has AA x3 -> clump {AA}.
        var one = AnalysisTools.FindClumps("AAAA", 2, 4, 3).Kmers;
        Assert.That(one, Is.EquivalentTo(new[] { "AA" }));

        // minOcc=4 -> AA only reaches 3 -> no clumps.
        var none = AnalysisTools.FindClumps("AAAA", 2, 4, 4).Kmers;
        Assert.That(none, Is.Empty);

        // "AAAACCCC", k=2, window=4, minOcc=3: "AAAA"->AA, "CCCC"->CC -> {AA, CC}.
        var two = AnalysisTools.FindClumps("AAAACCCC", 2, 4, 3).Kmers;
        Assert.That(two, Is.EquivalentTo(new[] { "AA", "CC" }));
    }
}
