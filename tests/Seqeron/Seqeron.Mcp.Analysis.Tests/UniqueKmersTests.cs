using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>unique_kmers</c> MCP tool.
/// Expected values from the singleton (count==1) definition, NOT the wrapper output.
/// </summary>
[TestFixture]
public class UniqueKmersTests
{
    [Test]
    public void UniqueKmers_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.UniqueKmers("ATGATG", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.UniqueKmers("", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.UniqueKmers(null!, 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.UniqueKmers("ATGATG", 0));
        Assert.Throws<ArgumentException>(() => AnalysisTools.UniqueKmers("ATGATG", -1));
    }

    [Test]
    public void UniqueKmers_Binding_InvokesSuccessfully()
    {
        // "ATGATG" k=3: ATG:2 (not unique), TGA:1, GAT:1 -> {TGA, GAT}.
        var uniq = AnalysisTools.UniqueKmers("ATGATG", 3).Kmers;
        Assert.That(uniq, Is.EqualTo(new[] { "GAT", "TGA" }), "lexicographic order");

        // "AAAA" k=2: AA:3 -> no singletons.
        var none = AnalysisTools.UniqueKmers("AAAA", 2).Kmers;
        Assert.That(none, Is.Empty);
    }

    [Test]
    public void UniqueKmers_JellyfishModes_DelegateToOptionAwareOverload()
    {
        // Jellyfish 2.3.1: count -C -m 3 + dump -c -L 1 -U 1 on ATGATG -> ATC:1 TCA:1;
        // plain (ACGT-only) count on acgtNNacgtacgRtTTGCAnA k=3 -> GCA GTA TAC TGC TTG TTT.
        Assert.That(AnalysisTools.UniqueKmers("ATGATG", 3, canonical: true).Kmers, Is.EqualTo(new[] { "ATC", "TCA" }));
        Assert.That(AnalysisTools.UniqueKmers("acgtNNacgtacgRtTTGCAnA", 3, acgtOnly: true).Kmers,
            Is.EqualTo(new[] { "GCA", "GTA", "TAC", "TGC", "TTG", "TTT" }));
        Assert.That(AnalysisTools.UniqueKmers("acgtNNacgtacgRtTTGCAnA", 3, canonical: true).Kmers,
            Is.EqualTo(new[] { "AAA", "CAA" }));
    }
}
