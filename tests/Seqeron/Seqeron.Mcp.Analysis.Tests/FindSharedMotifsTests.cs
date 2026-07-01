using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_shared_motifs</c> MCP tool.
/// Expected values from the "k-mer in >= minSequences" definition (prevalence =
/// fraction of sequences), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindSharedMotifsTests
{
    [Test]
    public void FindSharedMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindSharedMotifs(new[] { "ATGCATGC", "ATGCTTTT" }, 4, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindSharedMotifs(null!, 4, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindSharedMotifs(Array.Empty<string>(), 4, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindSharedMotifs(new[] { "XYZ" }, 4, 2));
    }

    [Test]
    public void FindSharedMotifs_Binding_InvokesSuccessfully()
    {
        // ATGC is in both sequences -> prevalence 1.0, indices {0,1}.
        var shared = AnalysisTools.FindSharedMotifs(new[] { "ATGCATGC", "ATGCTTTT" }, 4, 2).Items;
        var atgc = shared.Single(m => m.Sequence == "ATGC");
        Assert.Multiple(() =>
        {
            Assert.That(atgc.SequenceIndices, Is.EquivalentTo(new[] { 0, 1 }));
            Assert.That(atgc.Prevalence, Is.EqualTo(1.0).Within(1e-9));
        });

        // Disjoint sequences share no 4-mer.
        var none = AnalysisTools.FindSharedMotifs(new[] { "AAAA", "TTTT" }, 4, 2).Items;
        Assert.That(none, Is.Empty);
    }
}
