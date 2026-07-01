using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_orthologs</c> MCP tool.
/// Expected values from ComparativeGenomics' own unit test
/// (ComparativeGenomics_FindOrthologs_Tests.MutualBestHits: two identical-sequence pairs,
/// identity/coverage 1.0), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindOrthologsTests
{
    private const string AcgtRepeat14 = "ACGTACGTACGTAC";
    private const string TtBlock = "TTTTGGGGCCCCAAAA";

    private static GeneInput Gene(string id, string genome, string seq)
        => new(id, genome, 0, seq.Length, '+', seq);

    [Test]
    public void FindOrthologs_Schema_ValidatesCorrectly()
    {
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14) };
        Assert.DoesNotThrow(() => AnalysisTools.FindOrthologs(g1, g2));
    }

    [Test]
    public void FindOrthologs_Binding_InvokesSuccessfully()
    {
        var g1 = new[] { Gene("a1", "G1", AcgtRepeat14), Gene("a2", "G1", TtBlock) };
        var g2 = new[] { Gene("b1", "G2", AcgtRepeat14), Gene("b2", "G2", TtBlock) };

        var pairs = AnalysisTools.FindOrthologs(g1, g2).Items;
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Length.EqualTo(2));
            var byG1 = pairs.ToDictionary(p => p.Gene1Id, p => p.Gene2Id);
            Assert.That(byG1["a1"], Is.EqualTo("b1"));
            Assert.That(byG1["a2"], Is.EqualTo("b2"));
            Assert.That(pairs.Select(p => p.Identity), Is.All.EqualTo(1.0).Within(1e-10));
            Assert.That(pairs.Select(p => p.Coverage), Is.All.EqualTo(1.0).Within(1e-10));
        });
    }
}
