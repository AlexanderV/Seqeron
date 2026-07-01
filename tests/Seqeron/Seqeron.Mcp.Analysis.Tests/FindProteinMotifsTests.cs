using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_protein_motifs</c> MCP tool.
/// Expected values from ProteinMotifFinder's own unit test
/// (ProteinMotifFinder_FindCommonMotifs_Tests, PS00001 N-{P}-[ST]-{P} on "AAAANFTAAAA"),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindProteinMotifsTests
{
    [Test]
    public void FindProteinMotifs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindProteinMotifs("AAAANFTAAAA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindProteinMotifs(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindProteinMotifs(null!));
    }

    [Test]
    public void FindProteinMotifs_Binding_InvokesSuccessfully()
    {
        // PS00001 N-glycosylation matches the NFTA window at index 4 exactly once.
        var hits = AnalysisTools.FindProteinMotifs("AAAANFTAAAA").Items
            .Where(m => m.MotifName == "ASN_GLYCOSYLATION").ToArray();
        Assert.Multiple(() =>
        {
            Assert.That(hits, Has.Length.EqualTo(1));
            Assert.That(hits[0].Start, Is.EqualTo(4));
            Assert.That(hits[0].End, Is.EqualTo(7));
            Assert.That(hits[0].Sequence, Is.EqualTo("NFTA"));
        });

        // Poly-alanine has no N-glycosylation site.
        var none = AnalysisTools.FindProteinMotifs("AAAAAA").Items
            .Where(m => m.MotifName == "ASN_GLYCOSYLATION");
        Assert.That(none, Is.Empty);
    }
}
