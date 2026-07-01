using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_motif_by_prosite</c> MCP tool.
/// Expected coordinates from the PROSITE->regex conversion (N-{P}-[ST] -> N[^P][ST])
/// applied to MNSTV, NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindMotifByPrositeTests
{
    [Test]
    public void FindMotifByProsite_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindMotifByProsite("MNSTV", "N-{P}-[ST]"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByProsite("", "N-{P}-[ST]"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByProsite(null!, "N-{P}-[ST]"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByProsite("MNSTV", ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByProsite("MNSTV", null!));
    }

    [Test]
    public void FindMotifByProsite_Binding_InvokesSuccessfully()
    {
        // N-{P}-[ST] -> N[^P][ST]; NST at index 1 of MNSTV.
        var m = AnalysisTools.FindMotifByProsite("MNSTV", "N-{P}-[ST]", "NGlyc").Items;
        Assert.Multiple(() =>
        {
            Assert.That(m, Has.Length.EqualTo(1));
            Assert.That(m[0].Start, Is.EqualTo(1));
            Assert.That(m[0].End, Is.EqualTo(3));
            Assert.That(m[0].Sequence, Is.EqualTo("NST"));
            Assert.That(m[0].MotifName, Is.EqualTo("NGlyc"));
            Assert.That(m[0].Pattern, Is.EqualTo("N-{P}-[ST]"));
        });

        // No N in MKV -> no match.
        var none = AnalysisTools.FindMotifByProsite("MKV", "N-{P}-[ST]").Items;
        Assert.That(none, Is.Empty);
    }
}
