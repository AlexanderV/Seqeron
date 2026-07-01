using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_motif_by_pattern</c> MCP tool.
/// Expected coordinates from the regex-match definition (0-based start/end inclusive),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindMotifByPatternTests
{
    [Test]
    public void FindMotifByPattern_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindMotifByPattern("MKV", "K"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByPattern("", "K"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByPattern(null!, "K"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByPattern("MKV", ""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindMotifByPattern("MKV", null!));
    }

    [Test]
    public void FindMotifByPattern_Binding_InvokesSuccessfully()
    {
        // "K" matches at index 1 of "MKV".
        var k = AnalysisTools.FindMotifByPattern("MKV", "K", "Lys").Items;
        Assert.Multiple(() =>
        {
            Assert.That(k, Has.Length.EqualTo(1));
            Assert.That(k[0].Start, Is.EqualTo(1));
            Assert.That(k[0].End, Is.EqualTo(1));
            Assert.That(k[0].Sequence, Is.EqualTo("K"));
            Assert.That(k[0].MotifName, Is.EqualTo("Lys"));
        });

        // No W in the sequence.
        var none = AnalysisTools.FindMotifByPattern("MKV", "W").Items;
        Assert.That(none, Is.Empty);
    }
}
