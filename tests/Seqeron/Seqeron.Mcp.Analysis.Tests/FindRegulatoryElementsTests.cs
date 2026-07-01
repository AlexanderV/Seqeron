using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_regulatory_elements</c> MCP tool.
/// Expected values from the built-in TATA box consensus (TATAAA) catalog entry,
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindRegulatoryElementsTests
{
    [Test]
    public void FindRegulatoryElements_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindRegulatoryElements("GGTATAAAGG"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRegulatoryElements(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRegulatoryElements(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindRegulatoryElements("XYZ123"));
    }

    [Test]
    public void FindRegulatoryElements_Binding_InvokesSuccessfully()
    {
        // TATA box consensus TATAAA at position 2.
        var tata = AnalysisTools.FindRegulatoryElements("GGTATAAAGG").Items;
        var box = tata.Single(e => e.Name == "TATA Box");
        Assert.Multiple(() =>
        {
            Assert.That(box.Position, Is.EqualTo(2));
            Assert.That(box.Sequence, Is.EqualTo("TATAAA"));
            Assert.That(box.Pattern, Is.EqualTo("TATAAA"));
        });

        // Poly-A tract matches none of the catalog motifs.
        var none = AnalysisTools.FindRegulatoryElements("AAAAAAAA").Items;
        Assert.That(none, Is.Empty);
    }
}
