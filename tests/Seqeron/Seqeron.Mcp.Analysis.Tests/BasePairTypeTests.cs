using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>base_pair_type</c> MCP tool.
/// Expected classifications follow RnaSecondaryStructure.GetBasePairType:
/// A-U / G-C -> WatsonCrick, G-U -> Wobble, otherwise null.
/// </summary>
[TestFixture]
public class BasePairTypeTests
{
    [Test]
    public void BasePairType_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.BasePairType("A", "U"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.BasePairType("", "U"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.BasePairType("A", null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.BasePairType("AU", "U"));
    }

    [Test]
    public void BasePairType_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.BasePairType("A", "U").Type, Is.EqualTo("WatsonCrick"));
            Assert.That(AnalysisTools.BasePairType("U", "A").Type, Is.EqualTo("WatsonCrick"));
            Assert.That(AnalysisTools.BasePairType("G", "C").Type, Is.EqualTo("WatsonCrick"));
            Assert.That(AnalysisTools.BasePairType("G", "U").Type, Is.EqualTo("Wobble"));
            Assert.That(AnalysisTools.BasePairType("U", "G").Type, Is.EqualTo("Wobble"));
            // Non-pairing -> null.
            Assert.That(AnalysisTools.BasePairType("A", "G").Type, Is.Null);
            Assert.That(AnalysisTools.BasePairType("A", "C").Type, Is.Null);
            // Case-insensitive.
            Assert.That(AnalysisTools.BasePairType("a", "u").Type, Is.EqualTo("WatsonCrick"));
        });
    }
}
