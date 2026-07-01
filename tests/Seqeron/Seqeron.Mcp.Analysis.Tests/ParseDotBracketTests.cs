using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>parse_dot_bracket</c> MCP tool.
/// Expected pairs from stack-based dot-bracket parsing ((0,3),(1,2) for "(())"),
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class ParseDotBracketTests
{
    [Test]
    public void ParseDotBracket_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.ParseDotBracket("(())"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ParseDotBracket(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ParseDotBracket(null!));
    }

    [Test]
    public void ParseDotBracket_Binding_InvokesSuccessfully()
    {
        // "(())" -> pairs (0,3) and (1,2).
        var pairs = AnalysisTools.ParseDotBracket("(())").Pairs;
        Assert.Multiple(() =>
        {
            Assert.That(pairs, Has.Length.EqualTo(2));
            Assert.That(pairs.Select(p => (p.Position1, p.Position2)),
                Is.EquivalentTo(new[] { (0, 3), (1, 2) }));
        });

        // Fully unpaired -> no pairs.
        var none = AnalysisTools.ParseDotBracket("....").Pairs;
        Assert.That(none, Is.Empty);
    }
}
