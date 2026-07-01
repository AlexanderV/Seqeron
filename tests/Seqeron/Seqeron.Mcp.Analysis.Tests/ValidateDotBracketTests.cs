using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>validate_dot_bracket</c> MCP tool.
/// Expected values from the balanced-bracket definition (family-aware), NOT the wrapper output.
/// </summary>
[TestFixture]
public class ValidateDotBracketTests
{
    [Test]
    public void ValidateDotBracket_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.ValidateDotBracket("(())"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ValidateDotBracket(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.ValidateDotBracket(null!));
    }

    [Test]
    public void ValidateDotBracket_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Balanced nested and dotted structures.
            Assert.That(AnalysisTools.ValidateDotBracket("(())").Result, Is.True);
            Assert.That(AnalysisTools.ValidateDotBracket("(((...)))").Result, Is.True);
            // Unbalanced: missing closing bracket.
            Assert.That(AnalysisTools.ValidateDotBracket("(()").Result, Is.False);
            // Family mismatch: '(' closed by ']'.
            Assert.That(AnalysisTools.ValidateDotBracket("(]").Result, Is.False);
        });
    }
}
