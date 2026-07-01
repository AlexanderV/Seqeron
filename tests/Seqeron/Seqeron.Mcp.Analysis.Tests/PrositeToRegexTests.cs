using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>prosite_to_regex</c> MCP tool.
/// Expected regex strings derived by hand from the documented PROSITE->regex rules,
/// NOT the wrapper output.
/// </summary>
[TestFixture]
public class PrositeToRegexTests
{
    [Test]
    public void PrositeToRegex_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.PrositeToRegex("N-{P}-[ST]-{P}"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PrositeToRegex(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.PrositeToRegex(null!));
    }

    [Test]
    public void PrositeToRegex_Binding_InvokesSuccessfully()
    {
        // N-glycosylation PS00001: exclusion {P} -> [^P], separators dropped.
        Assert.That(AnalysisTools.PrositeToRegex("N-{P}-[ST]-{P}").Regex, Is.EqualTo("N[^P][ST][^P]"));

        // Variable gap x(2) -> .{2}; class preserved.
        Assert.That(AnalysisTools.PrositeToRegex("[AC]-x(2)-V").Regex, Is.EqualTo("[AC].{2}V"));

        // Single x -> '.'
        Assert.That(AnalysisTools.PrositeToRegex("A-x-C").Regex, Is.EqualTo("A.C"));
    }
}
