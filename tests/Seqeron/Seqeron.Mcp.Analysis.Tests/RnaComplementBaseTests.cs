using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>rna_complement_base</c> MCP tool.
/// Expected values from the RNA complement rule (A&lt;-&gt;U, G&lt;-&gt;C), NOT the wrapper output.
/// </summary>
[TestFixture]
public class RnaComplementBaseTests
{
    [Test]
    public void RnaComplementBase_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.RnaComplementBase("A"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaComplementBase(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaComplementBase(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaComplementBase("AU"));
    }

    [Test]
    public void RnaComplementBase_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.RnaComplementBase("A").Complement, Is.EqualTo('U'));
            Assert.That(AnalysisTools.RnaComplementBase("U").Complement, Is.EqualTo('A'));
            Assert.That(AnalysisTools.RnaComplementBase("G").Complement, Is.EqualTo('C'));
            Assert.That(AnalysisTools.RnaComplementBase("C").Complement, Is.EqualTo('G'));
        });
    }

    /// <summary>
    /// A1-12 / F25: documented IUPAC + T behaviour equals Biopython 1.88 <c>complement_rna</c>:
    /// complement_rna("ACGTURYSWKMBDHVN") = "UGCAAYRSWMKVHDBN".
    /// </summary>
    [Test]
    public void RnaComplementBase_IupacAndT_MatchBiopythonComplementRna()
    {
        const string input = "ACGTURYSWKMBDHVN";
        const string expected = "UGCAAYRSWMKVHDBN";
        for (int i = 0; i < input.Length; i++)
            Assert.That(AnalysisTools.RnaComplementBase(input[i].ToString()).Complement, Is.EqualTo(expected[i]), $"base {input[i]}");
    }
}
