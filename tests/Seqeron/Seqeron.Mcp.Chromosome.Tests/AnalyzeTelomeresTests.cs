using NUnit.Framework;
using Seqeron.Mcp.Chromosome.Tools;

namespace Seqeron.Mcp.Chromosome.Tests;

/// <summary>
/// Tests for the <c>analyze_telomeres</c> MCP tool.
/// Expected values derive from ChromosomeAnalyzer.AnalyzeTelomeres: the 3' end is matched against
/// the repeat unit and the 5' end against its reverse complement, walking inward while per-window
/// similarity stays >= 0.7.
/// </summary>
[TestFixture]
public class AnalyzeTelomeresTests
{
    private static string Repeat(string unit, int n) => string.Concat(Enumerable.Repeat(unit, n));

    [Test]
    public void AnalyzeTelomeres_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ChromosomeTools.AnalyzeTelomeres("chr1", "ACGTACGT"));

        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeTelomeres("chr1", ""));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeTelomeres("chr1", null!));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeTelomeres("", "ACGT"));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeTelomeres("chr1", "ACGT", telomereRepeat: ""));
        Assert.Throws<ArgumentException>(() => ChromosomeTools.AnalyzeTelomeres("chr1", "ACGT", searchLength: 0));
    }

    [Test]
    public void AnalyzeTelomeres_Binding_InvokesSuccessfully()
    {
        // 100 tandem TTAGGG = 600 bp: pure 3' telomere; 5' end starts with TTAGGG, not CCCTAA.
        var result = ChromosomeTools.AnalyzeTelomeres("chrT", Repeat("TTAGGG", 100));

        Assert.Multiple(() =>
        {
            Assert.That(result.Chromosome, Is.EqualTo("chrT"));
            Assert.That(result.Has3PrimeTelomere, Is.True);
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(600));
            Assert.That(result.RepeatPurity3Prime, Is.EqualTo(1.0).Within(1e-9));
            Assert.That(result.Has5PrimeTelomere, Is.False);
            Assert.That(result.TelomereLength5Prime, Is.EqualTo(0));
            // 600 < criticalLength (3000) and 3' telomere present -> critically short.
            Assert.That(result.IsCriticallyShort, Is.True);
        });
    }

    [Test]
    public void AnalyzeTelomeres_NonTelomericSequence_NoTelomere()
    {
        var result = ChromosomeTools.AnalyzeTelomeres("chr1", Repeat("ACGT", 300));

        Assert.Multiple(() =>
        {
            Assert.That(result.Has3PrimeTelomere, Is.False);
            Assert.That(result.TelomereLength3Prime, Is.EqualTo(0));
            Assert.That(result.Has5PrimeTelomere, Is.False);
        });
    }
}
