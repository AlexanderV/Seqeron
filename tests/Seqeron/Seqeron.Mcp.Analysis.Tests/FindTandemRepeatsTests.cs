using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>find_tandem_repeats</c> MCP tool.
/// Expected values from the consecutive-unit definition (unit, position, repetitions,
/// unit.Length*repetitions), NOT the wrapper output.
/// </summary>
[TestFixture]
public class FindTandemRepeatsTests
{
    [Test]
    public void FindTandemRepeats_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.FindTandemRepeats("ATGATGATG", 3, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeats("", 3, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeats(null!, 3, 2));
        Assert.Throws<ArgumentException>(() => AnalysisTools.FindTandemRepeats("XYZ", 3, 2));
        // Algorithm rejects minUnitLength < 1 and minRepetitions < 2.
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindTandemRepeats("ATGATGATG", 0, 2));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindTandemRepeats("ATGATGATG", 3, 1));
    }

    [Test]
    public void FindTandemRepeats_Binding_InvokesSuccessfully()
    {
        // "ATGATGATG" unit>=3, rep>=2 -> ATG at 0, 3 copies, span 9.
        var tri = AnalysisTools.FindTandemRepeats("ATGATGATG", 3, 2).Items;
        Assert.Multiple(() =>
        {
            Assert.That(tri, Has.Length.EqualTo(1));
            Assert.That(tri[0].Unit, Is.EqualTo("ATG"));
            Assert.That(tri[0].Position, Is.EqualTo(0));
            Assert.That(tri[0].Repetitions, Is.EqualTo(3));
            Assert.That(tri[0].TotalLength, Is.EqualTo(9));
        });

        // "AAAAA" unit>=1 -> A x5 (span 5) and, since unit length 2 is also scanned,
        // the AA x2 tiling (span 4).
        var mono = AnalysisTools.FindTandemRepeats("AAAAA", 1, 2).Items;
        var a = mono.Single(i => i.Unit == "A");
        var aa = mono.Single(i => i.Unit == "AA");
        Assert.Multiple(() =>
        {
            Assert.That(a.Repetitions, Is.EqualTo(5));
            Assert.That(a.TotalLength, Is.EqualTo(5));
            Assert.That(aa.Repetitions, Is.EqualTo(2));
            Assert.That(aa.TotalLength, Is.EqualTo(4));
        });
    }
}
