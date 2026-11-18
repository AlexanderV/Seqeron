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

        // "AAAAA" unit>=1 -> only the primitive unit A x5 (span 5): AA x2 is a non-primitive tiling of the same run
        // and is not reported (maximal-run / primitive-unit contract; Kolpakov & Kucherov 1999, MISA).
        var mono = AnalysisTools.FindTandemRepeats("AAAAA", 1, 2).Items;
        Assert.Multiple(() =>
        {
            Assert.That(mono, Has.Length.EqualTo(1));
            Assert.That(mono[0].Unit, Is.EqualTo("A"));
            Assert.That(mono[0].Position, Is.EqualTo(0));
            Assert.That(mono[0].Repetitions, Is.EqualTo(5));
            Assert.That(mono[0].TotalLength, Is.EqualTo(5));
        });

        // "ATATATAT": AT x4 @0 only (B04 / DUP_MAP section 17 request: not also ATAT x2).
        var at = AnalysisTools.FindTandemRepeats("ATATATAT", 1, 2).Items;
        Assert.That(at.Select(i => (i.Unit, i.Position, i.Repetitions)), Is.EqualTo(new[] { ("AT", 0, 4) }));
    }

    [Test]
    public void FindTandemRepeats_MaxUnitLength_RestrictsUnitLength()
    {
        // Oracle (scripts/analysis/generate_tandem_oracle.py): ACG x4 @0, GT x2 @11, TTAGC x3 @14, T x2 @14/19/24.
        const string seq = "ACGACGACGACGTGTTAGCTTAGCTTAGC";
        var all = AnalysisTools.FindTandemRepeats(seq, 1, 2).Items;
        Assert.That(all.Select(i => $"{i.Unit}:{i.Position}:{i.Repetitions}"),
            Is.EqualTo(new[] { "T:14:2", "T:19:2", "T:24:2", "GT:11:2", "ACG:0:4", "TTAGC:14:3" }));
        var upTo3 = AnalysisTools.FindTandemRepeats(seq, 1, 2, 3).Items;
        Assert.That(upTo3.Select(i => i.Unit), Is.EqualTo(new[] { "T", "T", "T", "GT", "ACG" }));
        var only5 = AnalysisTools.FindTandemRepeats(seq, 5, 2, 5).Items;
        Assert.That(only5.Select(i => $"{i.Unit}:{i.Position}:{i.Repetitions}:{i.TotalLength}"), Is.EqualTo(new[] { "TTAGC:14:3:15" }));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.FindTandemRepeats(seq, 3, 2, 2));
    }
}
