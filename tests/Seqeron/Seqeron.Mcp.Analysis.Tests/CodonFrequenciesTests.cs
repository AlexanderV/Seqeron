using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>codon_frequencies</c> MCP tool.
/// Expected values taken from SEQ-CODON-FREQ-001 / the algorithm's own unit tests
/// (SequenceStatistics_CalculateCodonFrequencies_Tests, Kazusa CUTG count/total),
/// NOT the wrapper's output.
/// </summary>
[TestFixture]
public class CodonFrequenciesTests
{
    [Test]
    public void CodonFrequencies_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.CodonFrequencies("ATGATGAAA"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CodonFrequencies(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CodonFrequencies(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CodonFrequencies("ATGATGAAA", 3));
        Assert.Throws<ArgumentException>(() => AnalysisTools.CodonFrequencies("ATGATGAAA", -1));
    }

    [Test]
    public void CodonFrequencies_Binding_InvokesSuccessfully()
    {
        // Frame 0 of ATGATGAAA: ATG, ATG, AAA -> ATG=2/3, AAA=1/3.
        var f0 = AnalysisTools.CodonFrequencies("ATGATGAAA", 0).Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(f0["ATG"], Is.EqualTo(2.0 / 3.0).Within(1e-10));
            Assert.That(f0["AAA"], Is.EqualTo(1.0 / 3.0).Within(1e-10));
            Assert.That(f0.Keys, Is.EquivalentTo(new[] { "ATG", "AAA" }));
        });
    }

    [Test]
    public void CodonFrequencies_FrameShiftsMultiset()
    {
        // Frame 1: TGA, TGA (trailing AA ignored) -> TGA=1.0.
        var f1 = AnalysisTools.CodonFrequencies("ATGATGAAA", 1).Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(f1["TGA"], Is.EqualTo(1.0).Within(1e-10));
            Assert.That(f1.Keys, Is.EquivalentTo(new[] { "TGA" }));
        });

        // Frame 2: GAT, GAA (trailing A ignored) -> each 1/2.
        var f2 = AnalysisTools.CodonFrequencies("ATGATGAAA", 2).Frequencies;
        Assert.Multiple(() =>
        {
            Assert.That(f2["GAT"], Is.EqualTo(0.5).Within(1e-10));
            Assert.That(f2["GAA"], Is.EqualTo(0.5).Within(1e-10));
            Assert.That(f2.Keys, Is.EquivalentTo(new[] { "GAT", "GAA" }));
        });
    }
}
