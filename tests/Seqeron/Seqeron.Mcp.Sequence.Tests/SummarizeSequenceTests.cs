using NUnit.Framework;
using Seqeron.Mcp.Sequence.Tools;

namespace Seqeron.Mcp.Sequence.Tests;

[TestFixture]
public class SummarizeSequenceTests
{
    [Test]
    public void SummarizeSequence_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => SequenceTools.SummarizeSequence("ATGCGATCGATCG"));
        Assert.Throws<ArgumentException>(() => SequenceTools.SummarizeSequence(""));
        Assert.Throws<ArgumentException>(() => SequenceTools.SummarizeSequence(null!));
    }

    [Test]
    public void SummarizeSequence_Binding_InvokesSuccessfully()
    {
        var result = SequenceTools.SummarizeSequence("ATGCGATCGATCG");
        Assert.That(result.Length, Is.EqualTo(13));
        Assert.That(result.GcContent, Is.InRange(0, 100));
        Assert.That(result.Entropy, Is.GreaterThanOrEqualTo(0));
        Assert.That(result.Complexity, Is.InRange(0, 1));
        Assert.That(result.MeltingTemperature, Is.TypeOf<double>());
        Assert.That(result.Composition, Contains.Key("A"));
        Assert.That(result.Composition, Contains.Key("T"));
        Assert.That(result.Composition, Contains.Key("G"));
        Assert.That(result.Composition, Contains.Key("C"));
    }

    // FIN-B03 F33: MW field (Biopython 1.88 molecular_weight executed) and the opt-in IUPAC-weighted Tm
    // (Biopython 1.88 Tm_Wallace("acgtNNSWacgt-ry", strict=False) = 42.0; A/C/G/T-only default 24.0).
    [Test]
    public void SummarizeSequence_MolecularWeightAndWeightedTm_MatchBiopython()
    {
        var dna = SequenceTools.SummarizeSequence("ATGCGATCGATCG");
        var rna = SequenceTools.SummarizeSequence("AUGCAUGC");
        Assert.Multiple(() =>
        {
            Assert.That(dna.MolecularWeight, Is.EqualTo(4054.5834).Within(1e-6));
            Assert.That(rna.MolecularWeight, Is.EqualTo(2589.5321).Within(1e-6));
            Assert.That(SequenceTools.SummarizeSequence("acgtNNSWacgt-ry").MeltingTemperature, Is.EqualTo(24.0));
            Assert.That(SequenceTools.SummarizeSequence("acgtNNSWacgt-ry", weightAmbiguousBasesInTm: true).MeltingTemperature,
                Is.EqualTo(42.0).Within(1e-9));
        });
    }
}
