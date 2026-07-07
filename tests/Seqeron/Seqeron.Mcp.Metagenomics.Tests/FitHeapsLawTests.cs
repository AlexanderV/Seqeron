using NUnit.Framework;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps PanGenomeAnalyzer.FitHeapsLaw (genomes overload). Model n(N) = Intercept * N^(-Alpha);
// IsOpen when Alpha < 1 (Tettelin 2008; micropan heaps). Reference behaviour from
// Seqeron.Genomics.Tests PanGenomeAnalyzer_FitHeapsLaw_Tests M2: a constant new-gene curve
// (each added genome contributes exactly one brand-new cluster) fits Alpha = 0, Intercept = 1,
// open. With permutations = 1 the fit uses the natural genome order and is exactly reproducible.
[TestFixture]
public class FitHeapsLawTests
{
    private const string Core = "ATGCGATCGATCGATCGATCGATCGATCGA";
    private const string A = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
    private const string B = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC";

    private static GenomeInput Genome(string id, params (string GeneId, string Seq)[] genes)
        => new(id, genes.Select(g => new GeneInput(g.GeneId, g.Seq)).ToList());

    // g1{core}, g2{core,a}, g3{core,b}: each later genome adds exactly one new cluster -> y=[1,1].
    private static GenomeInput[] ConstantCurveGenomes() => new[]
    {
        Genome("g1", ("c1", Core)),
        Genome("g2", ("c2", Core), ("a", A)),
        Genome("g3", ("c3", Core), ("b", B)),
    };

    [Test]
    public void FitHeapsLaw_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.FitHeapsLaw(ConstantCurveGenomes(), permutations: 1));

        // Empty genome set is defined (degenerate closed fit), not an error.
        Assert.DoesNotThrow(() =>
            MetagenomicsTools.FitHeapsLaw(System.Array.Empty<GenomeInput>()));
    }

    [Test]
    public void FitHeapsLaw_Binding_InvokesSuccessfully()
    {
        // Constant new-gene curve y=[1,1] -> Alpha = 0, Intercept = 1, open.
        var fit = MetagenomicsTools.FitHeapsLaw(ConstantCurveGenomes(), permutations: 1);

        Assert.Multiple(() =>
        {
            Assert.That(fit.Alpha, Is.EqualTo(0.0).Within(1e-9),
                "Constant new-gene count is best fit by Alpha = 0 (no decay).");
            Assert.That(fit.Intercept, Is.EqualTo(1.0).Within(1e-9),
                "Intercept equals the constant new-gene count (mean = 1).");
            Assert.That(fit.IsOpen, Is.True, "Alpha < 1 -> open pan-genome (Tettelin 2008).");
        });

        // Empty input -> degenerate closed fit (0, 0, closed).
        var empty = MetagenomicsTools.FitHeapsLaw(System.Array.Empty<GenomeInput>());
        Assert.Multiple(() =>
        {
            Assert.That(empty.Alpha, Is.EqualTo(0.0));
            Assert.That(empty.Intercept, Is.EqualTo(0.0));
            Assert.That(empty.IsOpen, Is.False);
        });
    }
}
