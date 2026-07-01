using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class PairwiseDistanceTests
{
    [Test]
    public void PairwiseDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.PairwiseDistance("ACGT", "TCGT", "Hamming"));

        // Null sequences → ArgumentException.
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.PairwiseDistance(null!, "ACGT"));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.PairwiseDistance("ACGT", null!));

        // Unequal lengths → ArgumentException (Seqeron.Genomics.Tests PHYLO-DIST-001 M10).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.PairwiseDistance("ACGT", "ACGTACGT"));

        // Unknown method → ArgumentException.
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.PairwiseDistance("ACGT", "TCGT", "bogus"));
    }

    [Test]
    public void PairwiseDistance_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Hamming: 2 raw mismatches (PHYLO-DIST-001 M06).
            Assert.That(PhylogeneticsTools.PairwiseDistance("ACGTACGT", "TCGTACGA", "Hamming").Distance,
                Is.EqualTo(2.0));

            // p-distance: 1 diff in 8 sites = 0.125 (PHYLO-DIST-001 M05).
            Assert.That(PhylogeneticsTools.PairwiseDistance("ACGTACGT", "TCGTACGT", "PDistance").Distance,
                Is.EqualTo(0.125).Within(1e-10));

            // Jukes-Cantor: p=1/8 → d = -3/4·ln(5/6) ≈ 0.13674 (PHYLO-DIST-001 M08).
            Assert.That(PhylogeneticsTools.PairwiseDistance("ACGTACGT", "TCGTACGT", "JukesCantor").Distance,
                Is.EqualTo(0.13674).Within(1e-4));

            // K2P pure transition (A→G at one of 4 sites): S=0.25, V=0 → -1/2·ln(1/2) ≈ 0.34657
            // (PHYLO-DIST-001 M09).
            Assert.That(PhylogeneticsTools.PairwiseDistance("ACGT", "GCGT", "Kimura2Parameter").Distance,
                Is.EqualTo(0.34657).Within(1e-4));

            // K2P pure transversion (A→C): S=0, V=0.25 → ≈ 0.31713 (PHYLO-DIST-001 M09).
            Assert.That(PhylogeneticsTools.PairwiseDistance("ACGT", "CCGT", "Kimura2Parameter").Distance,
                Is.EqualTo(0.31713).Within(1e-4));

            // JC saturation: all-different (p=1.0) → +Infinity (PHYLO-DIST-001 M13).
            Assert.That(double.IsPositiveInfinity(
                PhylogeneticsTools.PairwiseDistance("AAAA", "CCCC", "JukesCantor").Distance), Is.True);

            // Identical sequences → 0 for all methods (PHYLO-DIST-001 M01).
            Assert.That(PhylogeneticsTools.PairwiseDistance("ACGTACGT", "ACGTACGT", "JukesCantor").Distance,
                Is.EqualTo(0.0));
        });
    }
}
