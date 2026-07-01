using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class DistanceMatrixTests
{
    [Test]
    public void DistanceMatrix_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() =>
            PhylogeneticsTools.DistanceMatrix(new[] { "ACGT", "TCGT" }, "Hamming"));

        // Null / empty array → ArgumentException.
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.DistanceMatrix(null!));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.DistanceMatrix(Array.Empty<string>()));

        // Unknown method → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.DistanceMatrix(new[] { "ACGT", "TCGT" }, "bogus"));

        // Unequal lengths → ArgumentException (surfaces from CalculatePairwiseDistance).
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.DistanceMatrix(new[] { "ACGT", "ACGTACGT" }, "Hamming"));
    }

    [Test]
    public void DistanceMatrix_Binding_InvokesSuccessfully()
    {
        // Hamming distances (Seqeron.Genomics.Tests PHYLO-DIST-001 S02):
        //   AAAA vs AAAC = 1 ; AAAA vs CCCC = 4 ; AAAC vs CCCC = 3.
        var result = PhylogeneticsTools.DistanceMatrix(
            new[] { "AAAA", "AAAC", "CCCC" }, "Hamming");

        var m = result.Matrix;
        Assert.Multiple(() =>
        {
            Assert.That(m.Length, Is.EqualTo(3));
            Assert.That(m[0].Length, Is.EqualTo(3));

            // Diagonal is zero.
            Assert.That(m[0][0], Is.EqualTo(0.0));
            Assert.That(m[1][1], Is.EqualTo(0.0));
            Assert.That(m[2][2], Is.EqualTo(0.0));

            // Exact off-diagonal Hamming values.
            Assert.That(m[0][1], Is.EqualTo(1.0));
            Assert.That(m[0][2], Is.EqualTo(4.0));
            Assert.That(m[1][2], Is.EqualTo(3.0));

            // Symmetry.
            Assert.That(m[1][0], Is.EqualTo(m[0][1]));
            Assert.That(m[2][0], Is.EqualTo(m[0][2]));
            Assert.That(m[2][1], Is.EqualTo(m[1][2]));
        });

        // p-distance: 1 difference in 8 sites = 0.125 (PHYLO-DIST-001 M05).
        var p = PhylogeneticsTools.DistanceMatrix(new[] { "ACGTACGT", "TCGTACGT" }, "PDistance");
        Assert.That(p.Matrix[0][1], Is.EqualTo(0.125).Within(1e-10));
    }
}
