using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class BuildTreeFromMatrixTests
{
    // Wikipedia UPGMA 5S rRNA working example (Seqeron.Genomics.Tests PHYLO-TREE-001 S01/S01b/S01d).
    private static readonly string[] UpgmaTaxa = { "a", "b", "c", "d", "e" };
    private static double[][] UpgmaMatrix() => new[]
    {
        new double[] { 0, 17, 21, 31, 23 },
        new double[] { 17, 0, 30, 34, 21 },
        new double[] { 21, 30, 0, 28, 39 },
        new double[] { 31, 34, 28, 0, 43 },
        new double[] { 23, 21, 39, 43, 0 },
    };

    // Wikipedia Neighbor-Joining additive matrix (Seqeron.Genomics.Tests PHYLO-TREE-001 S02e).
    private static readonly string[] NjTaxa = { "a", "b", "c", "d", "e" };
    private static double[][] NjMatrix() => new[]
    {
        new double[] { 0, 5, 9, 9, 8 },
        new double[] { 5, 0, 10, 10, 9 },
        new double[] { 9, 10, 0, 8, 7 },
        new double[] { 9, 10, 8, 0, 3 },
        new double[] { 8, 9, 7, 3, 0 },
    };

    [Test]
    public void BuildTreeFromMatrix_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.BuildTreeFromMatrix(UpgmaTaxa, UpgmaMatrix(), "UPGMA"));

        // Null taxa → ArgumentException.
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.BuildTreeFromMatrix(null!, UpgmaMatrix()));

        // Single taxon → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BuildTreeFromMatrix(new[] { "a" }, new[] { new double[] { 0 } }));

        // Null matrix → ArgumentException.
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.BuildTreeFromMatrix(UpgmaTaxa, null!));

        // Matrix size mismatch → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BuildTreeFromMatrix(new[] { "a", "b" }, UpgmaMatrix()));

        // Unknown tree method → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BuildTreeFromMatrix(UpgmaTaxa, UpgmaMatrix(), "bogus"));
    }

    [Test]
    public void BuildTreeFromMatrix_Binding_InvokesSuccessfully()
    {
        // UPGMA exact ultrametric Newick (Seqeron.Genomics.Tests PHYLO-TREE-001 S01d).
        var upgma = PhylogeneticsTools.BuildTreeFromMatrix(UpgmaTaxa, UpgmaMatrix(), "UPGMA");
        Assert.Multiple(() =>
        {
            Assert.That(upgma.Method, Is.EqualTo("UPGMA"));
            Assert.That(upgma.Newick, Is.EqualTo(
                "(((a:8.5000,b:8.5000):2.5000,e:11.0000):5.5000,(c:14.0000,d:14.0000):2.5000);"));
            Assert.That(upgma.Taxa, Is.EquivalentTo(UpgmaTaxa));
        });

        // NeighborJoining exact trifurcating Newick (Seqeron.Genomics.Tests PHYLO-TREE-001 S02e).
        var nj = PhylogeneticsTools.BuildTreeFromMatrix(NjTaxa, NjMatrix(), "NeighborJoining");
        Assert.Multiple(() =>
        {
            Assert.That(nj.Method, Is.EqualTo("NeighborJoining"));
            Assert.That(nj.Newick, Is.EqualTo(
                "(((a:2.0000,b:3.0000):3.0000,c:4.0000):2.0000,d:2.0000,e:1.0000);"));
        });
    }
}
