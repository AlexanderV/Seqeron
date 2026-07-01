using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class RobinsonFouldsDistanceTests
{
    [Test]
    public void RobinsonFouldsDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() =>
            PhylogeneticsTools.RobinsonFouldsDistance("((A,B),(C,D));", "((A,B),(C,D));"));

        // Empty / null / whitespace Newick → ArgumentException (surfaces from ParseNewick).
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.RobinsonFouldsDistance("", "((A,B),(C,D));"));
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.RobinsonFouldsDistance("((A,B),(C,D));", null!));
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.RobinsonFouldsDistance("   ", "((A,B),(C,D));"));
    }

    [Test]
    public void RobinsonFouldsDistance_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Identical trees → RF = 0 (Seqeron.Genomics.Tests PHYLO-COMP-001 RF-M01).
            Assert.That(PhylogeneticsTools.RobinsonFouldsDistance(
                "((A,B),(C,D));", "((A,B),(C,D));").Distance, Is.EqualTo(0));

            // Three-taxa different topology ((A,B),C) vs ((A,C),B) → RF = 2 (RF-S02).
            Assert.That(PhylogeneticsTools.RobinsonFouldsDistance(
                "((A,B),C);", "((A,C),B);").Distance, Is.EqualTo(2));

            // Four-taxa maximally different ((A,B),(C,D)) vs ((A,C),(B,D)) → RF = 4 (RF-M06).
            Assert.That(PhylogeneticsTools.RobinsonFouldsDistance(
                "((A,B),(C,D));", "((A,C),(B,D));").Distance, Is.EqualTo(4));

            // Symmetry: RF(T1,T2) = RF(T2,T1) (RF-M02).
            Assert.That(PhylogeneticsTools.RobinsonFouldsDistance(
                "((A,C),(B,D));", "((A,B),(C,D));").Distance, Is.EqualTo(4));

            // Collapsing an internal edge removes exactly one clade → RF = 1 (RF-MULTI-ROOTED).
            Assert.That(PhylogeneticsTools.RobinsonFouldsDistance(
                "(((A,B),C),(D,E));", "((A,B,C),(D,E));").Distance, Is.EqualTo(1));
        });
    }
}
