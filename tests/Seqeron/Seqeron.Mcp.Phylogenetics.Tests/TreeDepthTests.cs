using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class TreeDepthTests
{
    [Test]
    public void TreeDepth_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.TreeDepth("((A,B),(C,D));"));

        // Empty / null / whitespace → ArgumentException (surfaces from ParseNewick).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeDepth(""));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeDepth(null!));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeDepth("   "));
    }

    [Test]
    public void TreeDepth_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Balanced 4-taxon tree ((A,B),(C,D)): height 2 edges (Seqeron.Genomics.Tests PHYLO-STATS-001 M7).
            Assert.That(PhylogeneticsTools.TreeDepth("((A:1,B:1):1,(C:1,D:1):1);").Depth, Is.EqualTo(2));

            // Caterpillar (A,(B,(C,D))): height 3 edges (PHYLO-STATS-001 M8).
            Assert.That(PhylogeneticsTools.TreeDepth("(A:1,(B:1,(C:1,D:1):0.5):0.5);").Depth, Is.EqualTo(3));

            // Two-leaf tree (A,B): height 1 edge (PHYLO-STATS-001 S2).
            Assert.That(PhylogeneticsTools.TreeDepth("(A:1,B:1);").Depth, Is.EqualTo(1));

            // Single leaf: height 0 (PHYLO-STATS-001 M9).
            Assert.That(PhylogeneticsTools.TreeDepth("A;").Depth, Is.EqualTo(0));
        });
    }
}
