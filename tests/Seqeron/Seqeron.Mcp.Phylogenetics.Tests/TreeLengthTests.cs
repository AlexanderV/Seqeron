using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class TreeLengthTests
{
    [Test]
    public void TreeLength_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.TreeLength("((A:1,B:1):1,(C:1,D:1):1);"));

        // Empty / null / whitespace → ArgumentException (surfaces from ParseNewick).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeLength(""));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeLength(null!));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeLength("   "));
    }

    [Test]
    public void TreeLength_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Balanced tree: 4 leaf edges (1 each) + 2 internal (1 each) = 6
            // (Seqeron.Genomics.Tests PHYLO-STATS-001 M4).
            Assert.That(PhylogeneticsTools.TreeLength("((A:1,B:1):1,(C:1,D:1):1);").Length,
                Is.EqualTo(6.0).Within(1e-10));

            // Caterpillar: 4 leaf edges (1 each) + 2 internal (0.5 each) = 5 (PHYLO-STATS-001 M5).
            Assert.That(PhylogeneticsTools.TreeLength("(A:1,(B:1,(C:1,D:1):0.5):0.5);").Length,
                Is.EqualTo(5.0).Within(1e-10));

            // Subtree with a root edge: 1 + 1 + 0.5 = 2.5 (root edge counted) (PHYLO-STATS-001 M6).
            Assert.That(PhylogeneticsTools.TreeLength("(C:1,D:1):0.5;").Length,
                Is.EqualTo(2.5).Within(1e-10));

            // Polytomy: 0.1+0.2+0.3+0.0+0.4+0.5+0.6 = 2.1 (PHYLO-COMP-001 STATS-POLYTOMY).
            Assert.That(PhylogeneticsTools.TreeLength("((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6);").Length,
                Is.EqualTo(2.1).Within(1e-9));
        });
    }
}
