using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class TreeLeavesTests
{
    [Test]
    public void TreeLeaves_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.TreeLeaves("((A,B),(C,D));"));

        // Empty / null / whitespace → ArgumentException (surfaces from ParseNewick).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeLeaves(""));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeLeaves(null!));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.TreeLeaves("   "));
    }

    [Test]
    public void TreeLeaves_Binding_InvokesSuccessfully()
    {
        // Balanced tree: leaves in left-to-right pre-order A,B,C,D
        // (Seqeron.Genomics.Tests PHYLO-STATS-001 M1/C1).
        var balanced = PhylogeneticsTools.TreeLeaves("((A,B),(C,D));");
        Assert.That(balanced.Leaves.Select(l => l.Name), Is.EqualTo(new[] { "A", "B", "C", "D" }));

        // Per-leaf branch lengths preserved (PHYLO-NEWICK-001 multifurcation branch lengths).
        var withLengths = PhylogeneticsTools.TreeLeaves("((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6);");
        var byName = withLengths.Leaves.ToDictionary(l => l.Name, l => l.BranchLength);
        Assert.Multiple(() =>
        {
            Assert.That(withLengths.Leaves, Has.Count.EqualTo(5));
            Assert.That(withLengths.Leaves.Select(l => l.Name), Is.EqualTo(new[] { "A", "B", "C", "D", "E" }));
            Assert.That(byName["A"], Is.EqualTo(0.1).Within(1e-9));
            Assert.That(byName["B"], Is.EqualTo(0.2).Within(1e-9));
            Assert.That(byName["C"], Is.EqualTo(0.3).Within(1e-9));
            Assert.That(byName["D"], Is.EqualTo(0.4).Within(1e-9));
            Assert.That(byName["E"], Is.EqualTo(0.5).Within(1e-9));
        });

        // Single leaf: the node itself (PHYLO-STATS-001 M3).
        var leaf = PhylogeneticsTools.TreeLeaves("X;");
        Assert.Multiple(() =>
        {
            Assert.That(leaf.Leaves, Has.Count.EqualTo(1));
            Assert.That(leaf.Leaves[0].Name, Is.EqualTo("X"));
        });
    }
}
