using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class ToNewickTests
{
    [Test]
    public void ToNewick_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.ToNewick("(A:0.1,B:0.2);"));

        // Empty / null / whitespace → ArgumentException (surfaces from ParseNewick).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.ToNewick(""));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.ToNewick(null!));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.ToNewick("   "));
    }

    [Test]
    public void ToNewick_Binding_InvokesSuccessfully()
    {
        // Multifurcation with branch lengths serialises with F4 formatting, InvariantCulture
        // (Seqeron.Genomics.Tests PHYLO-NEWICK-001 ParseNewick_MultifurcationWithBranchLengths).
        var withLengths = PhylogeneticsTools.ToNewick("(A:0.1,B:0.2,C:0.3);", includeBranchLengths: true);
        Assert.That(withLengths.Newick, Is.EqualTo("(A:0.1000,B:0.2000,C:0.3000);"));

        // Without branch lengths: no colons, topology preserved (PHYLO-NEWICK-001 ToNewick_WithoutBranchLengths).
        var noLengths = PhylogeneticsTools.ToNewick("(A:0.1,B:0.2,C:0.3);", includeBranchLengths: false);
        Assert.Multiple(() =>
        {
            Assert.That(noLengths.Newick, Is.EqualTo("(A,B,C);"));
            Assert.That(noLengths.Newick, Does.Not.Contain(":"));
        });

        // Nested tree ends with ';' and contains all leaf names.
        var nested = PhylogeneticsTools.ToNewick("((A,B),(C,D));", includeBranchLengths: false);
        Assert.Multiple(() =>
        {
            Assert.That(nested.Newick, Is.EqualTo("((A,B),(C,D));"));
            Assert.That(nested.Newick, Does.EndWith(";"));
        });
    }
}
