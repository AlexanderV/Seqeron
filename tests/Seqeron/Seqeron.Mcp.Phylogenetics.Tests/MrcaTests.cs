using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class MrcaTests
{
    private const string Tree = "((A,B),(C,D));";

    [Test]
    public void Mrca_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.Mrca(Tree, "A", "B"));

        // Empty / null / whitespace Newick → ArgumentException (surfaces from ParseNewick).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.Mrca("", "A", "B"));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.Mrca(null!, "A", "B"));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.Mrca("   ", "A", "B"));
    }

    [Test]
    public void Mrca_Binding_InvokesSuccessfully()
    {
        // Siblings A,B in ((A,B),(C,D)): MRCA subtree taxa are exactly {A,B}
        // (Seqeron.Genomics.Tests PHYLO-COMP-001 MRCA-M02).
        var ab = PhylogeneticsTools.Mrca(Tree, "A", "B");
        Assert.Multiple(() =>
        {
            Assert.That(ab.Found, Is.True);
            Assert.That(ab.Taxa, Is.EquivalentTo(new[] { "A", "B" }));
        });

        // Cross-clade A,C: MRCA is the root, containing all four taxa (MRCA-M03 / MRCA-S01).
        var ac = PhylogeneticsTools.Mrca(Tree, "A", "C");
        Assert.Multiple(() =>
        {
            Assert.That(ac.Found, Is.True);
            Assert.That(ac.Taxa, Is.EquivalentTo(new[] { "A", "B", "C", "D" }));
        });

        // Self-MRCA of an existing leaf is that leaf (MRCA-M01).
        var self = PhylogeneticsTools.Mrca(Tree, "A", "A");
        Assert.Multiple(() =>
        {
            Assert.That(self.Found, Is.True);
            Assert.That(self.Name, Is.EqualTo("A"));
            Assert.That(self.Taxa, Is.EquivalentTo(new[] { "A" }));
        });

        // Non-existent taxon → found=false with empty fields (MRCA-M06 / binding contract).
        var missing = PhylogeneticsTools.Mrca(Tree, "A", "Z");
        Assert.Multiple(() =>
        {
            Assert.That(missing.Found, Is.False);
            Assert.That(missing.Name, Is.Empty);
            Assert.That(missing.Taxa, Is.Empty);
        });
    }
}
