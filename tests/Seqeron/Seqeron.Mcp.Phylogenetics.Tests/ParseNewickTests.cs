using System;
using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class ParseNewickTests
{
    [Test]
    public void ParseNewick_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.ParseNewick("((A,B),(C,D));"));

        // Empty / null / whitespace → ArgumentException (Seqeron.Genomics.Tests PHYLO-NEWICK-001).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.ParseNewick(""));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.ParseNewick(null!));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.ParseNewick("   "));

        // Malformed (unbalanced parens / trailing garbage) → FormatException.
        Assert.Throws<FormatException>(() => PhylogeneticsTools.ParseNewick("(A,B"));
        Assert.Throws<FormatException>(() => PhylogeneticsTools.ParseNewick("(A,B);extra"));
    }

    [Test]
    public void ParseNewick_Binding_InvokesSuccessfully()
    {
        // Nested balanced 4-taxon tree (PHYLO-NEWICK-001 ParseNewick_LeafCountMatchesInput,
        // PHYLO-STATS-001 depth of ((A,B),(C,D)) = 2).
        var r = PhylogeneticsTools.ParseNewick("((A:1,B:1):1,(C:1,D:1):1);");
        Assert.Multiple(() =>
        {
            Assert.That(r.LeafCount, Is.EqualTo(4));
            Assert.That(r.Taxa, Is.EqualTo(new[] { "A", "B", "C", "D" }));  // pre-order
            Assert.That(r.Depth, Is.EqualTo(2));
            // Σ branch lengths = 1+1+1+1 (leaves) + 1+1 (internal) = 6 (PHYLO-STATS-001 M4).
            Assert.That(r.TotalLength, Is.EqualTo(6.0).Within(1e-10));
            // Canonical re-serialization round-trips the topology.
            Assert.That(r.Newick, Does.EndWith(";"));
            Assert.That(r.Newick, Does.Contain("A").And.Contain("B").And.Contain("C").And.Contain("D"));
        });

        // Two-taxon tree: 2 leaves, depth 1 (PHYLO-STATS-001 S2).
        var two = PhylogeneticsTools.ParseNewick("(A:1,B:1);");
        Assert.Multiple(() =>
        {
            Assert.That(two.LeafCount, Is.EqualTo(2));
            Assert.That(two.Depth, Is.EqualTo(1));
        });

        // Single leaf: 1 leaf, depth 0 (PHYLO-NEWICK-001 single-taxon; PHYLO-STATS-001 M9).
        var leaf = PhylogeneticsTools.ParseNewick("A;");
        Assert.Multiple(() =>
        {
            Assert.That(leaf.LeafCount, Is.EqualTo(1));
            Assert.That(leaf.Depth, Is.EqualTo(0));
            Assert.That(leaf.Taxa, Is.EqualTo(new[] { "A" }));
        });
    }
}
