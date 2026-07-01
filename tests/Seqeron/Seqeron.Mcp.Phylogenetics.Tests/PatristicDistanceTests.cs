using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class PatristicDistanceTests
{
    [Test]
    public void PatristicDistance_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.PatristicDistance("(A:1,B:2);", "A", "B"));

        // Empty / null / whitespace Newick → ArgumentException (surfaces from ParseNewick).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.PatristicDistance("", "A", "B"));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.PatristicDistance(null!, "A", "B"));
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.PatristicDistance("   ", "A", "B"));
    }

    [Test]
    public void PatristicDistance_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Two-taxon tree (A:1,B:2): PD(A,B) = 1 + 2 = 3 (Seqeron.Genomics.Tests PHYLO-COMP-001 PD-M02).
            Assert.That(PhylogeneticsTools.PatristicDistance("(A:1,B:2);", "A", "B").Distance,
                Is.EqualTo(3.0).Within(1e-10));

            // Same taxon → 0 (PD-M01).
            Assert.That(PhylogeneticsTools.PatristicDistance("(A:1,B:2);", "A", "A").Distance,
                Is.EqualTo(0.0));

            // Non-existent taxon → NaN (PD-M04).
            Assert.That(double.IsNaN(
                PhylogeneticsTools.PatristicDistance("(A:1,B:2);", "A", "Z").Distance), Is.True);
        });

        // Polytomy tree ((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6):
        //   PD(A,B) = 0.1 + 0.2 = 0.3 ; PD(A,C) = 0.1 + 0.3 = 0.4 (PHYLO-COMP-001 STATS-POLYTOMY).
        const string poly = "((A:0.1,B:0.2,C:0.3):0.0,(D:0.4,E:0.5):0.6);";
        Assert.Multiple(() =>
        {
            Assert.That(PhylogeneticsTools.PatristicDistance(poly, "A", "B").Distance,
                Is.EqualTo(0.3).Within(1e-9));
            Assert.That(PhylogeneticsTools.PatristicDistance(poly, "A", "C").Distance,
                Is.EqualTo(0.4).Within(1e-9));
        });
    }
}
