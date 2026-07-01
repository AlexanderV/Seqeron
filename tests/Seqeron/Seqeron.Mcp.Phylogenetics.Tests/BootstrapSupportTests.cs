using NUnit.Framework;
using Seqeron.Mcp.Phylogenetics.Tools;

namespace Seqeron.Mcp.Phylogenetics.Tests;

[TestFixture]
public class BootstrapSupportTests
{
    // Two well-separated invariant groups: A=B (all A), C=D (all G). Distances are invariant
    // under column resampling, so every UPGMA replicate recovers the {A,B},{C,D} topology.
    // Source: Seqeron.Genomics.Tests PHYLO-BOOT-001 (PhylogeneticAnalyzer_Bootstrap_Tests, M1).
    private static Dictionary<string, string> TwoGroupAlignment() => new()
    {
        ["A"] = "AAAAAAAAAA",
        ["B"] = "AAAAAAAAAA",
        ["C"] = "GGGGGGGGGG",
        ["D"] = "GGGGGGGGGG",
    };

    [Test]
    public void BootstrapSupport_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => PhylogeneticsTools.BootstrapSupport(TwoGroupAlignment(), replicates: 20));

        // Null sequences → ArgumentException (ArgumentNullException derives from it).
        Assert.Throws<ArgumentException>(() => PhylogeneticsTools.BootstrapSupport(null!));

        // Fewer than 2 sequences → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BootstrapSupport(new Dictionary<string, string> { ["A"] = "ACGT" }));

        // replicates < 1 → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BootstrapSupport(TwoGroupAlignment(), replicates: 0));

        // Unknown distance method → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BootstrapSupport(TwoGroupAlignment(), distanceMethod: "bogus"));

        // Unknown tree method → ArgumentException.
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticsTools.BootstrapSupport(TwoGroupAlignment(), treeMethod: "bogus"));
    }

    [Test]
    public void BootstrapSupport_Binding_InvokesSuccessfully()
    {
        // UPGMA on two invariant groups: both {A,B} and {C,D} clades recovered in every replicate.
        // Expected support exactly 1.0 (Seqeron.Genomics.Tests PHYLO-BOOT-001 M1).
        var result = PhylogeneticsTools.BootstrapSupport(
            TwoGroupAlignment(), replicates: 100, distanceMethod: "JukesCantor", treeMethod: "UPGMA");

        var byClade = result.Support.ToDictionary(x => x.Clade, x => x.Support);

        Assert.Multiple(() =>
        {
            Assert.That(byClade.ContainsKey("A|B"), Is.True, "Reference clade {A,B} must be scored");
            Assert.That(byClade.ContainsKey("C|D"), Is.True, "Reference clade {C,D} must be scored");
            Assert.That(byClade["A|B"], Is.EqualTo(1.0).Within(1e-10));
            Assert.That(byClade["C|D"], Is.EqualTo(1.0).Within(1e-10));
        });

        // Determinism: fixed internal seed (42) → identical support across invocations
        // (Seqeron.Genomics.Tests PHYLO-BOOT-001 M5).
        var second = PhylogeneticsTools.BootstrapSupport(
            TwoGroupAlignment(), replicates: 100, distanceMethod: "JukesCantor", treeMethod: "UPGMA");
        var byClade2 = second.Support.ToDictionary(x => x.Clade, x => x.Support);
        Assert.That(byClade2["A|B"], Is.EqualTo(byClade["A|B"]).Within(1e-12));
        Assert.That(byClade2["C|D"], Is.EqualTo(byClade["C|D"]).Within(1e-12));

        // NeighborJoining trifurcation ((A,B),C,D): only {A,B} is a rooted clade (PHYLO-BOOT-001 M7).
        var nj = PhylogeneticsTools.BootstrapSupport(
            TwoGroupAlignment(), replicates: 50, distanceMethod: "JukesCantor", treeMethod: "NeighborJoining");
        var njByClade = nj.Support.ToDictionary(x => x.Clade, x => x.Support);
        Assert.Multiple(() =>
        {
            Assert.That(njByClade.Keys, Is.EquivalentTo(new[] { "A|B" }));
            Assert.That(njByClade["A|B"], Is.EqualTo(1.0).Within(1e-10));
        });
    }
}
