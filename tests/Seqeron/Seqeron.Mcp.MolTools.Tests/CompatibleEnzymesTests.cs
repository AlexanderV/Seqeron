using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class CompatibleEnzymesTests
{
    [Test]
    public void CompatibleEnzymes_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.compatible_enzymes());
        var result = MolToolsTools.compatible_enzymes();
        Assert.That(result, Is.Not.Null);
        Assert.That(result.Pairs, Is.Not.Null);
    }

    [Test]
    public void CompatibleEnzymes_Binding_InvokesSuccessfully()
    {
        var result = MolToolsTools.compatible_enzymes();

        Assert.Multiple(() =>
        {
            // The 10 blunt cutters are pairwise compatible ("blunt"): C(10,2) = 45 blunt pairs.
            var bluntPairs = result.Pairs.Where(p => p.CompatibleEnd == "blunt").ToList();
            Assert.That(bluntPairs, Has.Count.EqualTo(45));

            // Every reported pair carries a non-empty compatible end and two distinct enzymes.
            Assert.That(result.Pairs.All(p => !string.IsNullOrEmpty(p.CompatibleEnd)), Is.True);
            Assert.That(result.Pairs.All(p => p.Enzyme1 != p.Enzyme2), Is.True);

            // BamHI (GGATCC,1,5) and BglII (AGATCT,1,5) both leave a 5' "GATC" overhang -> compatible.
            bool bamBgl = result.Pairs.Any(p =>
                p.CompatibleEnd == "GATC" &&
                new[] { p.Enzyme1, p.Enzyme2 }.Contains("BamHI") &&
                new[] { p.Enzyme1, p.Enzyme2 }.Contains("BglII"));
            Assert.That(bamBgl, Is.True, "BamHI/BglII should be GATC-compatible");

            // The complete listing, cross-checked against Biopython 1.88's `%` operator over all
            // 1521 ordered pairs of the 39 built-ins (0 disagreements): 64 pairs in 5 end groups.
            Assert.That(result.Pairs, Has.Count.EqualTo(64));
            Assert.That(result.Pairs.Count(p => p.CompatibleEnd == "GATC"), Is.EqualTo(6));
            Assert.That(result.Pairs.Count(p => p.CompatibleEnd == "CTAG"), Is.EqualTo(6));
            Assert.That(result.Pairs.Count(p => p.CompatibleEnd == "CG"), Is.EqualTo(6));
            Assert.That(result.Pairs.Count(p => p.CompatibleEnd == "TCGA"), Is.EqualTo(1));
            Assert.That(result.Pairs.Select(p => p.CompatibleEnd).Distinct().Count(), Is.EqualTo(5));
        });
    }

    /// <summary>
    /// RESTR-FILTER-001 / F32: with `enzyme_name` the tool answers "what can THIS enzyme be ligated
    /// to", the enzyme itself included — the query an enzyme with a unique overhang needs, since 15
    /// of the 39 built-ins appear in no pair of the unfiltered listing. Expected sets are
    /// Biopython 1.88's symmetrised `%` over the same 39 enzymes.
    /// </summary>
    [Test]
    public void CompatibleEnzymes_PerEnzyme_ReturnsTheEnzymesOwnCompatibleSet()
    {
        var bamHI = MolToolsTools.compatible_enzymes("BamHI");
        var sfiI = MolToolsTools.compatible_enzymes("SfiI");
        var ecoRV = MolToolsTools.compatible_enzymes("ecorv");   // case-insensitive

        Assert.Multiple(() =>
        {
            Assert.That(bamHI.Pairs.Select(p => p.Enzyme2),
                Is.EquivalentTo(new[] { "BamHI", "BglII", "MboI", "Sau3AI" }));
            Assert.That(bamHI.Pairs.All(p => p.Enzyme1 == "BamHI" && p.CompatibleEnd == "GATC"), Is.True);

            // SfiI's 3' NNN overhang is unique in this table, so it is compatible with itself only
            // and is absent from the unfiltered listing entirely.
            Assert.That(sfiI.Pairs.Select(p => p.Enzyme2), Is.EquivalentTo(new[] { "SfiI" }));
            Assert.That(sfiI.Pairs.Single().CompatibleEnd, Is.EqualTo("NNN"));
            Assert.That(MolToolsTools.compatible_enzymes().Pairs
                .Any(p => p.Enzyme1 == "SfiI" || p.Enzyme2 == "SfiI"), Is.False);

            // Every blunt cutter is compatible with every blunt cutter, itself included.
            Assert.That(ecoRV.Pairs.Select(p => p.Enzyme2), Is.EquivalentTo(new[]
            {
                "AluI", "DpnI", "EcoRV", "HaeIII", "HincII", "RsaI", "ScaI", "SmaI", "StuI", "SwaI"
            }));
            Assert.That(ecoRV.Pairs.All(p => p.CompatibleEnd == "blunt"), Is.True);

            Assert.Throws<ArgumentException>(() => MolToolsTools.compatible_enzymes("NotAnEnzyme"));
        });
    }
}
