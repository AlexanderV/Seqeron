using NUnit.Framework;
using Seqeron.Genomics.Metagenomics;
using Seqeron.Mcp.Metagenomics.Tools;

namespace Seqeron.Mcp.Metagenomics.Tests;

// Wraps MetagenomicsAnalyzer.GenerateTaxonomicProfile: per-rank abundances + species-level
// Shannon/Simpson; unclassified reads excluded from the classified count and denominator.
// Reference values from Seqeron.Genomics.Tests MetagenomicsAnalyzer_TaxonomicProfile_Tests.
[TestFixture]
public class TaxonomicProfileTests
{
    private const double Ln2 = 0.6931471805599453;

    private static MetagenomicsAnalyzer.TaxonomicClassification Classified(string readId, string species)
        => new(readId, TaxonId: 100, TaxonName: species, Rank: "species", RtlScore: 100,
            Confidence: 0.9, MatchedKmers: 100, TotalKmers: 110,
            Kingdom: "Bacteria", Phylum: "Proteobacteria", Class: "", Order: "", Family: "",
            Genus: "g", Species: species);

    private static MetagenomicsAnalyzer.TaxonomicClassification Unclassified(string readId)
        => new(readId, TaxonId: 1, TaxonName: "root", Rank: "root", RtlScore: 0,
            Confidence: 0, MatchedKmers: 0, TotalKmers: 110,
            Kingdom: "Unclassified", Phylum: "", Class: "", Order: "", Family: "",
            Genus: "", Species: "");

    private static MetagenomicsAnalyzer.TaxonomicClassification[] Reads() => new[]
    {
        Classified("r1", "coli"),
        Classified("r2", "coli"),
        Classified("r3", "aureus"),
        Classified("r4", "aureus"),
        Unclassified("r5"),
    };

    [Test]
    public void TaxonomicProfile_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MetagenomicsTools.TaxonomicProfile(Reads()));

        // Empty input is defined (empty profile), not an error.
        Assert.DoesNotThrow(() => MetagenomicsTools.TaxonomicProfile(
            Array.Empty<MetagenomicsAnalyzer.TaxonomicClassification>()));
    }

    [Test]
    public void TaxonomicProfile_Binding_InvokesSuccessfully()
    {
        // 5 reads: 4 classified (coli x2, aureus x2), 1 unclassified.
        // Species abundances 0.5/0.5 -> Shannon = ln2, Simpson = 0.5.
        var profile = MetagenomicsTools.TaxonomicProfile(Reads());

        Assert.Multiple(() =>
        {
            Assert.That(profile.TotalReads, Is.EqualTo(5));
            Assert.That(profile.ClassifiedReads, Is.EqualTo(4),
                "Unclassified read excluded from the classified count.");
            var coli = profile.SpeciesAbundance.Single(a => a.Name == "coli");
            var aureus = profile.SpeciesAbundance.Single(a => a.Name == "aureus");
            Assert.That(coli.Fraction, Is.EqualTo(0.5).Within(1e-12),
                "2 of 4 classified reads -> 0.5 (unclassified excluded from denominator).");
            Assert.That(aureus.Fraction, Is.EqualTo(0.5).Within(1e-12));
            Assert.That(profile.ShannonDiversity, Is.EqualTo(Ln2).Within(1e-10),
                "Two even species -> Shannon = ln2.");
            Assert.That(profile.SimpsonDiversity, Is.EqualTo(0.5).Within(1e-10),
                "Two even species -> Simpson = 0.5.");
            var bacteria = profile.KingdomAbundance.Single(a => a.Name == "Bacteria");
            Assert.That(bacteria.Fraction, Is.EqualTo(1.0).Within(1e-12),
                "All classified reads are Bacteria.");
        });
    }
}
