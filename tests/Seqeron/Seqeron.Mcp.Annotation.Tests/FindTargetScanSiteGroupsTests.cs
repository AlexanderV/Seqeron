using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class FindTargetScanSiteGroupsTests
{
    private static readonly string[] KitSpecies = { "10090", "10116", "9031", "9606", "9615" };

    [Test]
    public void FindTargetScanSiteGroups_Schema_ValidatesCorrectly()
    {
        var utrs = new Dictionary<string, string> { ["9606"] = new string('A', 30) + "CUACCUCA" + new string('A', 20) };
        Assert.DoesNotThrow(() => AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", utrs));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups("", utrs));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups(null!, utrs));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups("GAGGUA", utrs));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups("GAGGUAN", utrs));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", new Dictionary<string, string>()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", new Dictionary<string, string> { [""] = "ACGU" }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", new Dictionary<string, string> { ["9606"] = null! }));
    }

    // Kit gene LIN28B, family let-7/98 (seed GAGGUAG), 5 species: targetscan_70_output.txt rows (UTR_Sequences_sample.txt).
    // Local group numbers 1-5 = the kit's 55-59 (the script's counter continues over genes and families).
    [Test]
    public void FindTargetScanSiteGroups_Binding_Lin28bLet7_EqualsTargetScan70Output()
    {
        var r = AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", TargetScanKitLin28b.AlignedUtrs(), KitSpecies);

        Assert.That(r.Groups.Select(g => g.GroupType),
            Is.EqualTo(new[] { "7mer-m8", "7mer-m8+8mer-1a", "8mer-1a", "7mer-m8", "7mer-1a" }));
        Assert.That(r.Groups.Select(g => g.GroupNumber), Is.EqualTo(new[] { 1, 2, 3, 4, 5 }));

        // (group, species, MSA start, MSA end, UTR start, UTR end, type) of the kit rows
        var expected = new (int G, string Sp, int Ms, int Me, int Us, int Ue, string Type)[]
        {
            (1, "10090", 2922, 2928, 2362, 2368, "Seed7merM8"), (1, "10116", 2922, 2928, 2258, 2264, "Seed7merM8"),
            (1, "9031", 2922, 2928, 2255, 2261, "Seed7merM8"), (1, "9606", 2922, 2928, 2471, 2477, "Seed7merM8"),
            (1, "9615", 2922, 2928, 2499, 2505, "Seed7merM8"),
            (2, "10090", 3117, 3123, 2534, 2540, "Seed7merM8"), (2, "10116", 3117, 3123, 2403, 2409, "Seed7merM8"),
            (2, "9031", 3117, 3124, 2420, 2427, "Seed8mer"), (2, "9606", 3117, 3124, 2653, 2660, "Seed8mer"),
            (2, "9615", 3117, 3124, 2671, 2678, "Seed8mer"),
            (3, "10090", 45, 52, 36, 43, "Seed8mer"), (3, "10116", 45, 52, 35, 42, "Seed8mer"),
            (3, "9031", 45, 52, 45, 52, "Seed8mer"), (3, "9606", 45, 52, 44, 51, "Seed8mer"), (3, "9615", 45, 52, 45, 52, "Seed8mer"),
            (4, "10090", 4740, 4746, 3950, 3956, "Seed7merM8"), (4, "10116", 4740, 4746, 3816, 3822, "Seed7merM8"),
            (4, "9031", 4740, 4746, 3874, 3880, "Seed7merM8"), (4, "9606", 4740, 4746, 4092, 4098, "Seed7merM8"),
            (4, "9615", 4740, 4746, 4132, 4138, "Seed7merM8"),
            (5, "10090", 4937, 4943, 4140, 4146, "Seed7merA1"), (5, "9031", 4937, 4943, 4064, 4070, "Seed7merA1"),
            (5, "9606", 4937, 4943, 4277, 4283, "Seed7merA1"), (5, "9615", 4937, 4943, 4323, 4329, "Seed7merA1"),
        };
        var produced = r.Groups.SelectMany(g => g.Sites.Select(s => (g.GroupNumber, s.Species, s.MsaStart, s.MsaEnd, s.UtrStart, s.UtrEnd, s.SiteType))).ToList();
        Assert.That(produced, Is.EquivalentTo(expected.Select(e => (e.G, e.Sp, e.Ms, e.Me, e.Us, e.Ue, e.Type))));
        Assert.That(r.Groups.SelectMany(g => g.Sites).All(s => !s.OrfOverlap), Is.True);

        // Heterogeneous group 2 (kit: 7mer-m8+8mer-1a): species with the 8mer = 9031 9606 9615, with the 7mer-m8 = all five.
        var g2 = r.Groups[1];
        Assert.Multiple(() =>
        {
            Assert.That(g2.Species, Is.EqualTo(new[] { "9031", "9606", "9615", "10090", "10116" }));
            Assert.That(g2.SpeciesBySiteType.Single(x => x.SiteType == "Seed8mer").Species, Is.EqualTo(new[] { "9031", "9606", "9615" }));
            Assert.That(g2.SpeciesBySiteType.Single(x => x.SiteType == "Seed7merM8").Species, Is.EqualTo(new[] { "9031", "9606", "9615", "10090", "10116" }));
            Assert.That(r.Groups[4].Species, Is.EqualTo(new[] { "9031", "9606", "9615", "10090" }));
        });
    }

    [Test]
    public void FindTargetScanSiteGroups_SpeciesWithMiRna_RestrictsSearch()
    {
        var r = AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", TargetScanKitLin28b.AlignedUtrs(), new[] { "9606" });
        Assert.That(r.Groups, Has.Count.EqualTo(5));
        Assert.That(r.Groups.SelectMany(g => g.Sites).Select(s => s.Species).Distinct(), Is.EqualTo(new[] { "9606" }));
    }

    [Test]
    public void FindTargetScanSiteGroups_DelegatesToLibrary_MatureMiRnaEqualsSeed()
    {
        var utrs = TargetScanKitLin28b.AlignedUtrs();
        var direct = MiRnaAnalyzer.FindTargetScanSiteGroups("UGAGGUAGUAGGUUGUAUAGUU", utrs, KitSpecies);
        var viaTool = AnnotationTools.FindTargetScanSiteGroups("UGAGGUAGUAGGUUGUAUAGUU", utrs, KitSpecies).Groups;
        Assert.That(viaTool, Has.Count.EqualTo(direct.Count));
        for (int i = 0; i < direct.Count; i++)
        {
            Assert.That(viaTool[i].GroupNumber, Is.EqualTo(direct[i].GroupNumber));
            Assert.That(viaTool[i].GroupType, Is.EqualTo(direct[i].GroupType));
            Assert.That(viaTool[i].Species, Is.EqualTo(direct[i].Species));
            Assert.That(viaTool[i].Sites.Select(s => (s.Species, s.MsaStart, s.MsaEnd, s.UtrStart, s.UtrEnd, s.SiteType, s.OrfOverlap)),
                Is.EqualTo(direct[i].Sites.Select(s => (s.Species, s.MsaStart, s.MsaEnd, s.UtrStart, s.UtrEnd, s.SiteType.ToString(), s.OrfOverlap))));
        }
        Assert.That(AnnotationTools.FindTargetScanSiteGroups("GAGGUAG", utrs, KitSpecies).Groups.Count, Is.EqualTo(direct.Count));
    }
}
