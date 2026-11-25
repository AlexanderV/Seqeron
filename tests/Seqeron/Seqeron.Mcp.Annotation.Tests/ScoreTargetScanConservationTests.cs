using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ScoreTargetScanConservationTests
{
    private static readonly string[] KitSpecies = { "10090", "10116", "9031", "9606", "9615" };

    [Test]
    public void ScoreTargetScanConservation_Schema_ValidatesCorrectly()
    {
        var utrs = new Dictionary<string, string> { ["9606"] = new string('A', 30) + "CUACCUCA" + new string('A', 20) };
        Assert.DoesNotThrow(() => AnnotationTools.ScoreTargetScanConservation("GAGGUAG", utrs));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetScanConservation("", utrs));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetScanConservation("GAGGUAG", null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetScanConservation("GAGGUAG", new Dictionary<string, string>()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetScanConservation("GAGGUAG", utrs, referenceSpecies: ""));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetScanConservation("GAGGUAG", utrs, utrBin: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetScanConservation("GAGGUAG", utrs, utrBin: 11));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetScanConservation("GAGGUAG", utrs, thresholds: "Bogus"));
    }

    // Kit gene LIN28B, family let-7/98: UTRs_median_BLs_bins.txt "LIN28B 1.12468545 1" and the rows of
    // targetscan_70_output.BL_PCT.txt (Branch length score, Pct, Conserved).
    [Test]
    public void ScoreTargetScanConservation_Binding_Lin28bLet7_EqualsTargetScanBlPct()
    {
        var r = AnnotationTools.ScoreTargetScanConservation("GAGGUAG", TargetScanKitLin28b.AlignedUtrs(), KitSpecies);

        Assert.That(r.UtrBin, Is.EqualTo(1));
        Assert.That(r.MedianBranchLength!.Value, Is.EqualTo(1.12468545).Within(1e-12));
        Assert.That(r.Sites, Has.Count.EqualTo(24));

        // (group, type) → kit Branch length score, Pct (Conserved column is empty for all: every Bls is below the TargetScanHuman7 thresholds 1.8 / 2.8 / 3.6)
        var kit = new Dictionary<(int, string), (double Bl, double Pct)>
        {
            [(1, "Seed7merM8")] = (1.6379, 0.0058),
            [(2, "Seed7merM8")] = (1.6379, 0.0058),
            [(2, "Seed8mer")] = (1.1469, 0.0),
            [(3, "Seed8mer")] = (1.6379, 0.0),
            [(4, "Seed7merM8")] = (1.6379, 0.0058),
            [(5, "Seed7merA1")] = (1.5614, 0.0),
        };
        Assert.Multiple(() =>
        {
            foreach (var s in r.Sites)
            {
                var (bl, pct) = kit[(s.GroupNumber, s.Site.SiteType)];
                string id = $"group {s.GroupNumber} {s.Site.SiteType} {s.Site.Species}";
                Assert.That(s.BranchLength, Is.EqualTo(bl), id + " Branch length score");
                Assert.That(s.Pct, Is.EqualTo(pct), id + " Pct");
                Assert.That(s.Conserved, Is.False, id + " Conserved");
            }
        });
        // species list of the branch-length step for the 8mers of group 2 (kit column "9031 9606 9615")
        Assert.That(r.Sites.First(s => s.GroupNumber == 2 && s.Site.SiteType == "Seed8mer").SpeciesWithSite,
            Is.EqualTo(new[] { "9031", "9606", "9615" }));
    }

    [Test]
    public void ScoreTargetScanConservation_CallerBin_MedianIsNull()
    {
        var r = AnnotationTools.ScoreTargetScanConservation("GAGGUAG", TargetScanKitLin28b.AlignedUtrs(), KitSpecies, utrBin: 1);
        Assert.That(r.UtrBin, Is.EqualTo(1));
        Assert.That(r.MedianBranchLength, Is.Null);
        Assert.That(r.Sites.Select(s => s.BranchLength), Is.EqualTo(
            AnnotationTools.ScoreTargetScanConservation("GAGGUAG", TargetScanKitLin28b.AlignedUtrs(), KitSpecies).Sites.Select(s => s.BranchLength)));
    }

    [Test]
    public void ScoreTargetScanConservation_Thresholds_Mouse_ConservesEightmers()
    {
        // Human 7 thresholds (8mer 1.8, 7mer-m8 2.8, 7mer-1a 3.6): no LIN28B let-7 site is conserved (kit: empty column).
        // Mouse 7 (8mer 0.6, 7mer-m8 1.8, 7mer-1a 2.5): the ten 8mer sites (Bls 1.1469 / 1.6379 ≥ 0.6) are conserved.
        var human = AnnotationTools.ScoreTargetScanConservation("GAGGUAG", TargetScanKitLin28b.AlignedUtrs(), KitSpecies);
        var mouse = AnnotationTools.ScoreTargetScanConservation("GAGGUAG", TargetScanKitLin28b.AlignedUtrs(), KitSpecies, thresholds: "targetscanmouse7");
        Assert.That(human.Sites.Count(s => s.Conserved), Is.EqualTo(0));
        Assert.That(mouse.Sites.Where(s => s.Conserved).Select(s => s.Site.SiteType).Distinct(), Is.EqualTo(new[] { "Seed8mer" }));
        Assert.That(mouse.Sites.Count(s => s.Conserved), Is.EqualTo(8));
        var th = MiRnaAnalyzer.BranchLengthConservationThresholds.TargetScanMouse7;
        Assert.That(mouse.Sites.Select(s => s.Conserved), Is.EqualTo(
            human.Sites.Select(s => MiRnaAnalyzer.IsConservedByBranchLength(s.BranchLength, Enum.Parse<MiRnaAnalyzer.TargetSiteType>(s.Site.SiteType), th))));
    }

    [Test]
    public void ScoreTargetScanConservation_DelegatesToLibrary()
    {
        var utrs = TargetScanKitLin28b.AlignedUtrs();
        var direct = MiRnaAnalyzer.ScoreTargetScanConservation("GAGGUAG", utrs, KitSpecies);
        var r = AnnotationTools.ScoreTargetScanConservation("GAGGUAG", utrs, KitSpecies);
        Assert.Multiple(() =>
        {
            Assert.That(r.UtrBin, Is.EqualTo(direct.UtrBin.Bin));
            Assert.That(r.MedianBranchLength, Is.EqualTo(direct.UtrBin.MedianBranchLength));
            Assert.That(r.Groups.Select(g => g.GroupType), Is.EqualTo(direct.Groups.Select(g => g.GroupType)));
            Assert.That(r.Sites.Select(s => (s.Site.Species, s.Site.MsaStart, s.GroupNumber, s.BranchLength, s.Pct, s.Conserved)),
                Is.EqualTo(direct.Sites.Select(s => (s.Site.Species, s.Site.MsaStart, s.GroupNumber, s.BranchLength, s.Pct, s.Conserved))));
        });
    }
}
