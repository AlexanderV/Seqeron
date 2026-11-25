using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class AssignTargetScanUtrBinTests
{
    // Constructed 11-row alignment; targetscan_70_BL_bins.pl gives "GENE3 2.55195311 3"
    // (mirrors MiRnaAnalyzer_TargetScanAlignment_Tests.ConstructedAlignment_UtrBin_EqualBinsScript).
    private static Dictionary<string, string> Alignment() => new()
    {
        ["9606"] = "ACGUACGUAC-GUACGUACGU", ["9598"] = "NCG-ACNUA--GNAC-UANGU", ["9544"] = "NCGUNCGUNC-GNACGNACGN",
        ["10090"] = "NCGUA-GUAC-GUAC-UACGN", ["10116"] = "ncguacnuac-gnacguangu", ["9913"] = "NCGUACG-AC-GUANGUACGU",
        ["9615"] = "NCGUACGUNC-GUACGNACGU", ["9031"] = "NCGUACGUA--GUACGUANGU", ["8364"] = "NCGUACGUAC-GUACGUACGN",
        ["7897"] = "NCGUACGUAC--", ["99999"] = "NCGUACGUAC-GNACGUACGU",
    };

    [Test]
    public void AssignTargetScanUtrBin_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.AssignTargetScanUtrBin(Alignment()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssignTargetScanUtrBin(null!));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssignTargetScanUtrBin(new()));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssignTargetScanUtrBin(Alignment(), ""));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssignTargetScanUtrBin(new() { ["9606"] = null! }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssignTargetScanUtrBin(new() { [""] = "ACGU" }));
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssignTargetScanUtrBin(new() { ["10090"] = "ACGU" }), "no reference row");
        Assert.Throws<ArgumentException>(() => AnnotationTools.AssignTargetScanUtrBin(Alignment(), "9606", "(A:1,B:"), "bad Newick");
    }

    [Test]
    public void AssignTargetScanUtrBin_Binding_EqualsTargetScanBinsScript()
    {
        var r = AnnotationTools.AssignTargetScanUtrBin(Alignment());
        var lib = MiRnaAnalyzer.AssignTargetScanUtrBin(Alignment());
        Assert.Multiple(() =>
        {
            Assert.That(r.Bin, Is.EqualTo(3));
            Assert.That(r.MedianBranchLength!.Value.ToString("G15", System.Globalization.CultureInfo.InvariantCulture), Is.EqualTo("2.55195311"));
            Assert.That(r.Bin, Is.EqualTo(lib.Bin));
            Assert.That(r.MedianBranchLength, Is.EqualTo(lib.MedianBranchLength));
            Assert.That(r.BinThresholds, Is.EqualTo(MiRnaAnalyzer.TargetScanUtrBinThresholds));
            Assert.That(r.BinThresholds, Has.Count.EqualTo(10));
        });
    }

    [Test]
    public void AssignTargetScanUtrBin_NoReferenceNucleotide_MedianNull_Bin1()
    {
        var r = AnnotationTools.AssignTargetScanUtrBin(new() { ["9606"] = "---" });
        Assert.That(r.Bin, Is.EqualTo(1));
        Assert.That(r.MedianBranchLength, Is.Null);
    }
}
