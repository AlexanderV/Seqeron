using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Phylogenetics;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// MIRNA-PCT-001 parity with TargetScan 7.0 <c>targetscan_70_BL_PCT.pl</c>, run with perl 5.38.2 +
/// BioPerl on the official BL/PCT kit (targetscan.org vert_70 <c>targetscan_70_BL_PCT.zip</c>: the
/// script reproduces the shipped <c>targetscan_70_output.BL_PCT.txt</c> byte for byte) and on small
/// constructed inputs. Every expected value below is a number the perl printed.
/// </summary>
[TestFixture]
[Category("MiRNA")]
public class MiRnaAnalyzer_PctTargetScanParity_Tests
{
    // PCT_parameters/Tree.bin_01.txt from the kit (UTR conservation bin 1; 84 species as NCBI
    // taxonomy ids). The file is multi-line Newick; whitespace is removed here because
    // PhylogeneticAnalyzer.ParseNewick rejects whitespace in front of a '(' (e.g. "(A:1, (B:2,C:3):1);"
    // throws FormatException; Biopython reads it; cross-batch request to B16).
    private const string TreeBin01 =
        "((((((((((((((((9606:0.004011,9598:0.00548785):0.00421085,9595:0.00420561):0.00738345,9601:0.0173839" +
        "):0.00308363,61853:0.0261166):0.0153816,(((9544:0.00579397,9541:0.00135107):0.00507929,9557:0.015993" +
        "9):0.00651994,60711:0.0211527):0.0308647):0.0168375,(9483:0.0287722,27679:0.0282006):0.0432467):0.06" +
        "35518,30611:0.214672):0.0160949,246437:0.237699):0.0137463,(((43179:0.19295,(51337:0.201652,((79684:" +
        "0.0991215,(10029:0.0517547,10036:0.0606729):0.029827):0.0307678,(10090:0.07483,10116:0.0764417):0.06" +
        "59564):0.109254):0.0641651):0.027301,(10181:0.0977639,(10141:0.124063,(34839:0.0708313,10160:0.12096" +
        "3):0.0143806):0.0220224):0.126832):0.033751,(9986:0.120566,9978:0.189835):0.13925):0.0392648):0.0435" +
        "588,(((9823:0.152078,((30538:0.0194633,419612:0.018398):0.112528,((9739:0.00883887,9733:0.00977372):" +
        "0.0594695,(59538:0.0146747,(9913:0.04802,(9940:0.0124486,9925:0.0154556):0.0226778):0.000493754):0.1" +
        "07692):0.0213164):0.00176269):0.0529095,(((9796:0.087093,9807:0.0682549):0.0445918,(9685:0.112188,(9" +
        "615:0.101467,(9669:0.098561,(9646:0.0653356,(9708:0.0336535,9713:0.0280281):0.0283866):0.00470832):0" +
        ".0139606):0.0188477):0.0653555):0.00861285,((9402:0.00725208,132908:0.0144228):0.122678,((225400:0.0" +
        "247071,59463:0.0182822):0.0243191,29078:0.0385308):0.0959853):0.042831):0.00502305):0.00905478,(9365" +
        ":0.240066,(42254:0.238584,143302:0.179091):0.0255468):0.0637411):0.0392947):0.0131075,(((((9785:0.06" +
        "87812,28737:0.23315):0.00972317,127582:0.0670721):0.024011,(185453:0.147422,9371:0.206074):0.0421832" +
        "):0.002216,1230840:0.14177):0.0620492,9361:0.192446):0.0410187):0.0758266,(13616:0.158057,(9305:0.08" +
        "83645,9315:0.0883367):0.0100554):0.272713):0.00643345,9258:0.43427):0.0604324,(((((((345164:0.035750" +
        "4,8954:0.0420214):0.156104,(((59894:0.149031,((44394:0.0899586,48883:0.0726105):0.0426818,59729:0.12" +
        "6817):0.0259122):0.0201173,181119:0.153095):0.0791255,(13146:0.0896585,(241585:0.0512277,176014:0.06" +
        "24011):0.0495787):0.0964529):0.0175122):0.0116039,8932:0.209144):0.0249005,(8839:0.12677,(9031:0.126" +
        "577,9103:0.0590957):0.0661465):0.062461):0.0806016,8496:0.278286):0.0496249,((8469:0.072099,8478:0.0" +
        "749067):0.00497536,(13735:0.062758,55534:0.0698858):0.0535005):0.0647565):0.024098,28377:0.366718):0" +
        ".146068):0.0510899,8364:0.411942):0.195744,7897:0.195744);";

    // Kit rows, columns 2-5 of PCT_parameters/{8mer,7mer_m8,7mer_1a}_PCT_parameters.txt.
    private static readonly PctSigmoidParameters Let7_8mer = new(-19.5274488277304, 20.4881825909231, 1.89261216546126, 1.38265500496093);
    private static readonly PctSigmoidParameters Let7_7merM8 = new(-0.221367319203409, 1.24081845475817, 1.17370692612852, 3.41806979609985);
    private static readonly PctSigmoidParameters Let7_7merA1 = new(-19.1757185248891, 20.1324465827483, 1.20519593607371, 0.149796042394506);
    private static readonly PctSigmoidParameters MiR18_7merM8 = new(-0.0131143374708748, 0.858880428089925, 1.07700649062619, 4.22046780427304);
    private static readonly PctSigmoidParameters MiR10_7merM8 = new(-4.78533199013563e-08, 0.649311358067879, 5.83882485462565, 22.7807245754396);
    private static readonly PctSigmoidParameters MiR101_7merA1 = new(-1.63736632783097, 2.13106652078132, 0.526349738794667, -0.22849587270131);
    private static readonly PctSigmoidParameters MiR1_7merA1 = new(-21.1769568945173, 21.8595369166622, 1.68242543248019, 0.782681001915711);
    private static readonly PctSigmoidParameters MiR30_8mer = new(-1.15890235160174, 1.61944269599563, 11.8180287443514, 36.260388031391);

    private static PhylogeneticAnalyzer.PhyloNode Bin01() => PhylogeneticAnalyzer.ParseNewick(TreeBin01);

    // Rows of the kit's targetscan_70_output.BL_PCT.txt for genes in bin 1 (species list = column 13
    // of the targetscan_70.pl input, or column 12 when the group has one site type).
    private static IEnumerable<TestCaseData> KitRows()
    {
        yield return new TestCaseData("10090 10116", TargetSiteType.Seed6mer, Let7_8mer, 0.1513, 0.0)
            .SetName("Kit_CDC2L6_let7_group1_6mer");
        yield return new TestCaseData("9031 9606 9615", TargetSiteType.Seed8mer, Let7_8mer, 1.1469, 0.0)
            .SetName("Kit_LIN28B_let7_group56_8mer");
        yield return new TestCaseData("9031 9606 9615 10090 10116", TargetSiteType.Seed8mer, Let7_8mer, 1.6379, 0.0)
            .SetName("Kit_LIN28B_let7_group57_8mer");
        yield return new TestCaseData("9031 9606 9615 10090 10116", TargetSiteType.Seed7merM8, Let7_7merM8, 1.6379, 0.0058)
            .SetName("Kit_LIN28B_let7_group55_7merM8");
        yield return new TestCaseData("9031 9606 9615 10090", TargetSiteType.Seed7merA1, Let7_7merA1, 1.5614, 0.0)
            .SetName("Kit_LIN28B_let7_group59_7merA1");
        yield return new TestCaseData("10090 10116 9031 9606 9615", TargetSiteType.Seed7merM8, MiR18_7merM8, 1.6379, 0.0547)
            .SetName("Kit_CDC2L6_miR18_group17_7merM8");
        // perl prints "-0.0000" (sprintf of −4.7e-8; the "< 0" test sees the string "-0.0000" as 0).
        yield return new TestCaseData("9606 9615 10090 10116", TargetSiteType.Seed7merM8, MiR10_7merM8, 0.9265, 0.0)
            .SetName("Kit_LPHN1_miR10_group94_7merM8_NegativeZero");
        yield return new TestCaseData("9606 9615 10090 10116", TargetSiteType.Seed7merA1, MiR101_7merA1, 0.9265, 0.0)
            .SetName("Kit_CDC2L6_miR101_group12_7merA1");
        yield return new TestCaseData("9606", TargetSiteType.Seed7merA1, MiR1_7merA1, 0.0, 0.0)
            .SetName("Kit_CDC2L6_miR1_group4_7merA1_SingleSpecies");
    }

    [TestCaseSource(nameof(KitRows))]
    public void KitRow_BranchLengthAndPct_EqualTargetScan(
        string species, TargetSiteType type, PctSigmoidParameters p, double perlBl, double perlPct)
    {
        double bl = ComputeBranchLengthScore(Bin01(), species.Split(' '));
        double pct = PctFromBranchLength(bl, p, type);

        Assert.Multiple(() =>
        {
            Assert.That(bl, Is.EqualTo(perlBl), "Branch length score column (sprintf %.4f)");
            Assert.That(pct, Is.EqualTo(perlPct), "Pct column");
        });
    }

    // Constructed input, tree ((9606:1.60004,10090:1.6):0.5,(9031:1.5,9615:3.0):4.0); perl output:
    //   {9606,10090} 8mer miR-30 row   → BL 3.2000, PCT 0.1787, Conserved x
    // Unrounded, Bls = 3.20004 and the sigmoid gives 0.178813 (→ 0.1788); TargetScan rounds the BL
    // first and gets 0.178702 (→ 0.1787).
    [Test]
    public void BranchLength_RoundedTo4DecimalsBeforeSigmoid_AsTargetScan()
    {
        var tree = PhylogeneticAnalyzer.ParseNewick("((9606:1.60004,10090:1.6):0.5,(9031:1.5,9615:3.0):4.0);");

        double bl = ComputeBranchLengthScore(tree, new[] { "9606", "10090" });
        double blWide = ComputeBranchLengthScore(tree, new[] { "9606", "9031" });

        Assert.Multiple(() =>
        {
            Assert.That(bl, Is.EqualTo(3.2), "perl: 3.2000");
            Assert.That(PctFromBranchLength(bl, MiR30_8mer), Is.EqualTo(0.1787), "perl: 0.1787");
            Assert.That(IsConservedByBranchLength(bl, TargetSiteType.Seed8mer), Is.True, "perl: x (3.2 ≥ 1.8)");
            Assert.That(blWide, Is.EqualTo(7.6), "perl: 7.6000");
            Assert.That(PctFromBranchLength(blWide, MiR30_8mer), Is.EqualTo(0.4605), "perl: 0.4605");
        });
    }

    // Species not in the tree are ignored (README_BL_PCT_70: "Any other species will be ignored
    // during BL calculations"). perl run, same tree: "9606 10090 99999" → 3.2000; "9606 99999" →
    // 0.0000; "9606 9606 10090" → 3.2000. Before the fix the unmatched name still made the walk
    // count edges above the matched species (2.10004 for "9606 99999").
    [Test]
    public void SpeciesNotInTree_Ignored_DuplicatesCountOnce()
    {
        var tree = PhylogeneticAnalyzer.ParseNewick("((9606:1.60004,10090:1.6):0.5,(9031:1.5,9615:3.0):4.0);");

        Assert.Multiple(() =>
        {
            Assert.That(ComputeBranchLengthScore(tree, new[] { "9606", "10090", "99999" }), Is.EqualTo(3.2));
            Assert.That(ComputeBranchLengthScore(tree, new[] { "9606", "99999" }), Is.EqualTo(0.0));
            Assert.That(ComputeBranchLengthScore(tree, new[] { "9606", "9606", "10090" }), Is.EqualTo(3.2));
        });
    }

    // Rooting does not change the connecting-subtree length: the same unrooted tree rooted on a
    // different edge gives the same Bls (getBranchLength climbs to the nearest included node).
    [Test]
    public void BranchLength_IndependentOfRootPosition()
    {
        var rootedA = PhylogeneticAnalyzer.ParseNewick("((A:1.0,B:2.0):0.5,(C:1.5,D:3.0):4.0);");
        var rootedB = PhylogeneticAnalyzer.ParseNewick("(A:1.0,B:2.0,(C:1.5,D:3.0):4.5);");

        foreach (var set in new[] { new[] { "A", "C" }, new[] { "B", "D" }, new[] { "C", "D" }, new[] { "A", "B", "C" } })
            Assert.That(ComputeBranchLengthScore(rootedB, set), Is.EqualTo(ComputeBranchLengthScore(rootedA, set)),
                string.Join(",", set));
    }

    // Constructed input, sigmoid (0.1, 0.8, 2, 4) for which the curve at 0 is 0.1144 > 0:
    //   single-species 7mer-m8 → BL 0, PCT 0.0 (getPCT, Bls == 0)
    //   {9606,10090} 6mer      → BL 3.2000, PCT 0.0 (getPCT, 6mer)
    //   4-species 7mer-1a      → BL 12.2000, PCT 0.9000, Conserved x
    [Test]
    public void GetPctRules_ZeroBlAnd6mer_GiveZero()
    {
        var p = new PctSigmoidParameters(0.1, 0.8, 2, 4);

        Assert.Multiple(() =>
        {
            Assert.That(PctFromBranchLength(0.0, p, TargetSiteType.Seed7merM8), Is.EqualTo(0.0));
            Assert.That(PctFromBranchLength(3.2, p, TargetSiteType.Seed6mer), Is.EqualTo(0.0));
            Assert.That(PctFromBranchLength(3.2, p), Is.EqualTo(0.8335), "non-6mer path: 0.833467 → 0.8335");
            Assert.That(PctFromBranchLength(12.2, p, TargetSiteType.Seed7merA1), Is.EqualTo(0.9));
        });
    }

    // selectConservationForThisSite with %siteTypeToConsThreshold (TargetScanHuman: 1.8 / 2.8 / 3.6;
    // 6mer never). perl run: 7mer-1a at 3.2000 → "" ; 7mer-1a at 12.2000 → x; 8mer at 3.2 → x.
    [Test]
    public void IsConservedByBranchLength_UsesSiteTypeThresholds()
    {
        Assert.Multiple(() =>
        {
            Assert.That(IsConservedByBranchLength(1.8, TargetSiteType.Seed8mer), Is.True, ">= 1.8");
            Assert.That(IsConservedByBranchLength(1.7999, TargetSiteType.Seed8mer), Is.False);
            Assert.That(IsConservedByBranchLength(2.8, TargetSiteType.Seed7merM8), Is.True);
            Assert.That(IsConservedByBranchLength(2.7999, TargetSiteType.Seed7merM8), Is.False);
            Assert.That(IsConservedByBranchLength(3.2, TargetSiteType.Seed7merA1), Is.False, "perl: 7mer-1a 3.2000 not conserved");
            Assert.That(IsConservedByBranchLength(3.6, TargetSiteType.Seed7merA1), Is.True);
            Assert.That(IsConservedByBranchLength(50.0, TargetSiteType.Seed6mer), Is.False, "6mer threshold 100 / never");
            Assert.That(IsConservedByBranchLength(0.6, TargetSiteType.Seed8mer, BranchLengthConservationThresholds.TargetScanMouse7),
                Is.True, "mouse block: 8mer 0.6");
            Assert.That(IsConservedByBranchLength(2.4999, TargetSiteType.Seed7merA1, BranchLengthConservationThresholds.TargetScanMouse7),
                Is.False, "mouse block: 7mer-1a 2.5");
        });
    }

    // getPCT accepts only a numeric BL (and returns NA for a negative one); the library rejects
    // values for which no PCT is defined.
    [Test]
    public void InvalidArguments_Rejected()
    {
        var p = new PctSigmoidParameters(0, 1, 1, 0);
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => PctFromBranchLength(-0.1, p));
            Assert.Throws<ArgumentOutOfRangeException>(() => PctFromBranchLength(double.NaN, p));
            Assert.Throws<ArgumentOutOfRangeException>(() => PctFromBranchLength(double.PositiveInfinity, p));
            Assert.Throws<ArgumentException>(() => PctFromBranchLength(1.0, new PctSigmoidParameters(double.NaN, 1, 1, 0)));
            Assert.Throws<ArgumentException>(() => PctFromBranchLength(1.0, new PctSigmoidParameters(0, 1, double.PositiveInfinity, 0)));
            Assert.Throws<ArgumentOutOfRangeException>(() => IsConservedByBranchLength(double.NaN, TargetSiteType.Seed8mer));
        });
    }
}
