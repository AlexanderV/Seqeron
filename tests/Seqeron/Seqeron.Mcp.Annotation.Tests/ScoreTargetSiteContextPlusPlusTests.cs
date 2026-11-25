using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Core;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ScoreTargetSiteContextPlusPlusTests
{
    // 181-nt worked example of the TargetScan 7.0 perl run (see Seqeron.Genomics.Tests
    // MiRnaAnalyzer_ContextPlusPlusTargetScanParity_Tests): targetscan_70.pl sites 8mer 31-38, 6mer 82-87,
    // 7mer-1a 125-131 (1-based); context++ scores -0.742 / -0.118 / -0.257.
    private const string Utr =
        "GCUAAAGACAAUUACAUAACAUACACGUCACUACCUCAGCACGAAACUUGUUGGCCCAGUGUGAAUCGCUUAAGGGUUGCGUACCUCGAAGUAAG" +
        "UGUGAUGCAUACGCCUUUACUUGCUGUGGUACCUCAUCCACCCCAUCGGACUGGCAUUUUUAUUACACUCAGAAACAGAACUCGGG";
    private const string Let7a = "UGAGGUAGUAGGUUGUAUAGUU";

    private IDisposable? _mode;

    // A partial context++ score (an omitted feature) is a guarded limitation (MIRNA-TARGET-001): the library throws
    // SeqeronLimitationException outside Permissive mode. These tests exercise the real scoring under Permissive;
    // the gate itself is covered by ScoreTargetSiteContextPlusPlus_PartialScore_GatedUnlessPermissive.
    [SetUp]
    public void UsePermissiveMode() => _mode = LimitationPolicy.UsePermissive();

    [TearDown]
    public void RestoreMode() => _mode?.Dispose();

    // ORF of the perl run: 1500 nt after gap removal, one let-7 8mer (targetscan_count_8mers.pl).
    private static readonly string WorkedExampleOrf =
        "ATG" + string.Concat(Enumerable.Repeat("GCTGCA", 200)) + "--ctacctca.." + "CUACCUCG"
        + string.Concat(Enumerable.Repeat("GGCAAT", 46)) + "GCTAA";

    // type, 0-based start, end, RNAplfold l=14 value read by the perl, SPS, perl context++ score
    private static IEnumerable<TestCaseData> PerlRows()
    {
        yield return new TestCaseData("Seed8mer", 30, 37, 0.009772092, -9.25, -0.742).SetName("Perl_8mer_31_38");
        yield return new TestCaseData("Seed6mer", 81, 86, 0.008715454, -6.72, -0.118).SetName("Perl_6mer_82_87");
        yield return new TestCaseData("Seed7merA1", 124, 130, 0.0005159016, -6.72, -0.257).SetName("Perl_7merA1_125_131");
    }

    [Test]
    public void ScoreTargetSiteContextPlusPlus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus("", Let7a, 30, 37, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(null!, Let7a, 30, 37, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, "", 30, 37, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, "UGAGGUA", 30, 37, "Seed8mer"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, -1, 6, "Seed8mer"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 174, Utr.Length, "Seed8mer"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 37, 30, "Seed8mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Bogus"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "99"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Offset6mer"));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 36, "Seed8mer"), "span length must match the type");
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer", orfLength: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer", orf8mCount: -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer", plfoldUnpairedProbability: 1.5));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer", sps: double.NaN));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer", ta: double.PositiveInfinity));
    }

    [Test]
    public void ScoreTargetSiteContextPlusPlus_ConservationValidation()
    {
        var species = new[] { "9606", "10090" };
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer",
            conservation: new ContextConservationInputDto(null, null, species, null)), "neither newick nor bin");
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer",
            conservation: new ContextConservationInputDto("(9606:1,10090:1);", 1, species, null)), "both");
        Assert.Throws<ArgumentOutOfRangeException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer",
            conservation: new ContextConservationInputDto(null, 11, species, null)));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer",
            conservation: new ContextConservationInputDto(null, 1, Array.Empty<string>(), null)));
        Assert.Throws<ArgumentException>(() => AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer",
            conservation: new ContextConservationInputDto("((9606:1,10090:1", null, species, null)), "malformed Newick");
    }

    [TestCaseSource(nameof(PerlRows))]
    public void ScoreTargetSiteContextPlusPlus_Binding_ReproducesTargetScanPerl_SuppliedInputs(
        string type, int start, int end, double plfold, double sps, double perlScore)
    {
        var r = AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, start, end, type, miRnaName: "hsa-let-7a",
            sps: sps, ta: 3.393, orfLength: 1500, orf8mCount: 1, plfoldUnpairedProbability: plfold);
        Assert.That(r.Score.SiteType, Is.EqualTo(type));
        Assert.That(r.Score.ContextScore!.Value, Is.EqualTo(perlScore).Within(1e-12));
        Assert.That(r.Score.OmittedFeatures, Has.None.StartsWith("SA"));
    }

    [TestCaseSource(nameof(PerlRows))]
    public void ScoreTargetSiteContextPlusPlus_Binding_BundledTaSps_OrfSequence_ReproducesTargetScanPerl(
        string type, int start, int end, double plfold, double sps, double perlScore)
    {
        // No sps/ta/orf numbers: SPS and TA from the bundled TA_SPS_by_seed_region table, Len_ORF/ORF8m from the ORF sequence.
        var r = AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, start, end, type,
            orfSequence: WorkedExampleOrf, plfoldUnpairedProbability: plfold);
        Assert.That(r.Score.ContextScore!.Value, Is.EqualTo(perlScore).Within(1e-12));
        Assert.That(r.Score.OmittedFeatures, Has.Count.EqualTo(1));
        Assert.That(r.Score.OmittedFeatures[0], Does.StartWith("PCT"), "only PCT (no conservation) is residual");
    }

    [Test]
    public void ScoreTargetSiteContextPlusPlus_8mer_Contributions_EqualPerlColumns()
    {
        // perl columns of the 8mer 31-38 row (3 decimals): intercept -0.589, 3P -0.048, local AU -0.115, Min_dist 0.004,
        // sRNA8 0.015, Len_3UTR -0.033, SA -0.073, Len_ORF 0.082, ORF8m -0.118, Off6m -0.020, TA 0.083, SPS 0.070.
        var c = AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, 30, 37, "Seed8mer",
            orfSequence: WorkedExampleOrf, plfoldUnpairedProbability: 0.009772092).Score;
        Assert.Multiple(() =>
        {
            Assert.That(c.Intercept, Is.EqualTo(-0.589));
            Assert.That(c.ThreePrimePairingContribution, Is.EqualTo(-0.048).Within(5.01e-4));
            Assert.That(c.LocalAuContribution, Is.EqualTo(-0.115).Within(5.01e-4));
            Assert.That(c.MinDistContribution, Is.EqualTo(0.004).Within(5.01e-4));
            Assert.That(c.SRna8Contribution, Is.EqualTo(0.015).Within(5.01e-4));
            Assert.That(c.Len3UtrContribution, Is.EqualTo(-0.033).Within(5.01e-4));
            Assert.That(c.SaContribution, Is.EqualTo(-0.073).Within(5.01e-4));
            Assert.That(c.LenOrfContribution, Is.EqualTo(0.082).Within(5.01e-4));
            Assert.That(c.Orf8mContribution, Is.EqualTo(-0.118).Within(5.01e-4));
            Assert.That(c.Off6mContribution, Is.EqualTo(-0.020).Within(5.01e-4));
            Assert.That(c.TaContribution, Is.EqualTo(0.083).Within(5.01e-4));
            Assert.That(c.SpsContribution, Is.EqualTo(0.070).Within(5.01e-4));
        });
    }

    [Test]
    public void ScoreTargetSiteContextPlusPlus_PartialScore_GatedUnlessPermissive()
    {
        string utr = new string('G', 60) + "CUACCUC" + "G" + new string('G', 60);
        var species = new[] { "9031", "9606", "9615", "10090", "10116" };
        var conservation = new ContextConservationInputDto(null, 1, species, null);
        using (LimitationPolicy.Use(LimitationMode.Moderate))
        {
            // SA estimate / ORF features / PCT missing → partial score → guarded limitation MIRNA-TARGET-001.
            Assert.Throws<SeqeronLimitationException>(() =>
                AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 60, 66, "Seed7merM8"));
            // Every optional feature supplied (ORF, RNAplfold value, conservation) → complete score, no gate.
            var full = AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 60, 66, "Seed7merM8",
                orfSequence: WorkedExampleOrf, plfoldUnpairedProbability: 0.05, conservation: conservation).Score;
            Assert.That(full.OmittedFeatures, Is.Empty);
            Assert.That(full.ContextScore, Is.Not.Null);
        }
        using (LimitationPolicy.UsePermissive())
            Assert.That(AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 60, 66, "Seed7merM8").Score.OmittedFeatures, Is.Not.Empty);
    }

    [Test]
    public void ScoreTargetSiteContextPlusPlus_NearCds_ContextScoreNull()
    {
        string utr = "GGGGGGG" + "CUACCUC" + "A" + new string('G', 80);   // 1-based start 8 < 15
        var r = AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 7, 14, "seed8mer");
        Assert.That(r.Score.ContextScore, Is.Null);
    }

    [Test]
    public void ScoreTargetSiteContextPlusPlus_Conservation_KitBranchLengthAndPct()
    {
        // Kit LIN28B let-7 7mer-m8 group (bin 1, five species): Branch length 1.6379, Pct 0.0058
        // (targetscan_70_output.BL_PCT.txt). A 7mer-m8 site: reverse complement of GAGGUAG = CUACCUC, not followed by A.
        string utr = new string('G', 60) + "CUACCUC" + "G" + new string('G', 60);
        var species = new[] { "9031", "9606", "9615", "10090", "10116" };
        var r = AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 60, 66, "Seed7merM8",
            conservation: new ContextConservationInputDto(null, 1, species, null)).Score;
        Assert.Multiple(() =>
        {
            Assert.That(r.BranchLengthScore, Is.EqualTo(1.6379));
            Assert.That(r.Pct, Is.EqualTo(0.0058));
            Assert.That(r.PctIsNa, Is.False);
            Assert.That(r.OmittedFeatures, Has.None.StartsWith("PCT"));
        });

        // Explicit kit row (7mer_m8_PCT_parameters.txt, let-7) gives the same.
        var explicitRow = AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 60, 66, "Seed7merM8",
            conservation: new ContextConservationInputDto(null, 1, species,
                new PctSigmoidParametersDto(-0.221367319203409, 1.24081845475817, 1.17370692612852, 3.41806979609985))).Score;
        Assert.That(explicitRow.Pct, Is.EqualTo(0.0058));
        Assert.That(explicitRow.PctContribution, Is.EqualTo(r.PctContribution));
    }

    [Test]
    public void ScoreTargetSiteContextPlusPlus_Conservation_NewickEqualsBundledBinTree()
    {
        string utr = new string('G', 60) + "CUACCUC" + "G" + new string('G', 60);
        var species = new[] { "9606", "10090" };
        var viaBin = AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 60, 66, "Seed7merM8",
            conservation: new ContextConservationInputDto(null, 3, species, null)).Score;
        var viaNewick = AnnotationTools.ScoreTargetSiteContextPlusPlus(utr, Let7a, 60, 66, "Seed7merM8",
            conservation: new ContextConservationInputDto(MiRnaAnalyzer.GetTargetScanBinTreeNewick(3), null, species, null)).Score;
        Assert.That(viaNewick.BranchLengthScore, Is.EqualTo(viaBin.BranchLengthScore));
        Assert.That(viaNewick.Pct, Is.EqualTo(viaBin.Pct));
        Assert.That(viaBin.BranchLengthScore, Is.GreaterThan(0));
    }

    [TestCaseSource(nameof(PerlRows))]
    public void ScoreTargetSiteContextPlusPlus_DelegatesToLibrary(
        string type, int start, int end, double plfold, double sps, double perlScore)
    {
        var miRna = MiRnaAnalyzer.CreateMiRna("miRNA", Let7a);
        var enumType = Enum.Parse<MiRnaAnalyzer.TargetSiteType>(type);
        var site = MiRnaAnalyzer.FindTargetSites(Utr, miRna, minScore: 0.0).Single(s => s.Type == enumType);
        Assert.That(site.Start, Is.EqualTo(start));
        var direct = MiRnaAnalyzer.ScoreTargetSiteContextPlusPlus(Utr, miRna, site,
            new MiRnaAnalyzer.ContextPlusPlusInputs(Sps: sps, Ta: 3.393, OrfLength: 1500, Orf8mCount: 1, PlfoldUnpairedProbability: plfold));
        var c = AnnotationTools.ScoreTargetSiteContextPlusPlus(Utr, Let7a, start, end, type,
            sps: sps, ta: 3.393, orfLength: 1500, orf8mCount: 1, plfoldUnpairedProbability: plfold).Score;
        Assert.Multiple(() =>
        {
            Assert.That(c.Intercept, Is.EqualTo(direct.Intercept));
            Assert.That(c.LocalAuContribution, Is.EqualTo(direct.LocalAuContribution));
            Assert.That(c.SRna1Contribution, Is.EqualTo(direct.SRna1Contribution));
            Assert.That(c.SRna8Contribution, Is.EqualTo(direct.SRna8Contribution));
            Assert.That(c.Site8Contribution, Is.EqualTo(direct.Site8Contribution));
            Assert.That(c.SaContribution, Is.EqualTo(direct.SaContribution));
            Assert.That(c.ThreePrimePairingContribution, Is.EqualTo(direct.ThreePrimePairingContribution));
            Assert.That(c.MinDistContribution, Is.EqualTo(direct.MinDistContribution));
            Assert.That(c.Len3UtrContribution, Is.EqualTo(direct.Len3UtrContribution));
            Assert.That(c.Off6mContribution, Is.EqualTo(direct.Off6mContribution));
            Assert.That(c.SpsContribution, Is.EqualTo(direct.SpsContribution));
            Assert.That(c.TaContribution, Is.EqualTo(direct.TaContribution));
            Assert.That(c.LenOrfContribution, Is.EqualTo(direct.LenOrfContribution));
            Assert.That(c.Orf8mContribution, Is.EqualTo(direct.Orf8mContribution));
            Assert.That(c.PctContribution, Is.EqualTo(direct.PctContribution));
            Assert.That(c.ContextScorePartial, Is.EqualTo(direct.ContextScorePartial));
            Assert.That(c.ContextScore, Is.EqualTo(direct.ContextScore));
            Assert.That(c.OmittedFeatures, Is.EqualTo(direct.OmittedFeatures));
        });
    }
}
