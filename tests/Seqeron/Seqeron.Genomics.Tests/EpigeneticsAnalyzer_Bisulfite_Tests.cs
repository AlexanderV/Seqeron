// EPIGEN-BISULF-001 — Bisulfite Sequencing Analysis
// Evidence: docs/Evidence/EPIGEN-BISULF-001-Evidence.md
// TestSpec: tests/TestSpecs/EPIGEN-BISULF-001.md
// Source: Frommer M et al. (1992). PNAS 89(5):1827–1831 (bisulfite conversion: unmeth C→U→T, 5mC stays C);
//         Krueger & Andrews (2011). Bioinformatics 27(11):1571–1572 (Bismark call: read C=meth, T=unmeth);
//         Bismark User Guide v0.15.0 (level = meth/(meth+unmeth));
//         Schultz MD et al. (2012). Trends Genet. 28(12):583–585 (weighted methylation level).

namespace Seqeron.Genomics.Tests;

using MethylationType = EpigeneticsAnalyzer.MethylationType;
using MethylationSite = EpigeneticsAnalyzer.MethylationSite;

[TestFixture]
public class EpigeneticsAnalyzer_Bisulfite_Tests
{
    #region SimulateBisulfiteConversion

    // M1 — Unmethylated cytosines convert to thymine (Frommer 1992: C→U, reads as T)
    [Test]
    public void Convert_UnmethylatedSequence_AllCytosinesBecomeT()
    {
        // ACGTCGAA: C@1 and C@4 both unmethylated → both become T.
        string result = EpigeneticsAnalyzer.SimulateBisulfiteConversion("ACGTCGAA");

        Assert.That(result, Is.EqualTo("ATGTTGAA"),
            "Unprotected cytosines convert to T; A/G/T unchanged (Frommer 1992)");
    }

    // M2 — Methylated (protected) cytosine remains C (Frommer 1992: 5mC nonreactive)
    [Test]
    public void Convert_MethylatedPositionProtected_StaysCytosine()
    {
        // ACGTCGAA with C@1 protected: C@1 stays C, C@4 unmethylated → T.
        string result = EpigeneticsAnalyzer.SimulateBisulfiteConversion(
            "ACGTCGAA", new HashSet<int> { 1 });

        Assert.That(result, Is.EqualTo("ACGTTGAA"),
            "5-methylcytosine is nonreactive and stays C; the other C converts (Frommer 1992)");
    }

    // M3 — Sequence with no cytosines is returned unchanged (only C reacts with bisulfite)
    [Test]
    public void Convert_NoCytosines_ReturnsUnchanged()
    {
        string result = EpigeneticsAnalyzer.SimulateBisulfiteConversion("AGTAGT");

        Assert.That(result, Is.EqualTo("AGTAGT"),
            "A, G, T are unaffected by bisulfite (Frommer 1992)");
    }

    // M4 — Lowercase unmethylated cytosine converts to lowercase t (case preserved)
    [Test]
    public void Convert_LowercaseUnmethylated_BecomesLowercaseT()
    {
        string result = EpigeneticsAnalyzer.SimulateBisulfiteConversion("acgt");

        Assert.That(result, Is.EqualTo("atgt"),
            "Lowercase c converts to lowercase t; case is preserved");
    }

    // S1 — Empty input yields empty output
    [Test]
    public void Convert_EmptyInput_ReturnsEmpty()
    {
        Assert.That(EpigeneticsAnalyzer.SimulateBisulfiteConversion(""), Is.EqualTo(""),
            "Empty sequence has nothing to convert");
    }

    // S2 — Null input yields empty output (documented validation contract)
    [Test]
    public void Convert_NullInput_ReturnsEmpty()
    {
        Assert.That(EpigeneticsAnalyzer.SimulateBisulfiteConversion(null!), Is.EqualTo(""),
            "Null sequence is treated as empty per the validation contract");
    }

    // C2 — Output length always equals input length (INV-01, property-based)
    [Test]
    public void Convert_AnyInput_PreservesLength()
    {
        string[] inputs = { "C", "CCCC", "ACGTACGTACGT", "GGGGTTTTAAAA", "cgCGcg" };

        Assert.Multiple(() =>
        {
            foreach (string input in inputs)
            {
                string result = EpigeneticsAnalyzer.SimulateBisulfiteConversion(input);
                Assert.That(result, Has.Length.EqualTo(input.Length),
                    $"Conversion is a per-base substitution, so length must be preserved (input '{input}')");
            }
        });
    }

    #endregion

    #region CalculateMethylationFromBisulfite

    // M5 — One C read and one T read at a CpG → level 0.5, coverage 2 (Bismark meth/(meth+unmeth))
    [Test]
    public void CalculateMethylation_HalfMethylatedCpG_ReturnsHalf()
    {
        // Reference ACGTACGT has CpG cytosines at index 1 and 5.
        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(
            "ACGTACGT",
            new[] { ("C", 1), ("T", 1) }).ToList();

        var site = sites.Single(s => s.Position == 1);
        Assert.Multiple(() =>
        {
            Assert.That(site.MethylationLevel, Is.EqualTo(0.5).Within(1e-10),
                "One methylated (C) and one unmethylated (T) call → 1/(1+1) = 0.5 (Bismark)");
            Assert.That(site.Coverage, Is.EqualTo(2),
                "Coverage is the number of valid C/T calls at the site");
            Assert.That(site.Type, Is.EqualTo(MethylationType.CpG),
                "Calling targets CpG context");
        });
    }

    // M6 — A single T read at a CpG → level 0.0 (Bismark: unmethylated)
    [Test]
    public void CalculateMethylation_AllUnmethylated_ReturnsZero()
    {
        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(
            "ACGTACGT",
            new[] { ("T", 5) }).ToList();

        var site = sites.Single(s => s.Position == 5);
        Assert.Multiple(() =>
        {
            Assert.That(site.MethylationLevel, Is.EqualTo(0.0).Within(1e-10),
                "T at a reference CpG C means unmethylated → 0/(0+1) = 0.0 (Bismark)");
            Assert.That(site.Coverage, Is.EqualTo(1),
                "One valid call → coverage 1");
        });
    }

    // M7 — A single C read at a CpG → level 1.0 (Bismark: methylated)
    [Test]
    public void CalculateMethylation_FullyMethylated_ReturnsOne()
    {
        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(
            "ACGTAA",
            new[] { ("C", 1) }).ToList();

        var site = sites.Single(s => s.Position == 1);
        Assert.Multiple(() =>
        {
            Assert.That(site.MethylationLevel, Is.EqualTo(1.0).Within(1e-10),
                "C at a reference CpG C means methylated → 1/(1+0) = 1.0 (Bismark)");
            Assert.That(site.Coverage, Is.EqualTo(1),
                "One valid call → coverage 1");
        });
    }

    // S3 — A CpG covered by no read is excluded (percentage undefined; INV-06)
    [Test]
    public void CalculateMethylation_NoReadsCoverSite_SiteExcluded()
    {
        // Cover only CpG@1; CpG@5 has no covering read.
        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(
            "ACGTACGT",
            new[] { ("C", 1) }).ToList();

        Assert.That(sites.Any(s => s.Position == 5), Is.False,
            "A CpG with zero coverage has undefined methylation % and is excluded (Bismark)");
    }

    // S5 — A read base that is neither C nor T at a CpG C is not a valid call and is ignored
    [Test]
    public void CalculateMethylation_NonCTReadBase_NotCounted()
    {
        // Read 'A' at CpG@1 is not a valid bisulfite call; a following T@1 is the only count.
        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(
            "ACGTAA",
            new[] { ("A", 1), ("T", 1) }).ToList();

        var site = sites.Single(s => s.Position == 1);
        Assert.Multiple(() =>
        {
            Assert.That(site.Coverage, Is.EqualTo(1),
                "Only the T call is valid; the A call is ignored (Krueger & Andrews 2011)");
            Assert.That(site.MethylationLevel, Is.EqualTo(0.0).Within(1e-10),
                "The single valid call (T) is unmethylated → 0.0");
        });
    }

    // C1 — Read extending past the reference end ignores the out-of-reference bases
    [Test]
    public void CalculateMethylation_ReadPastReferenceEnd_ExtraBasesIgnored()
    {
        // Reference length 6 (ACGTAA), CpG@1. Read 'CXX' from pos 1 covers CpG@1 with C;
        // bases beyond the reference are ignored.
        var sites = EpigeneticsAnalyzer.CalculateMethylationFromBisulfite(
            "ACGTAA",
            new[] { ("CGGGGG", 1) }).ToList();

        var site = sites.Single(s => s.Position == 1);
        Assert.That(site.Coverage, Is.EqualTo(1),
            "Only the in-reference C call counts; bases past the reference end are ignored");
    }

    #endregion

    #region GenerateMethylationProfile

    // M8 — Weighted per-context level under unequal coverage (Schultz 2012), distinct from unweighted mean
    [Test]
    public void GenerateProfile_UnequalCoverage_UsesWeightedLevel()
    {
        // CpG site A: level 1.0, coverage 10 (10 meth / 10 total)
        // CpG site B: level 0.0, coverage 30 (0 meth / 30 total)
        // Weighted = (1.0*10 + 0.0*30)/(10+30) = 0.25; unweighted mean would be 0.5.
        var sites = new[]
        {
            new MethylationSite(1, MethylationType.CpG, "CG", 1.0, 10),
            new MethylationSite(5, MethylationType.CpG, "CG", 0.0, 30),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.Multiple(() =>
        {
            Assert.That(profile.CpGMethylation, Is.EqualTo(0.25).Within(1e-10),
                "Weighted methylation = Σ(level·cov)/Σ(cov) = 10/40 = 0.25 (Schultz 2012), not the 0.5 unweighted mean");
            Assert.That(profile.TotalCpGSites, Is.EqualTo(2),
                "Both inputs are CpG sites");
        });
    }

    // M8b — Schultz et al. (2012) canonical worked example: a region with one site at
    // 90/100 methylated reads and one at 1/2. The paper contrasts the weighted methylation
    // level (read-pooled) with the unweighted mean of per-site levels.
    //   weighted = (90 + 1) / (100 + 2) = 91/102 = 0.8921568627...
    //   unweighted mean = (0.90 + 0.50) / 2 = 0.70  (NOT what the weighted level returns)
    // Source: Schultz MD, Schmitz RJ, Ecker JR (2012) Trends Genet. 28(12):583–585.
    [Test]
    public void GenerateProfile_SchultzWorkedExample_ReturnsWeightedLevel()
    {
        var sites = new[]
        {
            new MethylationSite(1, MethylationType.CpG, "CG", 0.90, 100),
            new MethylationSite(5, MethylationType.CpG, "CG", 0.50, 2),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.Multiple(() =>
        {
            Assert.That(profile.CpGMethylation, Is.EqualTo(91.0 / 102.0).Within(1e-12),
                "Weighted CpG methylation = (90+1)/(100+2) = 91/102 ≈ 0.89216 (Schultz 2012), not the 0.70 unweighted mean");
            Assert.That(profile.GlobalMethylation, Is.EqualTo(91.0 / 102.0).Within(1e-12),
                "All sites are CpG, so the global weighted level equals the CpG weighted level");
            Assert.That(profile.CpGMethylation, Is.Not.EqualTo(0.70).Within(1e-3),
                "Weighted level must differ from the unweighted mean 0.70 under unequal coverage (Schultz 2012)");
            Assert.That(profile.TotalCpGSites, Is.EqualTo(2),
                "Both inputs are CpG sites");
            Assert.That(profile.MethylatedCpGSites, Is.EqualTo(2),
                "Both per-site levels (0.90 and 0.50) are ≥ 0.5, so both count as methylated");
        });
    }

    // M9 — Per-context separation: CpG and CHG levels reported independently
    [Test]
    public void GenerateProfile_SeparatesContexts_ReportsPerContextLevels()
    {
        var sites = new[]
        {
            new MethylationSite(1, MethylationType.CpG, "CG", 1.0, 5),
            new MethylationSite(7, MethylationType.CHG, "CAG", 0.0, 5),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.Multiple(() =>
        {
            Assert.That(profile.CpGMethylation, Is.EqualTo(1.0).Within(1e-10),
                "The single fully methylated CpG site → CpG level 1.0");
            Assert.That(profile.CHGMethylation, Is.EqualTo(0.0).Within(1e-10),
                "The single unmethylated CHG site → CHG level 0.0");
        });
    }

    // M9b — CHH context is weighted independently (third Bismark context: Krueger 2011).
    // Two CHH sites: 2/2 methylated (level 1.0, cov 2) and 0/8 (level 0.0, cov 8).
    //   weighted CHH = (1.0*2 + 0.0*8)/(2+8) = 2/10 = 0.2 (Schultz read-pooled level)
    [Test]
    public void GenerateProfile_ChhContext_UsesWeightedLevel()
    {
        var sites = new[]
        {
            new MethylationSite(1, MethylationType.CHH, "CAA", 1.0, 2),
            new MethylationSite(9, MethylationType.CHH, "CTT", 0.0, 8),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.Multiple(() =>
        {
            Assert.That(profile.CHHMethylation, Is.EqualTo(0.2).Within(1e-12),
                "Weighted CHH = (1.0·2 + 0.0·8)/(2+8) = 2/10 = 0.2 (Schultz 2012)");
            Assert.That(profile.CpGMethylation, Is.EqualTo(0.0).Within(1e-12),
                "No CpG sites → CpG level 0");
            Assert.That(profile.TotalCpGSites, Is.EqualTo(0),
                "No CpG sites present");
        });
    }

    // M9c — A per-site level below the 0.5 methylated-site cutoff is not counted in
    // MethylatedCpGSites, but still contributes to the continuous weighted level.
    [Test]
    public void GenerateProfile_LevelBelowThreshold_NotCountedAsMethylatedSite()
    {
        var sites = new[]
        {
            new MethylationSite(1, MethylationType.CpG, "CG", 0.49, 100),
            new MethylationSite(5, MethylationType.CpG, "CG", 0.80, 100),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.Multiple(() =>
        {
            Assert.That(profile.MethylatedCpGSites, Is.EqualTo(1),
                "Only the 0.80 site is ≥ 0.5; the 0.49 site is below the cutoff");
            Assert.That(profile.TotalCpGSites, Is.EqualTo(2),
                "Both are CpG sites regardless of the methylated-site cutoff");
            Assert.That(profile.CpGMethylation, Is.EqualTo((0.49 + 0.80) / 2.0).Within(1e-12),
                "Equal coverage → weighted level equals the mean of per-site levels = 0.645");
        });
    }

    // S4 — Empty site list → all-zero profile
    [Test]
    public void GenerateProfile_NoSites_ReturnsZeroProfile()
    {
        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(Array.Empty<MethylationSite>());

        Assert.Multiple(() =>
        {
            Assert.That(profile.GlobalMethylation, Is.EqualTo(0.0).Within(1e-10),
                "No sites → zero global methylation");
            Assert.That(profile.TotalCpGSites, Is.EqualTo(0),
                "No sites → zero CpG count");
            Assert.That(profile.MethylationByPosition, Is.Empty,
                "No sites → empty per-position list");
        });
    }

    #endregion
}
