// EPIGEN-METHYL-001 — Methylation Analysis (context classification + profile)
// Evidence: docs/Evidence/EPIGEN-METHYL-001-Evidence.md
// TestSpec: tests/TestSpecs/EPIGEN-METHYL-001.md
// Source: Cornish-Bowden (1985) Nucl. Acids Res. 13:3021–3030 (IUPAC H = A,C,T);
//         Krueger & Andrews (2011) Bioinformatics 27(11):1571–1572 (Bismark contexts);
//         Lister et al. (2009) Nature 462:315–322; Schultz et al. (2012) Trends Genet. 28:583–585.

namespace Seqeron.Genomics.Tests;

using MethylationType = EpigeneticsAnalyzer.MethylationType;
using MethylationSite = EpigeneticsAnalyzer.MethylationSite;

[TestFixture]
public class EpigeneticsAnalyzer_Methylation_Tests
{
    #region GetMethylationContext

    // M1 — CpG: C immediately followed by G (Bismark; Lister 2009)
    [Test]
    public void GetMethylationContext_CpG_ReturnsCpG()
    {
        // A C G T : C at index 1, next base G → CpG
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("ACGT", 1);

        Assert.That(ctx, Is.EqualTo(MethylationType.CpG),
            "C followed by G is CpG per Bismark/Lister context definition");
    }

    // M2 — CHG: C, H(=A), G (H ∈ {A,C,T} per IUPAC)
    [Test]
    public void GetMethylationContext_CHG_ReturnsCHG()
    {
        // C A G : C at 0, H=A, then G → CHG
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("CAG", 0);

        Assert.That(ctx, Is.EqualTo(MethylationType.CHG),
            "C,H,G is CHG (H=A is a valid IUPAC H)");
    }

    // M3 — CHH: C, H(=A), H(=A)
    [Test]
    public void GetMethylationContext_CHH_ReturnsCHH()
    {
        // C A A : C at 0, H=A, H=A → CHH
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("CAA", 0);

        Assert.That(ctx, Is.EqualTo(MethylationType.CHH),
            "C,H,H is CHH (both H positions are A)");
    }

    // M4 — H excludes G: C followed by G is CpG, never CHG
    [Test]
    public void GetMethylationContext_NextBaseG_IsCpGNotCHG()
    {
        // C G G : next base is G → CpG (the G is not an H)
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("CGG", 0);

        Assert.That(ctx, Is.EqualTo(MethylationType.CpG),
            "Next base G makes this CpG; G is excluded from IUPAC H (Cornish-Bowden 1985)");
    }

    // M5 — CHG with H=C
    [Test]
    public void GetMethylationContext_CHGWithHC_ReturnsCHG()
    {
        // C C G : H=C, then G → CHG
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("CCG", 0);

        Assert.That(ctx, Is.EqualTo(MethylationType.CHG),
            "C,H=C,G is CHG (C is a valid IUPAC H)");
    }

    // M6 — CHG with H=T
    [Test]
    public void GetMethylationContext_CHGWithHT_ReturnsCHG()
    {
        // C T G : H=T, then G → CHG
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("CTG", 0);

        Assert.That(ctx, Is.EqualTo(MethylationType.CHG),
            "C,H=T,G is CHG (T is a valid IUPAC H)");
    }

    // M7 — Index not a cytosine → no context
    [Test]
    public void GetMethylationContext_NonCytosine_ReturnsNull()
    {
        // A at index 0 has no methylation context
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("ACGT", 0);

        Assert.That(ctx, Is.Null, "Only cytosines have a methylation context");
    }

    // M8 — Incomplete downstream context (non-CpG needs a third base)
    [Test]
    public void GetMethylationContext_IncompleteContext_ReturnsNull()
    {
        // C A : C at 0, next is H but no third base → cannot decide CHG vs CHH
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("CA", 0);

        Assert.That(ctx, Is.Null,
            "CHG/CHH need a third base; a truncated window is unclassified");
    }

    // M9 — Terminal CG is still an unambiguous CpG (only two bases needed)
    [Test]
    public void GetMethylationContext_TerminalCG_ReturnsCpG()
    {
        // A A C G : C at index 2, next G at the very end → CpG
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("AACG", 2);

        Assert.That(ctx, Is.EqualTo(MethylationType.CpG),
            "Terminal CG is CpG; CpG needs only the two bases C and G");
    }

    // S1 — Case-insensitive classification (INV-06)
    [Test]
    public void GetMethylationContext_Lowercase_ReturnsCpG()
    {
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("acgt", 1);

        Assert.That(ctx, Is.EqualTo(MethylationType.CpG),
            "Classification is case-insensitive (lowercase acgt → CpG at 1)");
    }

    // C1 — Non-ACGT base in context → undetermined
    [Test]
    public void GetMethylationContext_NonAcgtBase_ReturnsNull()
    {
        // C N G : N is not in {A,C,T}, so the H position is invalid → unclassified
        var ctx = EpigeneticsAnalyzer.GetMethylationContext("CNG", 0);

        Assert.That(ctx, Is.Null,
            "A non-ACGT base in the H position leaves the context undetermined");
    }

    // C2 — Index out of range → null
    [Test]
    public void GetMethylationContext_IndexOutOfRange_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(EpigeneticsAnalyzer.GetMethylationContext("ACGT", 10), Is.Null,
                "Index beyond length returns null");
            Assert.That(EpigeneticsAnalyzer.GetMethylationContext("ACGT", -1), Is.Null,
                "Negative index returns null");
            Assert.That(EpigeneticsAnalyzer.GetMethylationContext(null!, 0), Is.Null,
                "Null sequence returns null");
        });
    }

    #endregion

    #region FindMethylationSites

    // M10 — Enumerates all three contexts at exact 0-based positions
    // Evidence: CGACAGCAA → CpG@0 (CG), CHG@3 (CAG), CHH@6 (CAA)
    [Test]
    public void FindMethylationSites_AllContexts_ReturnsExactSites()
    {
        string sequence = "CGACAGCAA";
        // C G A C A G C A A
        // 0 1 2 3 4 5 6 7 8

        var sites = EpigeneticsAnalyzer.FindMethylationSites(sequence).ToList();

        var cpg = sites.Single(s => s.Type == MethylationType.CpG);
        var chg = sites.Single(s => s.Type == MethylationType.CHG);
        var chh = sites.Single(s => s.Type == MethylationType.CHH);

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(3),
                "CGACAGCAA has exactly 3 classifiable cytosines (one per context)");
            Assert.That(cpg.Position, Is.EqualTo(0), "CpG at position 0 (CG)");
            Assert.That(chg.Position, Is.EqualTo(3), "CHG at position 3 (C,A,G)");
            Assert.That(chh.Position, Is.EqualTo(6), "CHH at position 6 (C,A,A)");
        });
    }

    // M11 — Every reported position indexes a cytosine (INV-03)
    [Test]
    public void FindMethylationSites_PositionsAreZeroBasedAtCytosine()
    {
        string sequence = "CGACAGCAA";

        var sites = EpigeneticsAnalyzer.FindMethylationSites(sequence).ToList();

        Assert.That(sites.All(s => char.ToUpperInvariant(sequence[s.Position]) == 'C'),
            Is.True, "Every site position must index a cytosine (0-based)");
    }

    // M15 — Null and empty input → empty
    [Test]
    public void FindMethylationSites_NullOrEmpty_ReturnsEmpty()
    {
        Assert.Multiple(() =>
        {
            Assert.That(EpigeneticsAnalyzer.FindMethylationSites(null!).ToList(), Is.Empty,
                "Null sequence yields no sites");
            Assert.That(EpigeneticsAnalyzer.FindMethylationSites("").ToList(), Is.Empty,
                "Empty sequence yields no sites");
        });
    }

    // S2 — Lowercase input yields the same sites (INV-06)
    [Test]
    public void FindMethylationSites_Lowercase_ReturnsSameSites()
    {
        var sites = EpigeneticsAnalyzer.FindMethylationSites("cgacagcaa").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(sites, Has.Count.EqualTo(3), "Lowercase input yields the same 3 sites");
            Assert.That(sites.Single(s => s.Type == MethylationType.CpG).Position, Is.EqualTo(0));
            Assert.That(sites.Single(s => s.Type == MethylationType.CHG).Position, Is.EqualTo(3));
            Assert.That(sites.Single(s => s.Type == MethylationType.CHH).Position, Is.EqualTo(6));
        });
    }

    #endregion

    #region GenerateMethylationProfile

    // M12 — Weighted methylation level per Schultz (2012): Σmeth/Σtotal
    // CpG sites (0.8, cov 10) and (0.2, cov 10) → (8+2)/(10+10) = 0.5
    [Test]
    public void GenerateMethylationProfile_WeightedCpGLevel_ReturnsExact()
    {
        var sites = new[]
        {
            new MethylationSite(0, MethylationType.CpG, "CG", 0.8, 10),
            new MethylationSite(5, MethylationType.CpG, "CG", 0.2, 10),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.That(profile.CpGMethylation, Is.EqualTo(0.5).Within(1e-10),
            "Weighted CpG level = (0.8*10 + 0.2*10)/(10+10) = 10/20 = 0.5 (Schultz 2012)");
    }

    // S3 — Weighted level differs from unweighted mean under unequal coverage (INV-05)
    // (0.8, cov 90) and (0.2, cov 10) → (72+2)/(100) = 0.74; mean of fractions = 0.5
    [Test]
    public void GenerateMethylationProfile_UnequalCoverage_UsesWeightedLevel()
    {
        var sites = new[]
        {
            new MethylationSite(0, MethylationType.CpG, "CG", 0.8, 90),
            new MethylationSite(5, MethylationType.CpG, "CG", 0.2, 10),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.That(profile.CpGMethylation, Is.EqualTo(0.74).Within(1e-10),
            "Weighted level = (0.8*90 + 0.2*10)/100 = 74/100 = 0.74, not the 0.5 unweighted mean");
    }

    // M13 — CpG counts and the 0.5 descriptive cutoff
    [Test]
    public void GenerateMethylationProfile_CpGCounts_AreExact()
    {
        var sites = new[]
        {
            new MethylationSite(0, MethylationType.CpG, "CG", 0.8, 10), // ≥0.5 → methylated
            new MethylationSite(5, MethylationType.CpG, "CG", 0.2, 10), // <0.5 → not
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.Multiple(() =>
        {
            Assert.That(profile.TotalCpGSites, Is.EqualTo(2), "Two CpG sites present");
            Assert.That(profile.MethylatedCpGSites, Is.EqualTo(1),
                "One CpG site has level ≥ 0.5 (descriptive cutoff)");
        });
    }

    // S3b — Per-context levels are computed independently
    [Test]
    public void GenerateMethylationProfile_PerContext_AreSeparate()
    {
        var sites = new[]
        {
            new MethylationSite(0, MethylationType.CpG, "CG", 0.9, 10),
            new MethylationSite(3, MethylationType.CHG, "CAG", 0.5, 10),
            new MethylationSite(6, MethylationType.CHH, "CAA", 0.1, 10),
        };

        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(sites);

        Assert.Multiple(() =>
        {
            Assert.That(profile.CpGMethylation, Is.EqualTo(0.9).Within(1e-10), "CpG level = 0.9");
            Assert.That(profile.CHGMethylation, Is.EqualTo(0.5).Within(1e-10), "CHG level = 0.5");
            Assert.That(profile.CHHMethylation, Is.EqualTo(0.1).Within(1e-10), "CHH level = 0.1");
            // Global weighted level over all sites: (9+5+1)/30 = 0.5
            Assert.That(profile.GlobalMethylation, Is.EqualTo(0.5).Within(1e-10),
                "Global weighted level = (0.9+0.5+0.1)*10 / 30 = 15/30 = 0.5");
        });
    }

    // M14 — Empty site list → all-zero profile
    [Test]
    public void GenerateMethylationProfile_NoSites_ReturnsZeros()
    {
        var profile = EpigeneticsAnalyzer.GenerateMethylationProfile(Array.Empty<MethylationSite>());

        Assert.Multiple(() =>
        {
            Assert.That(profile.GlobalMethylation, Is.EqualTo(0.0), "Global = 0 for empty input");
            Assert.That(profile.CpGMethylation, Is.EqualTo(0.0), "CpG = 0");
            Assert.That(profile.CHGMethylation, Is.EqualTo(0.0), "CHG = 0");
            Assert.That(profile.CHHMethylation, Is.EqualTo(0.0), "CHH = 0");
            Assert.That(profile.TotalCpGSites, Is.EqualTo(0), "0 CpG sites");
            Assert.That(profile.MethylatedCpGSites, Is.EqualTo(0), "0 methylated CpG sites");
            Assert.That(profile.MethylationByPosition, Is.Empty, "No positions");
        });
    }

    #endregion
}
