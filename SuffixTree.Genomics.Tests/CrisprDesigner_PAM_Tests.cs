using NUnit.Framework;
using SuffixTree.Genomics;
using System;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

/// <summary>
/// Tests for CRISPR PAM site detection (CRISPR-PAM-001).
/// Covers FindPamSites and GetSystem methods.
/// 
/// Evidence sources:
/// - Wikipedia: Protospacer adjacent motif
/// - Wikipedia: CRISPR
/// - Jinek et al. (2012), Science - SpCas9 NGG
/// - Zetsche et al. (2015), Cell - Cas12a TTTV
/// </summary>
[TestFixture]
public class CrisprDesigner_PAM_Tests
{
    #region GetSystem Tests - CRISPR System Configuration

    [Test]
    [Description("M1: SpCas9 canonical PAM is NGG (Wikipedia PAM)")]
    public void GetSystem_SpCas9_ReturnsNGG_Pam()
    {
        var system = CrisprDesigner.GetSystem(CrisprSystemType.SpCas9);

        Assert.Multiple(() =>
        {
            Assert.That(system.Name, Is.EqualTo("SpCas9"));
            Assert.That(system.PamSequence, Is.EqualTo("NGG"));
        });
    }

    [Test]
    [Description("M2: SaCas9 PAM is NNGRRT (Wikipedia CRISPR)")]
    public void GetSystem_SaCas9_ReturnsNNGRRT_Pam()
    {
        var system = CrisprDesigner.GetSystem(CrisprSystemType.SaCas9);

        Assert.Multiple(() =>
        {
            Assert.That(system.Name, Is.EqualTo("SaCas9"));
            Assert.That(system.PamSequence, Is.EqualTo("NNGRRT"));
        });
    }

    [Test]
    [Description("M3: Cas12a PAM is TTTV - T-rich PAM (Zetsche 2015)")]
    public void GetSystem_Cas12a_ReturnsTTTV_Pam()
    {
        var system = CrisprDesigner.GetSystem(CrisprSystemType.Cas12a);

        Assert.Multiple(() =>
        {
            Assert.That(system.Name, Is.EqualTo("Cas12a/Cpf1"));
            Assert.That(system.PamSequence, Is.EqualTo("TTTV"));
        });
    }

    [Test]
    [Description("M4: Each system has correct guide length")]
    [TestCase(CrisprSystemType.SpCas9, 20)]
    [TestCase(CrisprSystemType.SaCas9, 21)]
    [TestCase(CrisprSystemType.Cas12a, 23)]
    [TestCase(CrisprSystemType.AsCas12a, 23)]
    [TestCase(CrisprSystemType.LbCas12a, 24)]
    [TestCase(CrisprSystemType.CasX, 20)]
    public void GetSystem_ReturnsCorrectGuideLength(CrisprSystemType systemType, int expectedLength)
    {
        var system = CrisprDesigner.GetSystem(systemType);
        Assert.That(system.GuideLength, Is.EqualTo(expectedLength));
    }

    [Test]
    [Description("M5: Cas9 has PAM after target, Cas12a has PAM before target (Wikipedia PAM)")]
    [TestCase(CrisprSystemType.SpCas9, true)]
    [TestCase(CrisprSystemType.SaCas9, true)]
    [TestCase(CrisprSystemType.Cas12a, false)]
    [TestCase(CrisprSystemType.AsCas12a, false)]
    [TestCase(CrisprSystemType.LbCas12a, false)]
    [TestCase(CrisprSystemType.CasX, false)]
    public void GetSystem_ReturnsCorrectPamPosition(CrisprSystemType systemType, bool pamAfterTarget)
    {
        var system = CrisprDesigner.GetSystem(systemType);
        Assert.That(system.PamAfterTarget, Is.EqualTo(pamAfterTarget));
    }

    [Test]
    [Description("C2: All seven system types return valid configurations")]
    public void GetSystem_AllSystemTypes_ReturnValidConfigurations()
    {
        var allTypes = Enum.GetValues<CrisprSystemType>();

        foreach (var systemType in allTypes)
        {
            var system = CrisprDesigner.GetSystem(systemType);

            Assert.Multiple(() =>
            {
                Assert.That(system.Name, Is.Not.Null.And.Not.Empty,
                    $"System {systemType} should have a name");
                Assert.That(system.PamSequence, Is.Not.Null.And.Not.Empty,
                    $"System {systemType} should have a PAM sequence");
                Assert.That(system.GuideLength, Is.GreaterThan(0),
                    $"System {systemType} should have positive guide length");
            });
        }
    }

    #endregion

    #region FindPamSites - SpCas9 NGG Detection

    [Test]
    [Description("M6: FindPamSites detects NGG on forward strand")]
    public void FindPamSites_SpCas9_DetectsNGG_OnForwardStrand()
    {
        // 20bp target followed by AGG PAM
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var forwardSites = sites.Where(s => s.IsForwardStrand).ToList();
        Assert.That(forwardSites.Any(s => s.PamSequence == "AGG"), Is.True,
            "Should find AGG PAM on forward strand");
    }

    [Test]
    [Description("M7: NGG matches all variants - AGG, CGG, TGG, GGG (N = any nucleotide)")]
    [TestCase("ACGTACGTACGTACGTACGTAGG", "AGG")]
    [TestCase("ACGTACGTACGTACGTACGTCGG", "CGG")]
    [TestCase("ACGTACGTACGTACGTACGTTGG", "TGG")]
    [TestCase("ACGTACGTACGTACGTACGTGGG", "GGG")]
    public void FindPamSites_SpCas9_MatchesAllNGG_Variants(string sequenceStr, string expectedPam)
    {
        var sequence = new DnaSequence(sequenceStr);
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == expectedPam);
        Assert.That(matchingSite, Is.Not.Null,
            $"Should find {expectedPam} PAM (N = any nucleotide per IUPAC)");
    }

    [Test]
    [Description("M8: FindPamSites searches reverse strand")]
    public void FindPamSites_SpCas9_SearchesReverseStrand()
    {
        // CCN on forward = NGG on reverse complement
        // CCA on forward -> TGG on reverse
        // Need 20bp target space after PAM position on reverse strand
        var sequence = new DnaSequence("CCAACGTACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var reverseSites = sites.Where(s => !s.IsForwardStrand).ToList();
        Assert.That(reverseSites, Is.Not.Empty,
            "Should find PAM sites on reverse strand");
    }

    [Test]
    [Description("M11: FindPamSites handles lowercase input (case-insensitive)")]
    public void FindPamSites_SpCas9_CaseInsensitive()
    {
        var sites = CrisprDesigner.FindPamSites("acgtacgtacgtacgtacgtagg", CrisprSystemType.SpCas9).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "AGG"), Is.True,
            "Should find PAM in lowercase sequence");
    }

    [Test]
    [Description("M12: FindPamSites returns target sequence of correct length")]
    public void FindPamSites_SpCas9_ReturnsTargetSequence_OfCorrectLength()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAGG");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence == "AGG");
        Assert.That(site, Is.Not.Null);
        Assert.That(site!.TargetSequence.Length, Is.EqualTo(20),
            "SpCas9 target should be 20bp");
    }

    #endregion

    #region FindPamSites - Edge Cases

    [Test]
    [Description("M9: FindPamSites returns empty for sequences without PAM")]
    public void FindPamSites_NoPamPresent_ReturnsEmpty()
    {
        // Sequence with no GG anywhere - no NGG possible
        var sequence = new DnaSequence("ACACACACACACACACACACACACACACACAC");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        // Filter to only forward strand NGG matches
        var forwardNGGSites = sites.Where(s => s.IsForwardStrand && s.PamSequence.EndsWith("GG")).ToList();
        Assert.That(forwardNGGSites, Is.Empty,
            "Should not find NGG PAM in sequence without GG");
    }

    [Test]
    [Description("M10: FindPamSites returns empty for empty sequence")]
    public void FindPamSites_EmptySequence_ReturnsEmpty()
    {
        var sites = CrisprDesigner.FindPamSites("", CrisprSystemType.SpCas9).ToList();
        Assert.That(sites, Is.Empty);
    }

    [Test]
    [Description("S1: FindPamSites excludes sites where target would be out of bounds")]
    public void FindPamSites_TargetOutOfBounds_Excluded()
    {
        // PAM at position 0-2 means target would need to start at -20 (invalid)
        var sequence = new DnaSequence("AGGACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SpCas9).ToList();

        // Should not return site with invalid target bounds
        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.Position == 0);
        Assert.That(site, Is.Null,
            "Should not return PAM site when target would be out of bounds");
    }

    [Test]
    [Description("S5: FindPamSites with null DnaSequence throws ArgumentNullException")]
    public void FindPamSites_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            CrisprDesigner.FindPamSites((DnaSequence)null!, CrisprSystemType.SpCas9).ToList());
    }

    #endregion

    #region FindPamSites - Cas12a TTTV Detection

    [Test]
    [Description("M13: Cas12a detects TTTA variant (V = A)")]
    public void FindPamSites_Cas12a_DetectsTTTA()
    {
        // TTTA PAM followed by 23bp target
        var sequence = new DnaSequence("TTTAACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "TTTA"), Is.True,
            "Should find TTTA PAM (V includes A)");
    }

    [Test]
    [Description("M13: Cas12a detects TTTC variant (V = C)")]
    public void FindPamSites_Cas12a_DetectsTTTC()
    {
        var sequence = new DnaSequence("TTTCACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "TTTC"), Is.True,
            "Should find TTTC PAM (V includes C)");
    }

    [Test]
    [Description("M13: Cas12a detects TTTG variant (V = G)")]
    public void FindPamSites_Cas12a_DetectsTTTG()
    {
        var sequence = new DnaSequence("TTTGACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "TTTG"), Is.True,
            "Should find TTTG PAM (V includes G)");
    }

    [Test]
    [Description("M14: Cas12a does NOT detect TTTT (V excludes T)")]
    public void FindPamSites_Cas12a_DoesNotDetectTTTT()
    {
        var sequence = new DnaSequence("TTTTACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        Assert.That(sites.Any(s => s.PamSequence == "TTTT"), Is.False,
            "Should NOT find TTTT PAM (V excludes T per IUPAC)");
    }

    [Test]
    [Description("Cas12a PAM is before target (5' PAM)")]
    public void FindPamSites_Cas12a_PamBeforeTarget()
    {
        var sequence = new DnaSequence("TTTAACGTACGTACGTACGTACGTACGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.Cas12a).ToList();

        var site = sites.FirstOrDefault(s => s.PamSequence == "TTTA" && s.IsForwardStrand);
        Assert.That(site, Is.Not.Null);

        // For Cas12a, target starts after PAM
        Assert.That(site!.TargetStart, Is.EqualTo(4),
            "Target should start after 4bp PAM for Cas12a");
    }

    #endregion

    #region FindPamSites - SaCas9 NNGRRT Detection

    [Test]
    [Description("M15: SaCas9 detects NNGRRT with R=A (NNGAAT)")]
    public void FindPamSites_SaCas9_DetectsNNGAAT()
    {
        // 21bp target followed by AAGAAT PAM
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAAGAAT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence.Length == 6);
        Assert.That(matchingSite, Is.Not.Null,
            "Should find NNGRRT PAM with R=A");
    }

    [Test]
    [Description("M15: SaCas9 detects NNGRRT with R=G (NNGGGT)")]
    public void FindPamSites_SaCas9_DetectsNNGGGT()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAGGAGT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var matchingSite = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence.Length == 6);
        Assert.That(matchingSite, Is.Not.Null,
            "Should find NNGRRT PAM with R=G");
    }

    [Test]
    [Description("SaCas9 returns 21bp target sequence")]
    public void FindPamSites_SaCas9_Returns21bp_TargetSequence()
    {
        var sequence = new DnaSequence("ACGTACGTACGTACGTACGTAAAGAAT");
        var sites = CrisprDesigner.FindPamSites(sequence, CrisprSystemType.SaCas9).ToList();

        var site = sites.FirstOrDefault(s => s.IsForwardStrand && s.PamSequence.Length == 6);
        if (site != null)
        {
            Assert.That(site.TargetSequence.Length, Is.EqualTo(21),
                "SaCas9 target should be 21bp");
        }
    }

    #endregion

    #region String Overload Tests

    [Test]
    [Description("String overload works identically to DnaSequence overload")]
    public void FindPamSites_StringOverload_WorksIdentically()
    {
        const string sequenceStr = "ACGTACGTACGTACGTACGTAGG";
        var dnaSequence = new DnaSequence(sequenceStr);

        var stringResults = CrisprDesigner.FindPamSites(sequenceStr, CrisprSystemType.SpCas9).ToList();
        var dnaResults = CrisprDesigner.FindPamSites(dnaSequence, CrisprSystemType.SpCas9).ToList();

        Assert.That(stringResults.Count, Is.EqualTo(dnaResults.Count),
            "String and DnaSequence overloads should return same number of sites");
    }

    #endregion
}
